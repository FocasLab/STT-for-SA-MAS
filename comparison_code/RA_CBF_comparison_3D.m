% % Y. Lyu, W. Luo, and J. M. Dolan. Risk-aware safe control for decentralized multi-agent systems via dynamic responsibility allocation. In IEEE/RSJ International Conference on
% Implementation of RA-CBF (19) for 3D drone case with velocity input

clear; clc; close all

%%  USER SETTINGS
N_runs      = 25;
N_list      = [5 10 20];   % agent counts to sweep
dt          = 0.01;
T           = 40;

Rsafe       = 0.4;             
gamma_cbf   = 2.0;
k_gain      = 0.8;
goal_tol    = 0.80;

R_min       = 25;
R_max       = 30;

run_with_disturbance =1; % set zero to run without disturbance
dist_sigma  = 1.0;             

cond_labels = {'No Disturbance','With Disturbance'};
n_conds     = 1 + run_with_disturbance;

%% OUTER LOOP OVER AGENT COUNTS
for n_idx = 1:length(N_list)
    N     = N_list(n_idx);
    steps = T / dt;

    % Agent SVO personalities (theta) — auto-sized to N
    theta_base = [2 8 8 8 2 8 8 8];
    if N <= length(theta_base)
        theta = theta_base(1:N);
    else
        theta = [theta_base, 8*ones(1, N - length(theta_base))];
    end

    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
    fprintf('║   R-SVO CBF 3-D  ——  N = %-3d agents                             ║\n', N);
    fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Runs=%-3d  dt=%.3f  T=%-4.0f  Rsafe=%.2f  gamma=%.1f  k=%.2f      ║\n',...
            N_runs, dt, T, Rsafe, gamma_cbf, k_gain);
    fprintf('╚══════════════════════════════════════════════════════════════════╝\n');

    %% Storage for No-Dist vs With-Dist comparison
    R_nd = struct();
    R_wd = struct();

    %% ─── Condition loop (no dist / with dist) ────────────────────────
    for cond = 1:n_conds

        use_dist = (cond == 2);
        rng(42 + cond + n_idx*100);

        %% Per-run accumulators
        success_flag        = false(N_runs, 1);
        collision_time      = nan(N_runs, 1);
        min_dist_all_runs   = zeros(N_runs, 1);
        avg_ctrl_time_runs  = zeros(N_runs, 1);
        std_ctrl_time_runs  = zeros(N_runs, 1);
        reach_time_runs     = nan(N_runs, N);
        total_sim_time_runs = zeros(N_runs, 1);
        pct_agents_reached  = zeros(N_runs, 1);
        circle_rad_used     = zeros(N_runs, 1);

        %% ─── Monte-Carlo loop ────────────────────────────────────────
        for run = 1:N_runs

            R_circle = R_min + (R_max - R_min)*rand();
            circle_rad_used(run) = R_circle;
            bounds = [0 R_circle; 0 R_circle; 0 R_circle];   % cube like your example
            [s3, goal3] = generateDiagonal3D(N, bounds);
            x    = s3';       % 3 x N
            goal = goal3';  
            ctrl_time    = zeros(steps, 1);
            reach_time   = nan(N, 1);
            min_dist_run = inf;
            collision    = false;

            sim_t = tic;

            for t_step = 1:steps

                current_time = t_step * dt;
                u      = zeros(3, N);
                t_ctrl = tic;

                %% ── QP loop ──────────────────────────────────────────
                for i = 1:N
                    xi    = x(:, i);
                    u_nom = -k_gain * (xi - goal(:, i));
                    A = []; b = [];

                    for j = 1:N
                        if i == j; continue; end

                        xj   = x(:, j);
                        d    = xi - xj;
                        dist = norm(d);

                        %% Track minimum pairwise distance
                        if dist < min_dist_run
                            min_dist_run = dist;
                        end

                        %% CBF safety function
                        h_val = dist^2 - (2*Rsafe)^2;

                        %% 3-D deadlock avoidance (cross-product perp)
                        if dist < 2.5*Rsafe
                            perp = cross(d, [0; 0; 1]);
                            if norm(perp) < 1e-6
                                perp = cross(d, [0; 1; 0]);
                            end
                            perp  = perp / norm(perp);
                            u_nom = perp;
                        end

                        %% R-SVO responsibility weight
                        phi_i   = theta(i) / (theta(i) + theta(j)) * pi/2;
                        omega_i = cos(phi_i)^2;

                        %% CBF constraint row  (1 x 3)
                        A = [A; -2*(xi - xj)'];
                        b = [b;  omega_i * gamma_cbf * h_val];
                    end

                    opts = optimoptions('quadprog', ...
                        'Display',             'off', ...
                        'Algorithm',           'interior-point-convex', ...
                        'OptimalityTolerance', 1e-16, ...
                        'ConstraintTolerance', 1e-16, ...
                        'StepTolerance',       1e-16, ...
                        'MaxIterations',       20000);

                    [u(:,i), ~, ~] = quadprog(eye(3), -u_nom, A, b, ...
                                              [], [], [], [], [], opts);
                end

                ctrl_time(t_step) = toc(t_ctrl);

                %% Disturbance (3-D, per-axis)
                if use_dist
                    noise = dist_sigma * randn(3, N);
                else
                    noise = zeros(3, N);
                end

                %% Euler step  x_dot = u + dist
                x = x + (u + noise) * dt;

                %% Collision check (all pairs)
                for i = 1:N-1
                    for j = i+1:N
                        d_pair = norm(x(:,i) - x(:,j));
                        if d_pair < min_dist_run
                            min_dist_run = d_pair;
                        end
                        if ~collision && d_pair < 2*Rsafe
                            collision = true;
                            collision_time(run) = current_time;
                        end
                    end
                end

                %% Goal reach
                for i = 1:N
                    if isnan(reach_time(i)) && ...
                            norm(x(:,i) - goal(:,i)) < goal_tol
                        reach_time(i) = current_time;
                    end
                end

                if ~any(isnan(reach_time)); break; end

            end   % time-step loop

            total_sim_time_runs(run) = toc(sim_t);
            success_flag(run)        = ~collision;
            min_dist_all_runs(run)   = min_dist_run;
            avg_ctrl_time_runs(run)  = mean(ctrl_time(ctrl_time > 0));
            std_ctrl_time_runs(run)  = std(ctrl_time(ctrl_time > 0));
            reach_time_runs(run,:)   = reach_time;
            pct_agents_reached(run)  = 100 * mean(~isnan(reach_time));

        end   % run loop

        %% ─── Aggregate metrics ───────────────────────────────────────
        success_pct      = 100 * mean(success_flag);
        fail_pct         = 100 - success_pct;
        n_success        = sum(success_flag);
        n_fail           = sum(~success_flag);

        overall_avg_ctrl = mean(avg_ctrl_time_runs) * 1e3;
        overall_std_ctrl = std(avg_ctrl_time_runs)  * 1e3;
        overall_min_ctrl = min(avg_ctrl_time_runs)  * 1e3;
        overall_max_ctrl = max(avg_ctrl_time_runs)  * 1e3;

        all_reach_flat = reach_time_runs(:);
        reached_mask   = ~isnan(all_reach_flat);
        avg_reach      = mean(all_reach_flat(reached_mask));
        std_reach      = std(all_reach_flat(reached_mask));
        min_reach      = min(all_reach_flat(reached_mask));
        max_reach      = max(all_reach_flat(reached_mask));

        per_run_avg_reach = mean(reach_time_runs, 2, 'omitnan');

        avg_min_dist   = mean(min_dist_all_runs);
        std_min_dist   = std(min_dist_all_runs);
        worst_min_dist = min(min_dist_all_runs);

        avg_sim_time   = mean(total_sim_time_runs);
        std_sim_time   = std(total_sim_time_runs);
        avg_pct_reached = mean(pct_agents_reached);

        col_times_valid = collision_time(~isnan(collision_time));

        %% ─── Print report (identical structure to 2-D) ───────────────
        fprintf('\n');
        fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  R-SVO CBF 3-D  [N=%-3d]  –  %-38s║\n', N, cond_labels{cond});
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  Agents=%-3d  Runs=%-3d  dt=%.3f  T=%-4.0f  Rsafe=%.2f            ║\n',...
                N, N_runs, dt, T, Rsafe);
        if use_dist
        fprintf('║  Disturbance sigma = %-44.4f║\n', dist_sigma);
        end
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  ── SAFETY ─────────────────────────────────────────────────    ║\n');
        fprintf('║  Success rate (no collision)         : %5.1f%%  (%d/%d runs)      ║\n',...
                success_pct, n_success, N_runs);
        fprintf('║  Failure rate                        : %5.1f%%  (%d/%d runs)      ║\n',...
                fail_pct, n_fail, N_runs);
        if ~isempty(col_times_valid)
        fprintf('║  Avg time of first collision (s)     : %8.3f                  ║\n',...
                mean(col_times_valid));
        fprintf('║  Earliest collision (s)              : %8.3f                  ║\n',...
                min(col_times_valid));
        end
        fprintf('║  Min pairwise dist – avg (m)         : %8.4f                  ║\n', avg_min_dist);
        fprintf('║  Min pairwise dist – std (m)         : %8.4f                  ║\n', std_min_dist);
        fprintf('║  Worst-case min dist (m) [2r=%.3f]  : %8.4f                  ║\n',...
                2*Rsafe, worst_min_dist);
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  ── CONTROL COMPUTATION TIME ───────────────────────────────    ║\n');
        fprintf('║  Avg ctrl time/step  (ms)            : %8.5f ± %-8.5f      ║\n',...
                overall_avg_ctrl, overall_std_ctrl);
        fprintf('║  Min ctrl time/step  (ms)            : %8.5f                  ║\n', overall_min_ctrl);
        fprintf('║  Max ctrl time/step  (ms)            : %8.5f                  ║\n', overall_max_ctrl);
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  ── GOAL REACHING ──────────────────────────────────────────    ║\n');
        fprintf('║  Avg %% agents reached                : %5.1f%%                    ║\n', avg_pct_reached);
        fprintf('║  Avg goal reach time (s)             : %8.3f ± %-8.3f      ║\n', avg_reach, std_reach);
        fprintf('║  Min goal reach time (s)             : %8.3f                  ║\n', min_reach);
        fprintf('║  Max goal reach time (s)             : %8.3f                  ║\n', max_reach);
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  ── SIMULATION TIMING ──────────────────────────────────────    ║\n');
        fprintf('║  Avg wall-clock/run (s)              : %8.3f ± %-8.3f      ║\n', avg_sim_time, std_sim_time);
        fprintf('║  Total wall-clock   (s)              : %8.3f                  ║\n', sum(total_sim_time_runs));
        fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
        fprintf('║  ── PER-RUN BREAKDOWN ──────────────────────────────────────    ║\n');
        fprintf('║  %-3s  %-6s  %-7s  %-10s  %-10s  %-8s  %-6s  ║\n',...
                'Run','R(m)','Status','CtrlT(ms)','ReachT(s)','MinDist','%Reach');
        fprintf('║  %-3s  %-6s  %-7s  %-10s  %-10s  %-8s  %-6s  ║\n',...
                '---','------','-------','----------','----------','-------','------');
        for r = 1:N_runs
            st = 'OK';
            if ~success_flag(r); st = 'FAIL'; end
            fprintf('║  %-3d  %-6.2f  %-7s  %-10.5f  %-10.3f  %-8.4f  %-6.1f  ║\n',...
                r, circle_rad_used(r), st,...
                avg_ctrl_time_runs(r)*1e3,...
                per_run_avg_reach(r),...
                min_dist_all_runs(r),...
                pct_agents_reached(r));
        end
        fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

        %% Store for comparison table
        if cond == 1
            R_nd.success_pct     = success_pct;
            R_nd.avg_ctrl        = overall_avg_ctrl;
            R_nd.std_ctrl        = overall_std_ctrl;
            R_nd.avg_min_dist    = avg_min_dist;
            R_nd.worst_min_dist  = worst_min_dist;
            R_nd.avg_reach       = avg_reach;
            R_nd.std_reach       = std_reach;
            R_nd.avg_pct_reached = avg_pct_reached;
            R_nd.avg_sim_time    = avg_sim_time;
        else
            R_wd.success_pct     = success_pct;
            R_wd.avg_ctrl        = overall_avg_ctrl;
            R_wd.std_ctrl        = overall_std_ctrl;
            R_wd.avg_min_dist    = avg_min_dist;
            R_wd.worst_min_dist  = worst_min_dist;
            R_wd.avg_reach       = avg_reach;
            R_wd.std_reach       = std_reach;
            R_wd.avg_pct_reached = avg_pct_reached;
            R_wd.avg_sim_time    = avg_sim_time;
        end

    end   % condition loop

    %% ─── No Dist vs With Dist comparison table ───────────────────────
    if run_with_disturbance
        fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  R-SVO CBF 3-D [N=%-3d]  –  No Disturbance vs With Disturbance  ║\n', N);
        fprintf('╠══════════════════════════════╦═══════════════╦═══════════════════╣\n');
        fprintf('║  Metric                      ║  No Dist      ║  With Dist        ║\n');
        fprintf('╠══════════════════════════════╬═══════════════╬═══════════════════╣\n');
        fprintf('║  Success rate (%%)            ║  %9.1f    ║  %9.1f          ║\n', R_nd.success_pct,    R_wd.success_pct);
        fprintf('║  Avg min pairwise dist (m)   ║  %9.4f    ║  %9.4f          ║\n', R_nd.avg_min_dist,   R_wd.avg_min_dist);
        fprintf('║  Worst min dist (m)          ║  %9.4f    ║  %9.4f          ║\n', R_nd.worst_min_dist, R_wd.worst_min_dist);
        fprintf('║  Avg ctrl time (ms)          ║  %9.5f    ║  %9.5f          ║\n', R_nd.avg_ctrl,       R_wd.avg_ctrl);
        fprintf('║  Std ctrl time (ms)          ║  %9.5f    ║  %9.5f          ║\n', R_nd.std_ctrl,       R_wd.std_ctrl);
        fprintf('║  Avg goal reach time (s)     ║  %9.3f    ║  %9.3f          ║\n', R_nd.avg_reach,      R_wd.avg_reach);
        fprintf('║  Std goal reach time (s)     ║  %9.3f    ║  %9.3f          ║\n', R_nd.std_reach,      R_wd.std_reach);
        fprintf('║  Avg %% agents reached        ║  %9.1f    ║  %9.1f          ║\n', R_nd.avg_pct_reached,R_wd.avg_pct_reached);
        fprintf('║  Avg sim wall-clock (s)      ║  %9.3f    ║  %9.3f          ║\n', R_nd.avg_sim_time,   R_wd.avg_sim_time);
        fprintf('╚══════════════════════════════╩═══════════════╩═══════════════════╝\n\n');
    end

end   % N_list loop

%% ═══════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTION
%% ═══════════════════════════════════════════════════════════════════════

function [startPos, goalPos] = generateDiagonal3D(n_agents, bounds)
    % bounds = [xmin xmax; ymin ymax; zmin zmax]
    
    xmin = bounds(1,1); xmax = bounds(1,2);
    ymin = bounds(2,1); ymax = bounds(2,2);
    zmin = bounds(3,1); zmax = bounds(3,2);
    
    % Center of the box
    center = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];
    
    %% Generate structured grid points (corners + midpoints if needed)
    
    % All 8 cube corners
    corners = [
        xmin ymin zmin;
        xmax ymax zmax;
        xmin ymax zmin;
        xmax ymin zmax;
        xmin ymin zmax;
        xmax ymax zmin;
        xmax ymin zmin;
        xmin ymax zmax
    ];
    
    % If more agents needed → random inside box
    if n_agents <= size(corners,1)
        startPos = corners(1:n_agents,:);
    else
        extra = n_agents - size(corners,1);
        rand_pts = [xmin + (xmax-xmin)*rand(extra,1), ...
                    ymin + (ymax-ymin)*rand(extra,1), ...
                    zmin + (zmax-zmin)*rand(extra,1)];
        startPos = [corners; rand_pts];
    end
    
    %% Diagonal opposite goals
    goalPos = 2*center - startPos;

end