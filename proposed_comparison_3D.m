%% %% STT for SA-MAS for 3D drone case with velocity input
clc; clear; close all;


N_runs   = 25; % number of runs
n_a = 20;     % Number of agnets in SA-mAs

%% STT Parmeters
rad_max    = 1.2; % Maximum Allowable STT radius
rad_min  = 0.4;% Minimum Allowable STT radius
dT = 0.1; % target radius
h_1 = 0.15 * ones(1, n_a);  % gains for taget driven term in (12)
hat_h_2  = 0.04; % gains for agent avoidnace term in (12)
b =2; % rate at which SIF decrease 

%% Simulation Paramters
delt     = 2e-4;
tf       = 30;
tc_1     = 20;
t       = 0:delt:tf;
n_steps = length(t);
dist_sigma = 0;         % keep it non zero to run with disturbacne
b_max=30; % max size of cube where the tubes for each agent will be initialize
b_min = 25; % min size of cube where the tubes for each agent will be initialize
% choose social awareness value for each agents
s_a = [0.1; 0.9*ones(n_a-2,1); 0.1];


%%  RUN BOTH CONDITIONS
cond_labels = {'No Disturbance','With Disturbance'};

R_nd = struct();   % no-disturbance summary
R_wd = struct();   % with-disturbance summary

for cond = 1:2

    use_dist = (cond == 2);
    rng(42 + cond);

    %% Per-run storage
    success_flag        = false(N_runs, 1);
    collision_time      = nan(N_runs, 1);
    min_dist_runs       = zeros(N_runs, 1);
    avg_ctrl_time_runs  = zeros(N_runs, 1);
    std_ctrl_time_runs  = zeros(N_runs, 1);
    reach_time_runs     = nan(N_runs, n_a);
    total_sim_time_runs = zeros(N_runs, 1);
    pct_agents_reached  = zeros(N_runs, 1);
    bound_used          = zeros(N_runs, 1); 

    %% Monte-Carlo loop 
    for run = 1:N_runs
        %% Randomly draw box bound in [25, 30]
        bound = b_max + b_min*rand();
        bound_used(run) = bound;
        bounds = [0 bound; 0 bound; 0 bound];
        %% START OF BASE SIMULATION CODE (UNCHANGED)

        [cen_initial, eta] = generateDiagonal3D(n_a, bounds);
        robot = cen_initial; % initializing robots inside tube
        tc    = tc_1 * ones(1, n_a); % prescribed time of convergence for each agents
        rad   = rad_max * ones(n_a, 1);% initializing  radius of STT to maximum radius
        goal_time  = nan(n_a, 1);
        ctrl_times = zeros(n_steps, 1);
        min_dist   = inf;
        collision  = false;
        sim_t = tic;

        %% Main Simulation Loop
        for iter = 1:n_steps
            current_time = t(iter);
            t_ctrl = tic;
            %% -------- TUBE UPDATE --------
            for k = 1:n_a
                % first term calculation in (12)
                if current_time < tc(k)
                    gamma = -h_1(k)*tc(k)*(cen_initial(k,:)-eta(k,:)) / ...
                             (tc(k) - current_time);
                else
                    gamma = zeros(1,3);
                end
                d_cen_3 = zeros(1,3); % stores third term in (12)
                %% calculation of third term in (12)
                for l = 1:n_a
                    if l == k, continue; end
                    diff = cen_initial(k,:) - cen_initial(l,:);
                    dist_l = norm(diff);
                    d = dist_l - 2*rad_min;
                    if dist_l <= 2*rad_max && d > 0 % consider only those agent inside senisng radius
                        phi = phi_kl(s_a(k), s_a(l), t(iter), tc(k), b);
                        dir = diff / dist_l;
                        hat_m = hat_h_2 * (1/d^2) * dir;
                        hat_v = null(hat_m);
                        if norm(hat_v) > 0
                            hat_v = hat_v ;
                        end
                        d_cen_3 = d_cen_3 + ...
                            phi * ...
                            (hat_m + 0.3*hat_v(:,1)');
                    end
                end

                % Calculation of total d_cen in (12)
                d_cen_toal = gamma + d_cen_3;
                % euler update to get updated center of STT
                if norm(d_cen_toal) > 0
                    cen(k,:) = cen_initial(k,:) + delt*d_cen_toal;
                else
                    cen(k,:) = cen_initial(k,:);
                end
                % radius update for STT (17)
                sum_rad = 10^4;
                Ag_k_l_rad = cen_initial([1:k-1, k+1:end], :); % all agent except k-th agents
                s_a_values_rad = s_a([1:k-1, k+1:end]);
                for l = 1:size(Ag_k_l_rad,1)
                    phi = phi_kl(s_a(k), s_a_values_rad(l), t(iter), tc(k), b);
                    d_2k = norm(cen(k,:)-Ag_k_l_rad(l,:));
                    if d_2k <= 2*rad_max % consider only those agents inside sensing radius ie. l in N_a_k
                        W= (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi)+(d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                        sum_rad=smooth_min(sum_rad,rad_min+(d_2k-2*rad_min)*W);
                    end
                end
                rad(k) = smooth_min(rad_max, sum_rad);
            end

            %% Control Implementation using theorem (11)
            for k = 1:n_a
                u = control(cen(k,:), rad(k), robot(k,:));

                %% Disturbance (controlled by use_dist flag)
                if use_dist
                    noise = dist_sigma * randn(1,3);
                else
                    noise = zeros(1,3);
                end

                robot(k,:) = robot(k,:) + delt*(u + noise);

                %% Tube violation check
                if norm(cen(k,:) - robot(k,:)) > rad(k)
                    if ~collision
                        collision = true;
                        collision_time(run) = current_time;
                    end
                end

                %% Pairwise distance check (against all other tube centres)
                for l = 1:n_a
                    if l == k, continue; end
                    d_pair = norm(cen(k,:) - cen(l,:));
                    if d_pair < min_dist
                        min_dist = d_pair;
                    end
                    if (~collision && d_pair < 2*rad_min) ||(norm(cen(k,:)-robot(k,:))>rad(k))
                        collision = true;
                        collision_time(run) = current_time;
                    end
                end

                %% Goal check
                if isnan(goal_time(k)) && ...
                        norm(robot(k,:) - eta(k,:)) < dT
                    goal_time(k) = current_time;
                end
            end

            ctrl_times(iter) = toc(t_ctrl);
            cen_initial = cen;

            if ~any(isnan(goal_time)), break; end

        end   

        %% END OF BASE SIMULATION CODE

        total_sim_time_runs(run) = toc(sim_t);
        success_flag(run)        = ~collision;
        min_dist_runs(run)       = min_dist;
        avg_ctrl_time_runs(run)  = mean(ctrl_times(ctrl_times > 0));
        std_ctrl_time_runs(run)  = std(ctrl_times(ctrl_times > 0));
        reach_time_runs(run,:)   = goal_time';
        pct_agents_reached(run)  = 100 * mean(~isnan(goal_time));

    end   % run loop

    %% Aggregate metrics 
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

    avg_min_dist    = mean(min_dist_runs);
    std_min_dist    = std(min_dist_runs);
    worst_min_dist  = min(min_dist_runs);

    avg_sim_time    = mean(total_sim_time_runs);
    std_sim_time    = std(total_sim_time_runs);
    avg_pct_reached = mean(pct_agents_reached);

    col_times_valid = collision_time(~isnan(collision_time));

    %%  Print report 
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  STT 3D METHOD  –  %-48s║\n', cond_labels{cond});
    fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Agents=%-3d  Runs=%-3d  dt=%.4f  tf=%-4.0f  r_def=%.2f          ║\n',...
            n_a, N_runs, delt, tf, rad_max);
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
            2*rad_min, worst_min_dist);
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
            'Run','Bnd','Status','CtrlT(ms)','ReachT(s)','MinDist','%Reach');
    fprintf('║  %-3s  %-6s  %-7s  %-10s  %-10s  %-8s  %-6s  ║\n',...
            '---','------','-------','----------','----------','-------','------');
    for r = 1:N_runs
        st = 'OK';
        if ~success_flag(r), st = 'FAIL'; end
        fprintf('║  %-3d  %-6.2f  %-7s  %-10.5f  %-10.3f  %-8.4f  %-6.1f  ║\n',...
            r, bound_used(r), st,...
            avg_ctrl_time_runs(r)*1e3,...
            per_run_avg_reach(r),...
            min_dist_runs(r),...
            pct_agents_reached(r));
    end
    fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

    %% Store for comparison table
    if cond == 2   % No Disturbance
        R_nd.success_pct     = success_pct;
        R_nd.avg_ctrl        = overall_avg_ctrl;
        R_nd.std_ctrl        = overall_std_ctrl;
        R_nd.avg_min_dist    = avg_min_dist;
        R_nd.worst_min_dist  = worst_min_dist;
        R_nd.avg_reach       = avg_reach;
        R_nd.std_reach       = std_reach;
        R_nd.avg_pct_reached = avg_pct_reached;
        R_nd.avg_sim_time    = avg_sim_time;
    else           % With Disturbance
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

%%  COMPARISON TABLE
fprintf('\n╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STT 3D  Na=%-3d  Runs=%-3d  dt=%.4f  tf=%-4.0f  r_def=%.2f          ║\n',...
            n_a, N_runs, delt, tf, rad_max);
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

%%  LOCAL FUNCTIONS  

function [startPos, goalPos] = generateDiagonal3D(n_agents, bounds)
% bounds = [xmin xmax; ymin ymax; zmin zmax]

xmin = bounds(1,1); xmax = bounds(1,2);
ymin = bounds(2,1); ymax = bounds(2,2);
zmin = bounds(3,1); zmax = bounds(3,2);

center = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];

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

if n_agents <= size(corners,1)
    startPos = corners(1:n_agents,:);
else
    extra = n_agents - size(corners,1);
    rand_pts = [xmin + (xmax-xmin)*rand(extra,1), ...
                ymin + (ymax-ymin)*rand(extra,1), ...
                zmin + (zmax-zmin)*rand(extra,1)];
    startPos = [corners; rand_pts];
end

goalPos = 2*center - startPos;

end

%% --- Function to calculate SIF---
function phi = phi_kl(s_a_k, s_a_l, t, t_c_k, b)
    denom = s_a_l + s_a_k;
    if t < t_c_k
        phi = s_a_k / denom;
    else
        phi = (s_a_k / denom) * exp(-((t - t_c_k)^2) / (b^2));
    end
end

%% Smooth Minimum Function
function y = smooth_min(x1, x2)
    mu = 60;
    y = -(1/mu) * log( exp(-mu*x1) + exp(-mu*x2) );

end