%% %% STT for SA-MAS for 2D Mobile robot for Position swapping case for Comparison Stdy
clc; clear;

N_runs       = 25; % number of runs
n_a     = 5; %number of agents in SA-MAS
circle_rad   = 25;


%% STT Paramters
rad_max  = 0.9; % Maximum Allowable STT radius
rad_min = 0.6;% Minimum Allowable STT radius
dT = 0.3; % target radius
h_1 = [0.22 0.15 0.13 0.13 0.12 0.12 0.12 0.12 0.12*ones(1,n_a-8)]; % gains for taget driven term in (12)
h_1 = h_1(1:n_a); 
hat_h_2  = 0.04;% gains for agent avoidnace term in (12)
b=1; % rate at which SIF decreases

%% Simulation Prameters
delt         = 2e-4;
tf           = 30;
tc_1         = 20; % prescribed time of convergence for each agents
t = 0:delt:tf;
n_steps = length(t);
dist_sigma   = 2;     % disturbance std-dev
%  choose social awarness value of each agent between (0,1)
s_a = [0.1 0.99*ones(1,3) 0.1 0.99*ones(1,n_a-5)]';% Social Awareness Index Vaue for each agent
s_a = s_a(1:n_a);


%%  RUN BOTH Nominal and Disturbance Case
cond_labels = {'No Disturbance','With Disturbance'};

for cond = 1:2

    use_dist = (cond == 2);
    rng(42 + cond);

    %% Per-run storage
    success_flag        = false(N_runs,1);
    collision_time      = nan(N_runs,1);
    min_dist_runs       = zeros(N_runs,1);
    avg_ctrl_time_runs  = zeros(N_runs,1);
    std_ctrl_time_runs  = zeros(N_runs,1);
    reach_time_runs     = nan(N_runs,n_a);
    total_sim_time_runs = zeros(N_runs,1);
    pct_agents_reached  = zeros(N_runs,1);
    circle_rad_used     = zeros(N_runs,1);

    for run = 1:N_runs

        %% Random circle radius for initialization
        R_circle = circle_rad + 5*rand();
        circle_rad_used(run) = R_circle;
        center   = [10 10];

        [cen_initial, eta] = generateAgentsPolygon(n_a, R_circle, center); % tube center initialization and setting the goal for each agent

        robot     = cen_initial; % initialising each agent inside tube
        cen      = cen_initial;
        rad       = rad_max * ones(n_a,1); % initializing tube radius to maximum radius
        rad_p     = rad;
        tc        = tc_1 * ones(1,n_a); % presecribed time
        obstacles = [75 75 0.6 -0.8 -0.8];   % [x y r vx vy]
        goal_time   = nan(n_a,1); % array to store goal reaching time
        ctrl_times  = zeros(n_steps,1); % array to store average control computation time
        min_dist    = inf;
        collision   = false;
        sim_t = tic;

        % Main Simulation LOOP
        for iter = 1:n_steps
            current_time = t(iter);
            t_ctrl = tic;
            %% Tube update
            for k = 1:n_a
                %%  calculation of first term in (12)
                if t(iter) < tc(k)
                    gamma = -h_1(k)*tc(k)*(cen_initial(k,:)-eta(k,:)) ...
                            / (tc(k)-t(iter));
                else
                    gamma = [0 0];
                end
                d_cen_agent = [0 0]; % stores third term from (12)
                sum_dist    = 1e4;
                other_index    = [1:k-1, k+1:n_a];
                Ag_k_l = [cen_initial(other_index,:), rad_p(other_index)]; % consder all other agent excep k
                s_a_values_k   = s_a(other_index);
                for l = 1:size(Ag_k_l,1)
                    agent_pos = Ag_k_l(l,1:2);
                    d_agent   = norm(cen_initial(k,:)-agent_pos) - 2*rad_min;
                    phi = phi_kl(s_a(k), s_a_values_k(l), t(iter), tc, b);
                    if norm(cen_initial(k,:)-agent_pos) <= 2*rad_max % consder only those agnet inside sensing radius ie. l in N_a_k
                        hat_m = (1/norm(cen_initial(k,:)-agent_pos) ...
                            - 1/(2*rad_max)) * hat_h_2*(1/d_agent^2) * ...
                            ((cen_initial(k,:)-agent_pos) / ...
                             norm(cen_initial(k,:)-agent_pos));
                        hat_v = ([0 1;-1 0]*hat_m')';
                        d_cen_agent = d_cen_agent + ...
                            phi * ...
                            (hat_m + hat_v);
                    end
                end
                % calculate d_cen (12)
                d_cen = gamma + d_cen_agent;
                % using euler update to center of tube
                if norm(d_cen) > 0
                    cen(k,:) = cen_initial(k,:) + delt*d_cen;
                else
                    cen(k,:) = cen_initial(k,:);
                end
            end

            %% Radius Update in (17) and Control Implementation as perT heorem 11
            for k = 1:n_a
                sum_rad = 10^4;
                other_index = [1:k-1, k+1:n_a];
                Ag_k_l_rad = cen(other_index,:);
                s_a_values_rad = s_a(other_index);
                for l = 1:size(Ag_k_l_rad,1)
                    phi = phi_kl(s_a(k), s_a_values_rad(l), t(iter), tc, b);
                    d_2k = norm(cen(k,:)-Ag_k_l_rad(l,:));
                    if d_2k <= 2*rad_max % consider only those agents inside sensing radius ie. l in N_a_k
                        W= (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi)+(d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                        sum_dist=min(sum_dist,rad_min+(d_2k-2*rad_min)*W);
                    end
                end
                rad(k) = min(rad_max, sum_rad);

                %% Control call
                u = control(cen(k,:), rad(k), robot(k,:));

                %% Disturbance
                if use_dist
                    u = u + dist_sigma * randn(1,2);
                end

                robot(k,:) = robot(k,:) + delt*u;

                %% Pairwise distance check (robot space)
                for l = 1:size(Ag_k_l_rad,1)
                    d_pair = norm(cen(k,:) -Ag_k_l_rad(l,:));
                    if d_pair < min_dist
                        min_dist = d_pair;
                    end
                     
                    if (~collision && d_pair < 2*rad_min) ||(norm(cen(k,:)-robot(k,:))>rad(k))
                        
                        disp("fail")
                        collision          = true;
                        collision_time(run)= current_time;
                        break;
                    end
                end

                %% Goal check
                if isnan(goal_time(k)) && ...
                        norm(robot(k,:)-eta(k,:)) < dT
                    goal_time(k) = current_time;
                end
            end

            ctrl_times(iter) = toc(t_ctrl);
            cen_initial = cen;
            rad_p        = rad;

            if ~any(isnan(goal_time)); break; end
        end

        total_sim_time_runs(run) = toc(sim_t);
        success_flag(run)        = ~collision;
        min_dist_runs(run)       = min_dist;
        avg_ctrl_time_runs(run)  = mean(ctrl_times(ctrl_times>0));
        std_ctrl_time_runs(run)  = std(ctrl_times(ctrl_times>0));
        reach_time_runs(run,:)   = goal_time';
        pct_agents_reached(run)  = 100*mean(~isnan(goal_time));

    end  % run loop

    %% Aggregate metrics 
    success_pct      = 100*mean(success_flag);
    fail_pct         = 100 - success_pct;
    n_success        = sum(success_flag);
    n_fail           = sum(~success_flag);

    overall_avg_ctrl = mean(avg_ctrl_time_runs)*1e3;
    overall_std_ctrl = std(avg_ctrl_time_runs)*1e3;
    overall_min_ctrl = min(avg_ctrl_time_runs)*1e3;
    overall_max_ctrl = max(avg_ctrl_time_runs)*1e3;

    all_reach_flat  = reach_time_runs(:);
    reached_mask    = ~isnan(all_reach_flat);
    avg_reach       = mean(all_reach_flat(reached_mask));
    std_reach       = std(all_reach_flat(reached_mask));
    min_reach       = min(all_reach_flat(reached_mask));
    max_reach       = max(all_reach_flat(reached_mask));

    per_run_avg_reach = nanmean(reach_time_runs,2);

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
    fprintf('║  STT METHOD  –  %-50s║\n', cond_labels{cond});
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
            min_dist_runs(r),...
            pct_agents_reached(r));
    end
    fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

    %% Store for comparison
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

end  

%%  COMPARISON TABLE
fprintf('\n╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STT METHOD  –  COMPARISON: No Disturbance vs With Disturbance  ║\n');
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
function [startPos, goalPos] = generateAgentsPolygon(n, R, center)
    if nargin < 2; R = 12; end
    if nargin < 3; center = [0 0]; end
    if n > 100 || n < 3; error('n must be between 3 and 100'); end
    theta = linspace(0,2*pi,n+1); theta(end) = [];
    startPos = [center(1)+R*cos(theta(:)), center(2)+R*sin(theta(:))];
    goalPos  = [center(1)+R*cos(theta(:)+pi), center(2)+R*sin(theta(:)+pi)];
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