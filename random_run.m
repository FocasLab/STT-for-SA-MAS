%% STT for SA-MAS with user defined Agents Position, Social Awareness and Obstacle Position
clc; clear all; close all;

%% USER CONFIGURATION
n_agents       = 2;   % select number of agnets
n_snaps        = 4;   % number of snapshot

%% Add Obstacle information in the form: [x, y, radius, vx, vy]
obstacles = [250 220 18 -4.4  0; ...
              70 180 18  2.0  0; ...
             175 310 18 -4.0  0; ...
              92  75 18  3.0  0];

%% Set STT Parameters
dT      = 30;% Target radius
h_1     = 0.05;% gain for target reaching
h_2     = 1;%gain for obstacle avoidance 
hat_h_2 = 8;%  %gain for inter agent avoidance 
rad_max = 0.9 * dT;% maximum allowable STT radius
rad_min = 0.7 * dT;%minimum allowable STT radius
b       = 10;% rate of decrease of SIF
dO      = 3;

%% STT center Initialisation
arena_x = 272; % arena for simulation
arena_y = 382.5;
if n_agents == 2
    cen_initial = [arena_x - dT,  dT; ... % tube initialization
                   dT,  arena_y - dT];
    eta        = [dT,  arena_y - dT; ... % Goal position
                   arena_x - dT,  dT];
else
    cen_initial = zeros(n_agents, 2);
    eta        = zeros(n_agents, 2);
    for k = 1:n_agents
        frac = (k - 1) / max(n_agents - 1, 1);
        cen_initial(k, :) = [dT,  dT + frac * (arena_y - 2*dT)];% tube initialization
        eta(k, :)        = [arena_x - dT,  arena_y - dT - frac * (arena_y - 2*dT)];% Goal position
    end
end

rad         = rad_max * ones(n_agents, 1);% radius of tube initialization
rad_initial = rad;
%% Agnet initialization
s_a = linspace(0.1, 0.9, n_agents);% set social awareness for each agents
robot       = cen_initial; % initializing system inside tube
tc   = 0.9 * 75; % prescribed time of convergence

%% Simulation Parameters
delt = 1e-2;
tf   = 75;
t    = 0:delt:tf;
control_agents = zeros(length(t), 2, n_agents);

%% Snapshot Time Configuration
snap_times = linspace(tf*0.05, tf*0.90, n_snaps);
snap_iter  = min(round(snap_times / delt) + 1, length(t));
snap_ready = false(n_snaps, 1);

%% Colour Maps
agent_colors = [0 1 1; 1 1 0; 1 0.5 0; 0.5 1 0; 0.8 0 0.8; 0 0.8 0.4; 1 0.3 0.3; 0.3 0.3 1];
robot_colors = [0 0.5 1; 1 0 1; 1 0.6 0; 0.4 0.9 0; 0.7 0 0.7; 0 0.7 0.3; 0.9 0.2 0.2; 0.2 0.2 0.9];
n_colors     = size(agent_colors, 1);
agent_colors = agent_colors(mod((0:n_agents-1), n_colors) + 1, :);
robot_colors = robot_colors(mod((0:n_agents-1), n_colors) + 1, :);

theta     = linspace(0, 2*pi, 600);
theta_obs = linspace(0, 2*pi, 400);

%% Figure: Snapshot Grid
n_cols = 2;
n_rows = 2;

fig2 = figure(2);
clf(fig2);
set(fig2, 'Position', [200, 80, 1100, 960], 'Color', 'k', 'Name', 'Trajectory Snapshots');

snap_axes = gobjects(n_snaps, 1);
for s = 1:n_snaps
    snap_axes(s) = subplot(n_rows, n_cols, s);
    set(snap_axes(s), 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 11);
    hold(snap_axes(s), 'on');
    box(snap_axes(s), 'on');
    axis(snap_axes(s), [0 arena_x 0 arena_y]);
    title(snap_axes(s), sprintf('t = %.1f s', snap_times(s)), 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
    for i = 1:n_agents
        fill(snap_axes(s), eta(i,1)+dT*cos(theta), eta(i,2)+dT*sin(theta), 'g', 'FaceAlpha', 0.18, 'EdgeColor', 'g', 'LineWidth', 1);
    end
end

%% Trail Storage
trail_tube  = cell(n_agents, 1);
trail_robot = cell(n_agents, 1);
trail_obs   = cell(size(obstacles,1), 1);

for i = 1:n_agents
    trail_tube{i}  = cen_initial(i,:);
    trail_robot{i} = robot(i,:);
end
for o = 1:size(obstacles,1)
    trail_obs{o} = obstacles(o,1:2);
end

%% STT Main Loop
tube           = cen_initial;
obstacle_tubes = cell(size(obstacles,1), 1);
disp('Simulation running...');
for iter = 1:length(t)
    current_time = t(iter);

    %% ── Update obstacle positions ────────────────────────────────────
    for o = 1:size(obstacles,1)
        obstacles(o,1) = obstacles(o,1) + delt * obstacles(o,4);
        obstacles(o,2) = obstacles(o,2) + delt * obstacles(o,5);

        if iter == 1
            obstacle_tubes{o} = [obstacles(o,1), obstacles(o,2), current_time];
        else
            obstacle_tubes{o} = [obstacle_tubes{o}; obstacles(o,1), obstacles(o,2), current_time];
        end
        trail_obs{o} = [trail_obs{o}; obstacles(o,1:2)];
    end

    %% ──  Compute first term of dcen (12) ─────────────────────────────────
    for k = 1:n_agents
        if t(iter) < tc
            gamma = -h_1 * tc * (cen_initial(k,:) - eta(k,:)) / (tc - t(iter));
        else
            gamma = [0, 0];
        end
        d_cen_2 = [0, 0]; % stores the second term of (12)
        d_cen_3 = [0, 0];  % stores the third term of (12)
        rad(k,1) = rad_max; % maximum radius
        sum_dist = 1e4;
        num_obs  = 0;
        %% second term computation for (12)
        for o = 1:size(obstacles,1)
            obs_pos    = obstacles(o,1:2);
            obs_radius = obstacles(o,3);
            d_1 = norm(cen_initial(k,:) - obs_pos) - obs_radius;
            d   = norm(cen_initial(k,:) - obs_pos) - (obs_radius + rad_min);
            if d_1 <= rad_max % Check for thoose obstacle inside sebsing range
                m_j = (1/d_1 - 1/rad_max) * h_2 * (1/d^2) * ((cen_initial(k,:)-obs_pos)/d);
                v_j = (mod(k,2) == 1) * ([0 -1;1 0]*m_j')' + (mod(k,2) == 0) * ([0  1;-1 0]*m_j')';
                d_cen_2  = d_cen_2 + m_j + v_j;
                num_obs  = num_obs + 1;
                sum_dist = smooth_min(sum_dist, d_1);
            end
        end
        if num_obs > 0 %% Update requird radius for obstacle avoidance
            rad(k,:) = smooth_min(rad_max, sum_dist);
        end
        % third term computation for (12)
        other_idx    = [1:k-1, k+1:n_agents];
        Ag_a_k       = [cen_initial(other_idx,:), rad_initial(other_idx,:)];
        s_a_values_k = s_a(other_idx);
        for l = 1:size(Ag_a_k,1)
            agent_pos = Ag_a_k(l,1:2);
            d_agent   = norm(cen_initial(k,:) - agent_pos) - 2*rad_min;
            phi       = phi_kl(s_a(k), s_a_values_k(l), t(iter), tc, b);
            if norm(cen_initial(k,:) - agent_pos) <= 2*rad_max%% consider only those agent inside the sensing radius
                hat_m_kl = (1/norm(cen_initial(k,:)-agent_pos) - 1/(2*rad_max)) * ...
                            hat_h_2 * (1/d_agent^2) * ((cen_initial(k,:)-agent_pos)/d_agent);
                hat_v_kl = (mod(k,2) == 1) * ([0 -1;1 0]*hat_m_kl')' + (mod(k,2) == 0) * ([0  1;-1 0]*hat_m_kl')';
                d_cen_3 = d_cen_3 + phi*(hat_m_kl + hat_v_kl);
            end
        end
        % Compute d_cen total (12) and then use euler update
        d_cen = gamma + d_cen_2 + d_cen_3;
        tube(k,:) = cen_initial(k,:) + delt * d_cen;
    end
    for k = 1:n_agents
        sum_dist     = rad(k);
        other_idx    = [1:k-1, k+1:n_agents];
        N_a_k_rad    = tube(other_idx,:);
        s_a_values_k = s_a(other_idx);
        
        for l = 1:size(N_a_k_rad,1)
            phi  = phi_kl(s_a(k), s_a_values_k(l), t(iter), tc, b);
            d_2k = norm(tube(k,:) - N_a_k_rad(l,:));
            if d_2k <= 2*rad_max
                w = (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi) + ...
                    (d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                sum_dist = smooth_min(sum_dist, rad_min+(d_2k-2*rad_min)*w);
            end
        end
        % update radius for inter agent collison avoidance
        rad(k) = smooth_min(rad_max, sum_dist);
        control_val = control(tube(k,:), rad(k), robot(k,:));
        control_agents(iter,:,k) = control_val;
        robot(k,:) = robot(k,:) + delt * control_val;
        trail_tube{k}  = [trail_tube{k};  tube(k,:)];
        trail_robot{k} = [trail_robot{k}; robot(k,:)];
    end

    %% ── Snapshot capture ─────────────────────────────────────────────
    for s = 1:n_snaps
        if ~snap_ready(s) && iter >= snap_iter(s)
            draw_snapshot_panel(snap_axes(s), s, snap_times, ...
                tube, robot, rad, obstacles, ...
                trail_tube, trail_robot, trail_obs, ...
                eta, agent_colors, robot_colors, ...
                n_agents, dT, dO, theta, theta_obs, arena_x, arena_y, s_a);
            snap_ready(s) = true;
        end
    end

    %% ── Termination ──────────────────────────────────────────────────
    if norm(tube - eta, 'fro') < 0.2
        fprintf('Goal reached at t = %.2f s\n', t(iter));
        break
    end
    cen_initial = tube;
    rad_initial = rad;
end
%% Final render for snapshots not yet captured
for s = 1:n_snaps
    if ~snap_ready(s)
        draw_snapshot_panel(snap_axes(s), s, snap_times, ...
            tube, robot, rad, obstacles, ...
            trail_tube, trail_robot, trail_obs, ...
            eta, agent_colors, robot_colors, ...
            n_agents, dT, dO, theta, theta_obs, arena_x, arena_y, s_a);
    end
end
drawnow;
disp('Simulation complete.');

%% Draw Snapshot Panel Function
function draw_snapshot_panel(ax, snap_idx, snap_times, ...
        tube_now, robot_now, rad_now, obs_now, ...
        trail_tube_, trail_robot_, trail_obs_, ...
        goal_, agent_colors_, robot_colors_, ...
        n_agents_, dT_, dO_, theta_, theta_obs_, arena_x_, arena_y_, s_a)

    cla(ax); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    axis(ax, [0 arena_x_ 0 arena_y_]);
    title(ax, sprintf('t = %.1f s', snap_times(snap_idx)), 'Color', 'w');

    for ii = 1:n_agents_
        fill(ax, goal_(ii,1)+dT_*cos(theta_), goal_(ii,2)+dT_*sin(theta_), 'g', 'FaceAlpha', 0.18, 'EdgeColor', 'g');
    end

    for oo = 1:size(obs_now,1)
        tr = trail_obs_{oo};
        plot(ax, tr(:,1), tr(:,2), '--', 'Color', [1 0.45 0.45 0.35]);
        fill(ax, obs_now(oo,1)+(obs_now(oo,3)-dO_)*cos(theta_obs_), obs_now(oo,2)+(obs_now(oo,3)-dO_)*sin(theta_obs_), 'r', 'FaceAlpha', 0.55);
    end

    for ii = 1:n_agents_
        ac = agent_colors_(ii,:); rc = robot_colors_(ii,:);
        tr_t = trail_tube_{ii};
        plot(ax, tr_t(:,1), tr_t(:,2), '-', 'Color', [ac, 0.6]);
        fill(ax, tube_now(ii,1)+rad_now(ii)*cos(theta_), tube_now(ii,2)+rad_now(ii)*sin(theta_), ac, 'FaceAlpha', 0.1, 'EdgeColor', ac);
        plot(ax, tube_now(ii,1), tube_now(ii,2), 'o', 'MarkerFaceColor', ac);
        plot(ax, robot_now(ii,1), robot_now(ii,2), 's', 'MarkerFaceColor', rc);
    end
end

%% Helper Functions
function phi = phi_kl(s_a_k, s_a_l, t, t_c_k, b)
    denom = s_a_l + s_a_k;
    if t < t_c_k
        phi = s_a_k / denom;
    else
        phi = (s_a_k / denom) * exp(-((t - t_c_k)^2) / (b^2));
    end
end

function y = smooth_min(x1, x2)
    mu = 1;
    y  = -(1/mu) * log(exp(-mu*x1) + exp(-mu*x2));
end

function u = control(tube_pos, rad_val, robot_pos)
    % Simple proportional control for the robot to follow the tube center
    K_p = 5;
    u = K_p * (tube_pos - robot_pos);
end