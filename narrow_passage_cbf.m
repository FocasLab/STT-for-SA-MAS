clear; clc; close all

%% ── System Parameters ────────────────────────────────────────────────────
n_agents  = 2;          % number of agents
dt        = 0.02;       % time step
T         = 15;         % simulation time
steps     = T/dt;
Q_world   = 20;         % world size [0,Q]^2
Ds        = 2*13;        % safety distance
gamma_cbf = 0.2;        % CBF class-K gain  κ(h) = γ·h³
k1        = 0.5;        % nominal PD gain (position)
k2        = 1.8;        % nominal PD gain (velocity)
alpha_max = 4000.0;        % max acceleration (all agents equal)
beta_max  = 2000.0;        % max velocity

%% ── Agent initialisation (circle swap) ───────────────────────────────────
R_c   = 6.0;
cx    = Q_world/2;  cy = Q_world/2;
theta = linspace(0, 2*pi, n_agents+1)';  theta(end) = [];

pos0  = [150 60;
                150 322.5];
goals = [150 322.5;
        150 60];
r_obs=112;
%% ── Static obstacles (treated as agents with u=0, v=0) ───────────────────
% Each row: [cx, cy, radius]
obs_static = [0   191.25 r_obs;
             272 191.25 r_obs];        % one obstacle at centre
n_obs = size(obs_static,1);

%% ── State initialisation ─────────────────────────────────────────────────
pos = pos0;
vel = zeros(n_agents, 2);

pos_hist = zeros(steps+1, n_agents, 2);
vel_hist = zeros(steps+1, n_agents, 2);
u_hist   = zeros(steps,   n_agents, 2);

pos_hist(1,:,:) = pos;
vel_hist(1,:,:) = vel;

dist_goal_hist = zeros(steps+1, n_agents);
dist_goal_hist(1,:) = vecnorm(pos - goals, 2, 2)';

goal_tol = 0.10;
alpha_i  = alpha_max * ones(n_agents,1);

%% Snapshot settings
snap_times = [10,15];
snap_iters = round(snap_times / dt);
snap_data  = struct('pos',{},'vel',{},'time',{});

%% ════════════════════════════════════════════════════════════════════════
%%  SIMULATION LOOP
%% ════════════════════════════════════════════════════════════════════════
fprintf('Running CBF-QP simulation...\n');
sim_t = tic;

for k = 1:steps

    %% ── Nominal PD controller for each agent ─────────────────────────
    u_nom = zeros(n_agents, 2);
    for i = 1:n_agents
        if norm(pos(i,:) - goals(i,:)) > goal_tol
            u_nom(i,:) =  -k1*(pos(i,:) - goals(i,:)) - k2*vel(i,:);
            nm = norm(u_nom(i,:));
            if nm > alpha_i(i)
                u_nom(i,:) = u_nom(i,:);
            end
        end
    end

    %% ── Build centralized CBF-QP ─────────────────────────────────────
    Ndv   = 2*n_agents;
    H_qp  = 2*eye(Ndv);
    f_qp  = -2*u_nom(:)';

    A_ineq = [];
    b_ineq = [];

    %% Agent-Agent safety constraints
    for i = 1:n_agents
        for j = i+1:n_agents
            Dp  = pos(i,:) - pos(j,:);
            Dv  = vel(i,:) - vel(j,:);
            d   = norm(Dp);
            if d < 1e-6; continue; end

            ai  = alpha_i(i);
            aj  = alpha_i(j);
            sum_a = ai + aj;

            inner = dot(Dp,Dv);

            radicand = 2*sum_a*(d - Ds);
            if radicand <= 0
                radicand = 1e-4;
            end
            h_ij = sqrt(radicand) + inner/d;

            rhs = gamma_cbf * h_ij * d ...
                  - (inner)^2 / d^2  ...
                  + norm(Dv)^2       ...
                  + sum_a * inner / sqrt(radicand + 1e-9);

            row = zeros(1, Ndv);
            idx_i = (i-1)*2 + (1:2);
            idx_j = (j-1)*2 + (1:2);
            row(idx_i) = -Dp;
            row(idx_j) =  Dp;

            A_ineq = [A_ineq; row];
            b_ineq = [b_ineq; rhs];
        end
    end

    %% Agent-Obstacle safety constraints
    for i = 1:n_agents
        for o = 1:n_obs
            obs_pos  = obs_static(o,1:2);
            obs_rad  = obs_static(o,3);
            Ds_obs   = Ds + obs_rad;

            Dp  = pos(i,:) - obs_pos;
            Dv  = vel(i,:);
            d   = norm(Dp);
            if d < 1e-6; continue; end

            ai = alpha_i(i);
            inner = dot(Dp, Dv);

            radicand = 2*ai*(d - Ds_obs);
            if radicand <= 0; radicand = 1e-4; end

            h_io = sqrt(radicand) + inner/d;

            rhs = gamma_cbf * h_io * d ...
                  - inner^2/d^2          ...
                  + norm(Dv)^2           ...
                  + ai * inner / sqrt(radicand + 1e-9);

            row = zeros(1, Ndv);
            idx_i = (i-1)*2 + (1:2);
            row(idx_i) = -Dp;

            A_ineq = [A_ineq; row];
            b_ineq = [b_ineq; rhs];
        end
    end

    %% Solve QP
    opts = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
    u_nom_vec =  reshape(u_nom.',[],1);

    [u_sol, ~, exitflag] = quadprog(H_qp, -2*u_nom_vec, ...
                                    A_ineq, b_ineq, ...
                                    [], [], [], [], [], opts);
    if exitflag <= 0
        u_sol = zeros(Ndv,1);
        for i = 1:n_agents
            v_i = vel(i,:);
            nv  = norm(v_i);
            if nv > 1e-6
                idx = (i-1)*2+(1:2);
                u_sol(idx) = -alpha_i(i) * v_i'/nv;
            end
        end
    end

    u_mat = reshape(u_sol, 2, n_agents)';

    %% Integrate double integrator
    vel = vel + dt * u_mat;
    for i = 1:n_agents
        nv = norm(vel(i,:));
        if nv > beta_max
            vel(i,:) = vel(i,:);
        end
    end
    pos = pos + dt * vel;

    %% Store
    pos_hist(k+1,:,:) = pos;
    vel_hist(k+1,:,:) = vel;
    u_hist(k,:,:)     = u_mat;
    dist_goal_hist(k+1,:) = vecnorm(pos - goals, 2, 2)';

    %% Save snapshots
    for s = 1:length(snap_iters)
        if k == snap_iters(s)
            snap_data(s).pos  = pos;
            snap_data(s).vel  = vel;
            snap_data(s).time = k * dt;
        end
    end

    if all(vecnorm(pos - goals, 2, 2) < goal_tol)
        fprintf('All agents reached goal at t = %.2f s\n', k*dt);
        steps = k;
        % Fill any unsaved snapshots before breaking
        for s = 1:length(snap_iters)
            if numel(snap_data) < s || isempty(snap_data(s).time)
                snap_data(s).pos  = pos;
                snap_data(s).vel  = vel;
                snap_data(s).time = k * dt;
            end
        end
        break
    end
end

total_time = toc(sim_t);
fprintf('Wall-clock time: %.3f s\n', total_time);

% Trim
pos_hist       = pos_hist(1:steps+1,:,:);
vel_hist       = vel_hist(1:steps+1,:,:);
dist_goal_hist = dist_goal_hist(1:steps+1,:);
times          = (0:steps)*dt;

%% ════════════════════════════════════════════════════════════════════════
%%  SNAPSHOT PLOTS
%% ════════════════════════════════════════════════════════════════════════
theta_c    = linspace(0,2*pi,60);
r_agent    = 21;
r_goal_plt = 30;
snap_labels = {'Early Phase', 'Late Phase'};

figure(1)
set(gcf, 'Position', [100, 0, 1600, 1125]);
set(gcf, 'Color', 'w');

for s = 1:length(snap_times)
    subplot(1, 2, s);
    hold on; grid on; box on;

    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 18);
    axis([0 272 0 382.5]);

    % Draw static obstacles
    for o = 1:n_obs
        fill(obs_static(o,1) + obs_static(o,3)*cos(theta_c), ...
             obs_static(o,2) + obs_static(o,3)*sin(theta_c), ...
             'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
    end

    % Draw goals
    for i = 1:n_agents
        xc = goals(i,1); yc = goals(i,2);
        xg = xc + r_goal_plt*cos(theta_c);
        yg = yc + r_goal_plt*sin(theta_c);
        fill(xg, yg, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'g');
    end

    % Colors matching STT code
    col1 = [0.3010 0.7450 0.9330];   % agent 1: light blue
    col2 = [0.9290 0.6940 0.1250];   % agent 2: orange-yellow

    % Draw trajectories up to this snapshot
    n_pts = min(snap_iters(s)+1, size(pos_hist,1));
    xy1 = squeeze(pos_hist(1:n_pts, 1, :));
    xy2 = squeeze(pos_hist(1:n_pts, 2, :));
    plot(xy1(:,1), xy1(:,2), '-', 'Color', col1, 'LineWidth', 2.0);
    plot(xy2(:,1), xy2(:,2), '-', 'Color', col2, 'LineWidth', 2.0);

    % Draw agent circles and centre dots at snapshot position
    sd = snap_data(s);
    for i = 1:n_agents
        p = sd.pos(i,:);
        if i == 1
            plot(p(1) + r_agent*cos(theta_c), p(2) + r_agent*sin(theta_c), ...
                 '--', 'Color', col1, 'LineWidth', 2.5);
            plot(p(1), p(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        else
            plot(p(1) + r_agent*cos(theta_c), p(2) + r_agent*sin(theta_c), ...
                 '--', 'Color', col2, 'LineWidth', 2.5);
            plot(p(1), p(2), 'o', 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', 8);
        end
    end

    title(sprintf('t = %.1f s  (%s)', sd.time, snap_labels{s}), ...
          'Color', 'k', 'FontSize', 16);
    xlabel('X'); ylabel('Y');
    hold off;
end

%% ════════════════════════════════════════════════════════════════════════
%%  POST-SIMULATION PLOTS
%% ════════════════════════════════════════════════════════════════════════
colors = lines(n_agents);

%% Distance to goal
figure('Name','Distance to Goal');
hold on; grid on;
for i = 1:n_agents
    plot(times, dist_goal_hist(:,i), 'Color', colors(i,:), 'LineWidth',1.4);
end
yline(goal_tol,'k--','LineWidth',1.2,'DisplayName','Goal tolerance');
xlabel('Time (s)'); ylabel('Distance (m)');
title('Distance to Goal — All Agents');
legend(arrayfun(@(i)sprintf('Agent %d',i),1:n_agents,'UniformOutput',false),...
       'Location','northeast');

%% Nearest-neighbour distance (collision check)
dist_nn = compute_nn_dist(pos_hist);
figure('Name','Nearest-Neighbour Distance');
hold on; grid on;
for i = 1:n_agents
    plot(times, dist_nn(:,i), 'Color', colors(i,:), 'LineWidth',1.2);
end
yline(Ds,'k--','LineWidth',1.4,'DisplayName',sprintf('D_s = %.2f',Ds));
xlabel('Time (s)'); ylabel('Distance (m)');
title('Nearest-Neighbour Distance (must stay ≥ D_s)');
legend('Location','northeast');

%% Control effort
u_norm = zeros(steps, n_agents);
for k = 1:steps
    for i = 1:n_agents
        u_norm(k,i) = norm(squeeze(u_hist(k,i,:)));
    end
end
figure('Name','Control Effort');
plot(times(1:steps), u_norm, 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('||u_i||');
title('Control Input Magnitude');
legend(arrayfun(@(i)sprintf('Agent %d',i),1:n_agents,'UniformOutput',false));
grid on;

%% ── Helper ───────────────────────────────────────────────────────────────
function dist_nn = compute_nn_dist(pos_hist)
    T = size(pos_hist,1);
    n = size(pos_hist,2);
    dist_nn = zeros(T,n);
    for t = 1:T
        pt = squeeze(pos_hist(t,:,:));
        for i = 1:n
            d = vecnorm(pt - pt(i,:), 2, 2);
            d(i) = inf;
            dist_nn(t,i) = min(d);
        end
    end
end