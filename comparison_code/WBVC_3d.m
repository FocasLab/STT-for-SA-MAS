%% Weighted Buffered Voronoi Cells (WBVC) – 3-D Extension
%  Dynamics : ẋ = u  (single-integrator, 3-D)
%  Scenario  : n-agent diagonal swap via generateDiagonal3D
%
%  Key changes vs. 2-D version
%    • All positions / velocities / cells live in R³
%    • Halfplane test uses 3-D dot products
%    • WBVC polyhedron visualised via QP boundary sampling + convhull
%    • Closest point uses quadprog (QP projection onto H-rep polyhedron)
%    • No external toolboxes required (uses only Optimization Toolbox)

clc; clear; close all;

%% ── Parameters ────────────────────────────────────────────────────────────
n_agents = 10;           % agents (≤8 use cube corners; >8 adds random pts)

% Workspace bounding box  [xmin xmax; ymin ymax; zmin zmax]
bounds = [0 20; 0 20; 0 20];
Q      = bounds(:,2);   % upper corner (lower corner is 0 by construction)

r      = 0.40;          % agent radius (m)
k_gain = 0.6;           % proportional gain
dt     = 0.04;          % time-step (s)
eps_dead = 0.5;        % deadlock threshold (m)
max_steps = 2000;       % max iterations
goal_tol  = 0.10;       % goal tolerance (m)

%% ── Initial & goal positions ──────────────────────────────────────────────
[positions, goals] = generateDiagonal3D(n_agents, bounds);
goals = goals+0.1*rand(n_agents, 1);
radii  = r * ones(n_agents, 1);

% SVO values γ ∈ [0,1]: 0 = altruistic, 1 = egoistic
% Assign thirds: egoistic / prosocial / altruistic
third  = floor(n_agents / 3);
gammas = [ones(1,third), 0.5*ones(1,third), zeros(1, n_agents - 2*third)];
gammas = gammas(randperm(n_agents));   % shuffle

fprintf('\n======  3-D WBVC: %d-Agent Diagonal Swap  ======\n', n_agents);
fprintf('Bounds : [0 %.0f] × [0 %.0f] × [0 %.0f]\n', Q(1), Q(2), Q(3));
fprintf('γ values: '); fprintf('%.2f  ', gammas); fprintf('\n');

%% ── Run simulation ────────────────────────────────────────────────────────
[history, times, dist_goal, ctrl_time_accum, ctrl_calls, ...
 goal_reach_iter, total_sim_time] = ...
    run_wbvc_3d(positions, goals, gammas, radii, bounds, ...
                k_gain, dt, eps_dead, max_steps, goal_tol);

print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
             goal_reach_iter, times, total_sim_time, '3-D Diagonal Swap');

%% ── Visualise ─────────────────────────────────────────────────────────────
T      = size(history,1);
snaps  = unique([1, max(1,round(T/3)), max(1,round(2*T/3)), T]);
labels = arrayfun(@(s) sprintf('t = %.1f s', times(s)), snaps, 'UniformOutput',false);

fig1 = figure('Name','3-D WBVC – Snapshots','NumberTitle','off');
fig1.Position = [60 60 1400 900];

for k = 1:min(4, numel(snaps))
    ax = subplot(2,2,k);
    pos_snap = squeeze(history(snaps(k),:,:));
    setup_axes_3d(ax, bounds, labels{k});
    draw_wbvc_cells_3d(ax, pos_snap, radii, gammas, bounds);
    draw_agents_3d(ax, pos_snap, radii, gammas, goals);
    draw_trajectories_3d(ax, history, gammas, snaps(k));
    if k == 4; add_legend(ax); end
end
sgtitle(fig1, sprintf('3-D WBVC: %d-Agent Diagonal Swap – Snapshots', n_agents), ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ── Distance plots ─────────────────────────────────────────────────────────
dist_nn = compute_nn_dist_hist(history);

fig2 = figure('Name','3-D WBVC – Distances','NumberTitle','off');
fig2.Position = [80 80 1100 420];
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
for i = 1:n_agents
    c = [svo_color(gammas(i)), 0.85];
    plot(ax1, times, dist_goal(:,i), 'Color', c, 'LineWidth', 1.0);
    plot(ax2, times, dist_nn(:,i),   'Color', c, 'LineWidth', 1.0);
end
yline(ax2, 2*r, 'k--', 'LineWidth', 1.2, 'DisplayName', '2r (collision)');
xlabel(ax1,'Time (s)'); ylabel(ax1,'Distance (m)');
title(ax1,'Distance to Goal'); ax1.FontSize = 11;
xlabel(ax2,'Time (s)'); ylabel(ax2,'Distance (m)');
title(ax2,'Nearest-Neighbour Distance');
legend(ax2,'FontSize',9); ax2.FontSize = 11;
add_legend(ax1);
sgtitle(fig2,'3-D WBVC: Distance Metrics','FontSize',13,'FontWeight','bold');


%% ══════════════════════════════════════════════════════════════════════════
%%  POSITION GENERATOR
%% ══════════════════════════════════════════════════════════════════════════
function [startPos, goalPos] = generateDiagonal3D(n_agents, bounds)
    % bounds = [xmin xmax; ymin ymax; zmin zmax]
    xmin = bounds(1,1); xmax = bounds(1,2);
    ymin = bounds(2,1); ymax = bounds(2,2);
    zmin = bounds(3,1); zmax = bounds(3,2);

    center = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];

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

    if n_agents <= size(corners,1)
        startPos = corners(1:n_agents,:);
    else
        extra    = n_agents - size(corners,1);
        rand_pts = [xmin + (xmax-xmin)*rand(extra,1), ...
                    ymin + (ymax-ymin)*rand(extra,1), ...
                    zmin + (zmax-zmin)*rand(extra,1)];
        startPos = [corners; rand_pts];
    end

    % Diagonal-opposite goals
    goalPos = 2*center - startPos;
end


%% ══════════════════════════════════════════════════════════════════════════
%%  CORE 3-D WBVC RUNNER
%% ══════════════════════════════════════════════════════════════════════════
function [history, times, dist_goal_hist, ctrl_time_accum, ctrl_calls, ...
          goal_reach_iter, total_sim_time] = ...
    run_wbvc_3d(pos0, goals, gammas, radii, bounds, ...
                k_gain, dt, eps_dead, max_steps, goal_tol)

    n   = size(pos0,1);
    pos = pos0;

    history            = zeros(max_steps+1, n, 3);
    history(1,:,:)     = pos;
    times              = zeros(max_steps+1, 1);
    dist_goal_hist     = zeros(max_steps+1, n);
    dist_goal_hist(1,:)= vecnorm(pos - goals, 2, 2)';

    ctrl_time_accum = zeros(max_steps,1);
    ctrl_calls      = zeros(max_steps,1);
    goal_reach_iter = -ones(n,1);

    sim_start = tic;

    step = 0;
    for iter = 1:max_steps
        new_pos = pos;
        t_ctrl  = tic;

        for i = 1:n
            if norm(pos(i,:) - goals(i,:)) < goal_tol
                continue
            end

            %% ẋ = u  (single integrator)
            ui    = k_gain * (goals(i,:) - pos(i,:));
            p_bar = pos(i,:) + ui * dt;           % predicted next position

            % if point_in_wbvc_3d(p_bar, i, pos, radii, gammas)
            %     p_target = p_bar;
            % else
            %     p_hat = closest_point_in_cell_3d(goals(i,:), i, pos, ...
            %                                       radii, gammas, bounds);
            %     if norm(pos(i,:) - p_hat) < eps_dead
            %         p_hat = deadlock_resolve_3d(p_hat, i, pos, radii, ...
            %                                     gammas, bounds);
            %     end
            %     ui2      = k_gain * (p_hat - pos(i,:));
            %     p_target = pos(i,:) + ui2 * dt;
            % end
                if point_in_wbvc_3d(p_bar, i, pos, radii, gammas)
        
                    % Safe to move normally
                    p_target = p_bar;
                
                else
                    
                    % Project goal into WBVC
                    p_hat = closest_point_in_cell_3d(goals(i,:), i, pos, ...
                                                      radii, gammas, bounds);
                
                    % Deadlock handling
                    % if norm(pos(i,:) - p_hat) < eps_dead
                    %     p_hat = deadlock_resolve_3d(p_hat, i, pos, radii, ...
                    %                                 gammas, bounds);
                    % end
                
                    % CRITICAL FIX: move ONLY up to p_hat (no overshoot)
                    step_vec = p_hat - pos(i,:);
                    step_mag = norm(step_vec);
                
                    max_step = k_gain * dt;
                
                    if step_mag > max_step
                        step_vec = (step_vec / step_mag) * max_step;
                    end
                
                    p_target = pos(i,:) + step_vec;
                end

            % Clip to workspace
            % p_target = max(p_target, bounds(:,1)' + radii(i));
            % p_target = min(p_target, bounds(:,2)' - radii(i));

            new_pos(i,:) = p_target;

            if goal_reach_iter(i) == -1 && ...
                    norm(new_pos(i,:) - goals(i,:)) < goal_tol
                goal_reach_iter(i) = iter;
            end


            for k=1:n-1
                    for j=k+1:n
                        d = norm(pos(k,:) - pos(j,:));
                        if d < 2*0.4-2*0.01
                            step
                            d
                            disp("collided")
                            break
                        end
                    end
                end
        end

        ctrl_time_accum(iter) = toc(t_ctrl);
        ctrl_calls(iter)      = 1;
        pos    = new_pos;
        step   = step + 1;
        history(step+1,:,:)      = pos;
        times(step+1)            = iter * dt;
        dist_goal_hist(step+1,:) = vecnorm(pos - goals, 2, 2)';

        if all(vecnorm(pos - goals, 2, 2) < goal_tol)
            break
        end
    end

    total_sim_time = toc(sim_start);

    history        = history(1:step+1,:,:);
    times          = times(1:step+1);
    dist_goal_hist = dist_goal_hist(1:step+1,:);
end


%% ══════════════════════════════════════════════════════════════════════════
%%  3-D WBVC MATHS
%% ══════════════════════════════════════════════════════════════════════════

function w = wij_weight(pi, pj, ri, rj, gi, gj)
    dij = norm(pi - pj);
    if dij < 1e-12; w = 0; return; end
    dg  = gi - gj;
    w   = (1 + 0.5*dg)*(ri + rj)*dij - 0.5*dg*dij^2;
end

function inside = point_in_halfplane_3d(q, pi, pj, ri, rj, gi, gj)
    w      = wij_weight(pi, pj, ri, rj, gi, gj);
    inside = (dot(pi-q, pi-q) <= dot(pj-q, pj-q) - w + 1e-9);
end

function inside = point_in_wbvc_3d(q, i, pos, radii, gammas)
    n = size(pos,1);
    for j = 1:n
        if j == i; continue; end
        if ~point_in_halfplane_3d(q, pos(i,:), pos(j,:), ...
                                   radii(i), radii(j), gammas(i), gammas(j))
            inside = false; return;
        end
    end
    inside = true;
end

function [A_ineq, b_ineq] = wbvc_halfplanes(i, pos, radii, gammas, bounds)
    % Returns H-rep: A_ineq * x <= b_ineq  (includes box constraints)
    n  = size(pos,1);
    nd = 3;

    % Box constraints: +I*x <= Q,  -I*x <= 0
    A_box = [eye(nd); -eye(nd)];
    b_box = [bounds(:,2); -bounds(:,1)];

    A_wbvc = zeros(n-1, nd);
    b_wbvc = zeros(n-1, 1);
    row = 0;
    for j = 1:n
        if j == i; continue; end
        pj = pos(j,:); pi_ = pos(i,:);
        dij = norm(pi_ - pj);
        if dij < 1e-12; continue; end
        w     = wij_weight(pi_, pj, radii(i), radii(j), gammas(i), gammas(j));
        % Halfplane: 2*(pj-pi)'*x <= ||pj||² - ||pi||² - w
        a_row = 2*(pj - pi_);
        b_row = dot(pj,pj) - dot(pi_,pi_) - w;
        row   = row + 1;
        A_wbvc(row,:) = a_row;
        b_wbvc(row)   = b_row;
    end
    A_ineq = [A_box; A_wbvc(1:row,:)];
    b_ineq = [b_box; b_wbvc(1:row)];
end

function p_hat = closest_point_in_cell_3d(target, i, pos, radii, gammas, bounds)
    % Project 'target' onto the WBVC polyhedron using QP
    % min  0.5||x - target||²   s.t.  A_ineq*x <= b_ineq
    if point_in_wbvc_3d(target, i, pos, radii, gammas)
        p_hat = target; return;
    end

    [A_ineq, b_ineq] = wbvc_halfplanes(i, pos, radii, gammas, bounds);

    opts = optimoptions('quadprog','Display','off',...
                        'Algorithm','interior-point-convex',...
                        'OptimalityTolerance',1e-16,...
                        'ConstraintTolerance',1e-16,...
                        'StepTolerance',1e-16,...
                        'MaxIterations',20000);
    H    = eye(3);
    f    = -target(:);
    [x, ~, flag] = quadprog(H, f, A_ineq, b_ineq, [], [], [], [], [], opts);
    if flag > 0
        p_hat = x(:)';
    else
        % Fallback: stay at current position
        p_hat = pos(i,:);
    end
end

function new_pt = deadlock_resolve_3d(boundary_pt, i, pos, radii, gammas, bounds, fraction)
    if nargin < 7; fraction = 0.25; end
    % Perturb along a random direction inside the cell
    rng_state = rng; rng('shuffle');
    delta = randn(1,3); delta = delta / max(norm(delta),1e-9);
    rng(rng_state);

    new_pt = boundary_pt + fraction * delta;
    new_pt = max(new_pt, bounds(:,1)' + 0.01);
    new_pt = min(new_pt, bounds(:,2)' - 0.01);
end


%% ══════════════════════════════════════════════════════════════════════════
%%  NEAREST-NEIGHBOUR DISTANCE HISTORY
%% ══════════════════════════════════════════════════════════════════════════
function dist_nn = compute_nn_dist_hist(history)
    T = size(history,1);
    n = size(history,2);
    dist_nn = zeros(T, n);
    for t = 1:T
        pos_t = squeeze(history(t,:,:));
        for i = 1:n
            dists    = vecnorm(pos_t - pos_t(i,:), 2, 2);
            dists(i) = inf;
            dist_nn(t,i) = min(dists);
        end
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%%  TIMING REPORT
%% ══════════════════════════════════════════════════════════════════════════
function print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
                      goal_reach_iter, times, total_sim_time, label)
    n           = numel(gammas);
    avg_ctrl_ms = mean(ctrl_time_accum) * 1e3;

    fprintf('\n══════════════════════════════════════════════════════════════\n');
    fprintf('  TIMING REPORT  –  %s\n', label);
    fprintf('══════════════════════════════════════════════════════════════\n');
    fprintf('  Total wall-clock time            : %.4f s\n', total_sim_time);
    fprintf('  Avg ctrl time / step (all agents): %.4f ms\n', avg_ctrl_ms);
    fprintf('──────────────────────────────────────────────────────────────\n');
    fprintf('  %-6s  %-8s  %-18s\n','Agent','SVO(γ)','Time to goal (s)');
    fprintf('  %-6s  %-8s  %-18s\n','-----','-------','----------------');
    for i = 1:n
        if goal_reach_iter(i) > 0
            g_str = sprintf('%.3f', times(goal_reach_iter(i)));
        else
            g_str = 'Not reached';
        end
        fprintf('  %-6d  %-8.2f  %-18s\n', i, gammas(i), g_str);
    end
    fprintf('══════════════════════════════════════════════════════════════\n');
end


%% ══════════════════════════════════════════════════════════════════════════
%%  VISUALISATION HELPERS – 3-D
%% ══════════════════════════════════════════════════════════════════════════

function c = svo_color(gamma)
    c = [0.20 + 0.65*gamma, 0.40 - 0.25*gamma, 0.80 - 0.65*gamma];
    c = max(0, min(1, c));
end

function setup_axes_3d(ax, bounds, ttl)
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlim(ax, bounds(1,:)); ylim(ax, bounds(2,:)); zlim(ax, bounds(3,:));
    ax.Color    = [0.97 0.97 0.97];
    ax.FontSize = 10;
    title(ax, ttl, 'FontSize', 10);
    xlabel(ax,'$x$ (m)','Interpreter','latex');
    ylabel(ax,'$y$ (m)','Interpreter','latex');
    zlabel(ax,'$z$ (m)','Interpreter','latex');
    view(ax, 35, 25);
end

function draw_wbvc_cells_3d(ax, pos, radii, gammas, bounds)
    % Visualise each WBVC cell as a semi-transparent convex hull.
    % Strategy: solve a set of QPs to find the extreme point of the cell
    % in each of N_DIR uniformly sampled directions on the unit sphere.
    % The resulting point cloud is then wrapped by convhull and rendered.
    % No external toolbox needed – only MATLAB's quadprog.

    n       = size(pos,1);
    N_DIR   = 60;      % directions to probe (higher → smoother cell shape)
    opts    = optimoptions('quadprog','Display','off', ...
                           'Algorithm','interior-point-convex');

    % Uniform directions on S² via Fibonacci lattice
    golden  = (1 + sqrt(5)) / 2;
    idx     = (0:N_DIR-1)';
    theta   = acos(1 - 2*(idx+0.5)/N_DIR);          % polar
    phi     = 2*pi * idx / golden;                   % azimuth
    dirs    = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];  % N_DIR×3

    for i = 1:n
        [A_ineq, b_ineq] = wbvc_halfplanes(i, pos, radii, gammas, bounds);
        gc  = svo_color(gammas(i));

        % Probe extreme points: max d'x  s.t.  Ax<=b  →  min -d'x  s.t.  Ax<=b
        V = zeros(N_DIR, 3);
        ok = false(N_DIR,1);
        H  = zeros(3);          % zero Hessian → LP via quadprog
        lb = bounds(:,1)';
        ub = bounds(:,2)';
        for k = 1:N_DIR
            f = -dirs(k,:)';    % minimise negative projection = maximise projection
            [x, ~, flag] = quadprog(H, f, A_ineq, b_ineq, ...
                                    [], [], lb, ub, pos(i,:)', opts);
            if flag > 0
                V(k,:) = x';
                ok(k)  = true;
            end
        end
        V = V(ok,:);

        % Also include the agent centre so degenerate cells still render
        V = [V; pos(i,:)];     %#ok<AGROW>

        % Deduplicate (within 1e-4 m)
        V = uniquetol(V, 1e-4, 'ByRows', true);

        if size(V,1) < 4; continue; end   % need ≥4 non-coplanar pts

        try
            K = convhull(V(:,1), V(:,2), V(:,3), 'Simplify', true);
            trisurf(K, V(:,1), V(:,2), V(:,3), ...
                    'FaceColor', gc, 'FaceAlpha', 0.08, ...
                    'EdgeColor', gc*0.6, 'EdgeAlpha', 0.20, ...
                    'LineWidth', 0.4, 'HandleVisibility','off', ...
                    'Parent', ax);
        catch
            % convhull fails for near-planar point sets – skip silently
        end
    end
end

function draw_agents_3d(ax, pos, radii, gammas, goals)
    n = size(pos,1);
    [sx,sy,sz] = sphere(20);
    for i = 1:n
        gc = svo_color(gammas(i));
        % Goal marker
        if nargin >= 5 && ~isempty(goals)
            plot3(ax, goals(i,1), goals(i,2), goals(i,3), ...
                  '*', 'Color', gc, 'MarkerSize', 10, 'LineWidth', 1.5, ...
                  'HandleVisibility','off');
        end
        % Agent sphere
        surf(ax, pos(i,1)+radii(i)*sx, pos(i,2)+radii(i)*sy, pos(i,3)+radii(i)*sz, ...
             'FaceColor', gc, 'FaceAlpha', 0.75, 'EdgeColor','none', ...
             'HandleVisibility','off');
    end
end

function draw_trajectories_3d(ax, history, gammas, k_end)
    if nargin < 4; k_end = size(history,1); end
    if k_end < 2; return; end
    n = size(history,2);
    for i = 1:n
        xyz  = squeeze(history(1:k_end, i, :));
        if isvector(xyz); continue; end
        rgba = [svo_color(gammas(i)), 0.75];
        plot3(ax, xyz(:,1), xyz(:,2), xyz(:,3), '-', ...
              'Color', rgba, 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end

function add_legend(ax)
    patch(ax, nan, nan, nan, svo_color(1.0), 'DisplayName','Egoistic (\gamma=1)');
    patch(ax, nan, nan, nan, svo_color(0.5), 'DisplayName','Prosocial (\gamma=0.5)');
    patch(ax, nan, nan, nan, svo_color(0.0), 'DisplayName','Altruistic (\gamma=0)');
    legend(ax,'Location','northeast','FontSize',8);
end