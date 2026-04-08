%% Weighted Buffered Voronoi Cells (WBVC)
%  Pierson, Schwarting, Karaman, Rus – ICRA 2020
%
%  Three scenarios:
%    1. Four-agent diagonal swap (equal SVO & mixed SVO)
%    2. Circle swap  (n agents, mixed SVO)
%    3. Random navigation (n agents, random SVO)
%
%  Timing metrics printed at end of each scenario:
%    - Average control computation time per step per agent (ms)
%    - Time to reach goal per agent (s)
%    - Total simulation wall-clock time (s)

clc; clear; close all;

%% ── Run scenarios ─────────────────────────────────────────────────────────
% scenario_four_agent();
scenario_circle_swap(5);
% scenario_random(16);


%% ══════════════════════════════════════════════════════════════════════════
%%  SCENARIO 1 – Four-agent diagonal swap
%% ══════════════════════════════════════════════════════════════════════════
function scenario_four_agent()
    fprintf('\n======  Scenario 1: Four-Agent Swap  ======\n');
    Q  = 20;
    r  = 0.2;
    [positions,goals] = generateAgentsPolygon(4,8,[10 10]); 
    % positions = [0.3 0.3; 1.3 0.3; 1.3 1.3; 0.3 1.3];
    % goals     = [1.3 1.3; 0.3 1.3; 0.3 0.3; 1.3 0.3];
    radii     = r * ones(4,1);

    cases = { [0.5 0.5 0.5 0.5]', 'Equal SVO (BVC)'; ...
              [1.0 0.5 0.0 0.0]', 'Mixed SVO' };

    fig = figure('Name','Scenario 1','NumberTitle','off');
    fig.Position = [50 50 1300 800];

    for row = 1:2
        gammas = cases{row,1};
        label  = cases{row,2};

        [history, times, dist_goal, ctrl_time_accum, ctrl_calls, ...
         goal_reach_iter, total_sim_time] = ...
            run_wbvc(positions, goals, gammas, radii, Q, ...
                     2.0, 0.04, 0.06, 1500, 0.05);

        n_steps = size(history,1) - 1;

        %% ── Print timing ───────────────────────────────────────────────
        print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
                     goal_reach_iter, times, total_sim_time, label);

        %% ── Plot initial ───────────────────────────────────────────────
        ax1 = subplot(2,3, (row-1)*3 + 1);
        setup_axes(ax1, Q, [label ' – Initial']);
        draw_wbvc_cells(ax1, positions, radii, gammas, Q);
        draw_agents(ax1, positions, radii, gammas, goals);
        add_legend(ax1);

        %% ── Plot final + trajectories ──────────────────────────────────
        pos_final = squeeze(history(end,:,:));
        ax2 = subplot(2,3, (row-1)*3 + 2);
        setup_axes(ax2, Q, sprintf('%s – Final (t=%.2fs)', label, times(end)));
        draw_wbvc_cells(ax2, pos_final, radii, gammas, Q);
        draw_agents(ax2, pos_final, radii, gammas, goals);
        draw_trajectories(ax2, history, gammas);

        %% ── Distance to goal plot ──────────────────────────────────────
        ax3 = subplot(2,3, (row-1)*3 + 3);
        hold(ax3,'on'); grid(ax3,'on');
        for i = 1:4
            plot(ax3, times, dist_goal(:,i), 'Color', svo_color(gammas(i)), ...
                 'LineWidth', 1.6, 'DisplayName', sprintf('Agent %d',i));
        end
        yline(ax3, 0.05,'--','Color',[0.5 0.5 0.5],'LineWidth',0.9, ...
              'DisplayName','Goal tol.');
        xlabel(ax3,'Time (s)'); ylabel(ax3,'Distance to Goal (m)');
        title(ax3,'Distance to Goal','FontSize',10);
        legend(ax3,'FontSize',8,'Location','northeast');
        ax3.FontSize = 10;
    end
    sgtitle(fig, 'Scenario 1: Four-Agent Position Swap', ...
            'FontSize',14,'FontWeight','bold');
end


%%%Start Position and Goal Positions%%%
function [startPos, goalPos] = generateAgentsPolygon(n, R, center)

    if nargin < 2
        R = 12;
    end

    if nargin < 3
        center = [0 0];   % default center
    end

    if n > 100 || n < 3
        error('n must be between 3 and 10');
    end

    theta = linspace(0,2*pi,n+1);
    theta(end) = [];

    % Start positions on circle centered at "center"
    startPos = [center(1) + R*cos(theta(:)), ...
                center(2) + R*sin(theta(:))];

    % Goals opposite on the same circle
    goalPos = [center(1) + R*cos(theta(:)+pi), ...
               center(2) + R*sin(theta(:)+pi)];

end
%% ══════════════════════════════════════════════════════════════════════════
%%  SCENARIO 2 – Circle swap
%% ══════════════════════════════════════════════════════════════════════════
function scenario_circle_swap(n_agents)
    fprintf('\n======  Scenario 2: Circle Swap (%d agents)  ======\n', n_agents);
    Q    = 70.0;
    R_c  = 3.8;
    cx   = Q/2;  cy = Q/2;
    r    = 0.20;

    % SVO: 1/3 each
    third  = floor(n_agents/3);
    gammas = [0.9 0.1 0.1 0.1 0.9 0.1 0.1 0.1 0.1*ones(1,n_agents-8)];

    % angles    = linspace(0, 2*pi, n_agents+1)'; angles(end) = [];
    % positions = [cx + R_c*cos(angles), cy + R_c*sin(angles)];
    % goals     = [cx + R_c*cos(angles+pi), cy + R_c*sin(angles+pi)];
    [positions,goals] = generateAgentsPolygon(n_agents,28.69,[31 31]); 
    radii     = r * ones(n_agents,1);

    [history, times, dist_goal, ctrl_time_accum, ctrl_calls, ...
     goal_reach_iter, total_sim_time] = ...
        run_wbvc(positions, goals, gammas, radii, Q, ...
                 1.5, 0.05, 0.10, 3000, 0.08);

    print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
                 goal_reach_iter, times, total_sim_time, 'Circle Swap');

    %%--------------------------Four snapshots --------------------------
    T      = size(history,1);
    snaps  = unique([1, max(1,round(T/10)), max(1,round(T/2)), T]);
    labels = arrayfun(@(s) sprintf('t=%.1fs', times(s)), snaps, 'UniformOutput',false);

    fig1 = figure('Name','Scenario 2 – Snapshots','NumberTitle','off');
    fig1.Position = [80 80 1200 1100];
    for k = 1:min(4,numel(snaps))
        ax = subplot(2,2,k);
        pos_snap = squeeze(history(snaps(k),:,:));
        setup_axes(ax, Q, labels{k});
        draw_wbvc_cells(ax, pos_snap, radii, gammas, Q);
        draw_agents(ax, pos_snap, radii, gammas, goals);
%%%-------------- Plot final + trajectories--------------------------
        draw_trajectories(ax, history, gammas,snaps(k));

        if k == 4; add_legend(ax); end
    end
    sgtitle(fig1,'Scenario 2: Circle-Swap Snapshots','FontSize',13,'FontWeight','bold');

    %% --------------------------Distance metrics--------------------------
    dist_nn = compute_nn_dist_hist(history);

    fig2 = figure('Name','Scenario 2 – Distances','NumberTitle','off');
    fig2.Position = [100 100 1100 420];
    ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
    for i = 1:n_agents
        c = [svo_color(gammas(i)), 0.8];
        plot(ax1, times, dist_goal(:,i), 'Color',c,'LineWidth',0.9);
        plot(ax2, times, dist_nn(:,i),   'Color',c,'LineWidth',0.9);
    end
    yline(ax2, 2*r,'k--','LineWidth',1.2,'DisplayName','2r');
    xlabel(ax1,'Time (s)'); ylabel(ax1,'Distance (m)');
    title(ax1,'Distance to Goal'); ax1.FontSize=11;
    xlabel(ax2,'Time (s)'); ylabel(ax2,'Distance (m)');
    title(ax2,'Nearest Neighbour Distance'); legend(ax2,'FontSize',9); ax2.FontSize=11;
    add_legend(ax1);
    sgtitle(fig2,'Scenario 2: Distance Metrics','FontSize',13,'FontWeight','bold');
end


%% ══════════════════════════════════════════════════════════════════════════
%%  SCENARIO 3 – Random navigation
%% ══════════════════════════════════════════════════════════════════════════
function scenario_random(n_agents)
    fprintf('\n======  Scenario 3: Random Navigation (%d agents)  ======\n', n_agents);
    rng(42);
    Q   = 9.0;
    r   = 0.10;
    svo_levels = [0.0 0.2 0.4 0.6 0.8 1.0];

    % Random SVO per agent
    gammas    = svo_levels(randi(numel(svo_levels), n_agents, 1))';
    positions = random_collision_free(n_agents, r, Q);
    perm      = randperm(n_agents);
    goals     = positions(perm,:);
    radii     = r * ones(n_agents,1);
    dt =0.002;
    k_gain = 0.5;
    [history, times, dist_goal, ctrl_time_accum, ctrl_calls, ...
     goal_reach_iter, total_sim_time] = ...
        run_wbvc(positions, goals, gammas, radii, Q, ...
                 k_gain, dt, 0.08, 100/dt, 0.06);

    print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
                 goal_reach_iter, times, total_sim_time, 'Random Navigation');

    %% ── Initial & final snapshot ───────────────────────────────────────
    fig1 = figure('Name','Scenario 3 – Navigation','NumberTitle','off');
    fig1.Position = [120 120 1200 560];
    for col = 1:2
        ax = subplot(1,2,col);
        if col == 1
            pos_show = positions;
            ttl = 'Initial (t=0s)';
        else
            pos_show = squeeze(history(end,:,:));
            ttl = sprintf('Final (t=%.1fs)', times(end));
        end
        setup_axes(ax, Q, ttl);
        draw_wbvc_cells(ax, pos_show, radii, gammas, Q);
        draw_agents(ax, pos_show, radii, gammas, goals);
        if col == 2
            draw_trajectories(ax, history, gammas);
            add_legend(ax);
        end
    end
    sgtitle(fig1,'Scenario 3: Random Navigation','FontSize',13,'FontWeight','bold');

    %% ── Distance metrics ───────────────────────────────────────────────
    dist_nn = compute_nn_dist_hist(history);

    fig2 = figure('Name','Scenario 3 – Distances','NumberTitle','off');
    fig2.Position = [140 140 1100 420];
    ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
    for i = 1:n_agents
        c = [svo_color(gammas(i)), 0.8];
        plot(ax1, times, dist_goal(:,i), 'Color',c,'LineWidth',0.7);
        plot(ax2, times, dist_nn(:,i),   'Color',c,'LineWidth',0.7);
    end
    yline(ax2, 2*r,'k--','LineWidth',1.2,'DisplayName','2r');
    xlabel(ax1,'Time (s)'); ylabel(ax1,'Distance (m)');
    title(ax1,'Distance to Goal'); ax1.FontSize=11;
    xlabel(ax2,'Time (s)'); ylabel(ax2,'Distance (m)');
    title(ax2,'Nearest Neighbour Distance'); legend(ax2,'FontSize',9); ax2.FontSize=11;
    add_legend(ax1);
    sgtitle(fig2,'Scenario 3: Distance Metrics','FontSize',13,'FontWeight','bold');
end


%% ══════════════════════════════════════════════════════════════════════════
%%  CORE WBVC RUNNER  –  returns history + timing data
%% ══════════════════════════════════════════════════════════════════════════
function [history, times, dist_goal_hist, ctrl_time_accum, ctrl_calls, ...
          goal_reach_iter, total_sim_time] = ...
    run_wbvc(pos0, goals, gammas, radii, Q, k_gain, dt, eps_dead, max_steps, goal_tol)

    n = size(pos0,1);
    pos = pos0;

    % Pre-allocate history
    history       = zeros(max_steps+1, n, 2);
    history(1,:,:) = pos;
    times          = zeros(max_steps+1, 1);
    dist_goal_hist = zeros(max_steps+1, n);
    dist_goal_hist(1,:) = vecnorm(pos - goals, 2, 2)';

    % Timing
    ctrl_time_accum = zeros(max_steps,1);
    ctrl_calls      = zeros(max_steps,1);
    goal_reach_iter = -ones(n,1);

    sim_start = tic;   % ← wall-clock timer

    step = 0;
    for iter = 1:max_steps
        new_pos = pos;
        t_ctrl = tic;
        for i = 1:n
            if norm(pos(i,:) - goals(i,:)) < goal_tol
                continue
            end

            %% ── Time control computation ──────────────────────────────
            

            ui    = k_gain * (goals(i,:) - pos(i,:));
            p_bar = pos(i,:) + ui * dt;

            if point_in_wbvc(p_bar, i, pos, radii, gammas)
                p_target = p_bar;
            else
                p_hat = closest_point_in_cell(goals(i,:), i, pos, radii, gammas, Q);
                if norm(pos(i,:) - p_hat) < eps_dead
                    p_hat = deadlock_resolve(p_hat, i, pos, radii, gammas, Q);
                end
                ui2      = k_gain * (p_hat - pos(i,:));
                p_target = pos(i,:) + ui2 * dt;
            end
            p_target = max(p_target, radii(i));
            p_target = min(p_target, Q - radii(i));

            
            %% ─────────────────────────────────────────────────────────

            new_pos(i,:) = p_target;

            if goal_reach_iter(i) == -1 && ...
                    norm(new_pos(i,:) - goals(i,:)) < goal_tol
                goal_reach_iter(i) = iter;
            end
        end

        ctrl_time_accum(iter) = ctrl_time_accum(iter) + toc(t_ctrl);
        ctrl_calls(iter)      = ctrl_calls(iter) + 1;
        pos   = new_pos;
        step  = step + 1;
        history(step+1,:,:)     = pos;
        times(step+1)           = iter * dt;
        dist_goal_hist(step+1,:)= vecnorm(pos - goals, 2, 2)';

        if all(vecnorm(pos - goals, 2, 2) < goal_tol)
            break
        end
    end

    total_sim_time = toc(sim_start);   % ← stop wall-clock

    % Trim to actual length
    history        = history(1:step+1,:,:);
    times          = times(1:step+1);
    dist_goal_hist = dist_goal_hist(1:step+1,:);
end


%% ══════════════════════════════════════════════════════════════════════════
%%  TIMING REPORT
%% ══════════════════════════════════════════════════════════════════════════
function print_timing(gammas, ctrl_time_accum, ctrl_calls, ...
                      goal_reach_iter, times, total_sim_time, label)
    n = numel(gammas);
    avg_ctrl_ms = mean(ctrl_time_accum) * 1e3;

    fprintf('\n══════════════════════════════════════════════════════════════\n');
    fprintf('  TIMING REPORT  –  %s\n', label);
    fprintf('══════════════════════════════════════════════════════════════\n');
    fprintf('  Total simulation wall-clock time : %.4f s\n', total_sim_time);
    fprintf('──────────────────────────────────────────────────────────────\n');
    fprintf('  %-6s  %-8s  %-26s  %-18s\n', ...
            'Agent', 'SVO(γ)', 'Avg ctrl time/step (ms)', 'Time to goal (s)');
    fprintf('  %-6s  %-8s  %-26s  %-18s\n', ...
            '-----', '-------', '--------------------------', '------------------');
    % for i = 1:n
    %     if goal_reach_iter(i) > 0
    %         g_str = sprintf('%.3f', times(goal_reach_iter(i)));
    %     else
    %         g_str = 'Not reached';
    %     end
    %     fprintf('  %-6d  %-8.1f  %-26.4f  %-18s\n', ...
    %             i, gammas(i), avg_ctrl_ms(i), g_str);
    % end
    fprintf('──────────────────────────────────────────────────────────────\n');
    fprintf('  Overall avg ctrl time/step (ms)  : %.4f\n',(avg_ctrl_ms));
    fprintf('══════════════════════════════════════════════════════════════\n');
end


%% ══════════════════════════════════════════════════════════════════════════
%%  WBVC MATHS
%% ══════════════════════════════════════════════════════════════════════════

function w = wij_weight(pi, pj, ri, rj, gi, gj)
    dij = norm(pi - pj);
    if dij < 1e-12; w = 0; return; end
    dg = gi - gj;
    w  = (1 + 0.5*dg)*(ri + rj)*dij - 0.5*dg*dij^2;
end

function inside = point_in_halfplane(q, pi, pj, ri, rj, gi, gj)
    w      = wij_weight(pi, pj, ri, rj, gi, gj);
    inside = (dot(pi-q, pi-q) <= dot(pj-q, pj-q) - w + 1e-9);
end

function inside = point_in_wbvc(q, i, pos, radii, gammas)
    n = size(pos,1);
    for j = 1:n
        if j == i; continue; end
        if ~point_in_halfplane(q, pos(i,:), pos(j,:), ...
                               radii(i), radii(j), gammas(i), gammas(j))
            inside = false; return;
        end
    end
    inside = true;
end

function poly = compute_wbvc_polygon(i, pos, radii, gammas, Q)
    % Sutherland-Hodgman clipping with domain + WBVC half-planes
    poly = [0 0; Q 0; Q Q; 0 Q];
    hp_a = [-1 0; 1 0; 0 -1; 0 1];
    hp_b = [0; Q; 0; Q];
    for k = 1:4
        poly = clip_polygon(poly, hp_a(k,:), hp_b(k));
        if isempty(poly); return; end
    end
    n = size(pos,1);
    for j = 1:n
        if j == i; continue; end
        pj  = pos(j,:);   pi_ = pos(i,:);
        dij = norm(pi_ - pj);
        if dij < 1e-12; continue; end
        w = wij_weight(pi_, pj, radii(i), radii(j), gammas(i), gammas(j));
        a = 2*(pj - pi_);
        b = dot(pj,pj) - dot(pi_,pi_) - w;
        poly = clip_polygon(poly, a, b);
        if isempty(poly); return; end
    end
end

function poly_out = clip_polygon(poly, a, b)
    if isempty(poly); poly_out = []; return; end
    n   = size(poly,1);
    out = zeros(2*n, 2);
    cnt = 0;
    for k = 1:n
        curr = poly(k,:);
        nxt  = poly(mod(k,n)+1,:);
        c_in = (dot(a,curr) <= b + 1e-9);
        n_in = (dot(a,nxt)  <= b + 1e-9);
        if c_in; cnt = cnt+1; out(cnt,:) = curr; end
        if c_in ~= n_in
            denom = dot(a, nxt-curr);
            if abs(denom) > 1e-12
                t = (b - dot(a,curr)) / denom;
                cnt = cnt+1; out(cnt,:) = curr + t*(nxt-curr);
            end
        end
    end
    poly_out = out(1:cnt,:);
end

function p_hat = closest_point_in_cell(target, i, pos, radii, gammas, Q)
    if point_in_wbvc(target, i, pos, radii, gammas)
        p_hat = target; return;
    end
    poly = compute_wbvc_polygon(i, pos, radii, gammas, Q);
    if size(poly,1) < 2; p_hat = pos(i,:); return; end

    pi_  = pos(i,:);
    dir  = target - pi_;
    dn   = norm(dir);
    if dn < 1e-12; p_hat = pi_; return; end
    dir  = dir / dn;

    best_pt = []; best_t = inf;
    np = size(poly,1);
    for k = 1:np
        ea = poly(k,:);
        eb = poly(mod(k,np)+1,:);
        [pt, t] = ray_segment_intersect(pi_, dir, ea, eb);
        if ~isempty(pt) && t > 0 && t < best_t
            best_t = t; best_pt = pt;
        end
    end
    if isempty(best_pt); best_pt = mean(poly,1); end
    p_hat = best_pt;
end
% function p_hat = closest_point_in_cell(target, i, pos, radii, gammas, Q)
%     % True closest point via QP projection — accurate in all cases
%     % min  0.5||x - target||^2   s.t.  A_ineq*x <= b_ineq
%     if point_in_wbvc(target, i, pos, radii, gammas)
%         p_hat = target; return;
%     end
% 
%     n  = size(pos,1);
%     nd = 2;
% 
%     % Box constraints
%     A_box = [eye(nd); -eye(nd)];
%     b_box = [Q; Q; 0; 0];
% 
%     % WBVC halfplane constraints
%     A_wbvc = zeros(n-1, nd);
%     b_wbvc = zeros(n-1, 1);
%     row = 0;
%     for j = 1:n
%         if j == i; continue; end
%         pj  = pos(j,:); pi_ = pos(i,:);
%         dij = norm(pi_ - pj);
%         if dij < 1e-12; continue; end
%         w     = wij_weight(pi_, pj, radii(i), radii(j), gammas(i), gammas(j));
%         a_row = 2*(pj - pi_);
%         b_row = dot(pj,pj) - dot(pi_,pi_) - w;
%         row   = row + 1;
%         A_wbvc(row,:) = a_row;
%         b_wbvc(row)   = b_row;
%     end
% 
%     A_ineq = [A_box; A_wbvc(1:row,:)];
%     b_ineq = [b_box; b_wbvc(1:row)];
% 
%     % opts = optimoptions('quadprog', 'Display', 'off', ...
%     %                     'Algorithm', 'interior-point-convex');
%     opts = optimoptions('quadprog','Display','off',...
%                         'Algorithm','interior-point-convex',...
%                         'OptimalityTolerance',1e-8,...
%                         'ConstraintTolerance',1e-8,...
%                         'StepTolerance',1e-8,...
%                         'MaxIterations',2000);
%     H = eye(2);
%     f = -target(:);
%     [x, ~, flag] = quadprog(H, f, A_ineq, b_ineq, [], [], [], [], [], opts);
%     if flag > 0
%         p_hat = x(:)';
%     else
%         p_hat = pos(i,:);   % fallback: stay put
%     end
% end

function [pt, t] = ray_segment_intersect(origin, dir, a, b)
    v     = b - a;
    denom = dir(1)*v(2) - dir(2)*v(1);
    if abs(denom) < 1e-12; pt = []; t = inf; return; end
    w = origin - a;
    t = (w(2)*v(1)   - w(1)*v(2))   / denom;
    s = (w(2)*dir(1) - w(1)*dir(2)) / denom;
    if t > -1e-9 && s >= -1e-9 && s <= 1+1e-9
        pt = origin + t*dir;
    else
        pt = []; t = inf;
    end
end

function new_pt = deadlock_resolve(boundary_pt, i, pos, radii, gammas, Q, fraction)
    if nargin < 7; fraction = 0.25; end
    poly = compute_wbvc_polygon(i, pos, radii, gammas, Q);
    if size(poly,1) < 2; new_pt = boundary_pt; return; end
    np = size(poly,1);
    best_d = inf; best_k = 1;
    for k = 1:np
        d = point_to_seg_dist(boundary_pt, poly(k,:), poly(mod(k,np)+1,:));
        if d < best_d; best_d = d; best_k = k; end
    end
    ea = poly(best_k,:);
    eb = poly(mod(best_k,np)+1,:);
    new_pt = boundary_pt + fraction*(eb - ea);
    new_pt = max(new_pt, 0.01);
    new_pt = min(new_pt, Q - 0.01);
end

function d = point_to_seg_dist(p, a, b)
    ab = b - a; l2 = dot(ab,ab);
    if l2 < 1e-12; d = norm(p-a); return; end
    t = max(0, min(1, dot(p-a, ab)/l2));
    d = norm(p - (a + t*ab));
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
            dists = vecnorm(pos_t - pos_t(i,:), 2, 2);
            dists(i) = inf;
            dist_nn(t,i) = min(dists);
        end
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%%  VISUALISATION HELPERS
%% ══════════════════════════════════════════════════════════════════════════

function c = svo_color(gamma)
    % blue (altruistic) → purple (prosocial) → red (egoistic)
    c = [0.20 + 0.65*gamma, ...
         0.40 - 0.25*gamma, ...
         0.80 - 0.65*gamma];
    c = max(0, min(1, c));
end

function setup_axes(ax, Q, ttl)
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    axis(ax, [0 Q 0 Q]);
    ax.Color = [0.97 0.97 0.97];
    ax.FontSize = 10;
    title(ax, ttl, 'FontSize', 10);
    xlabel(ax,'$x$ (m)','Interpreter','latex');
    ylabel(ax,'$y$ (m)','Interpreter','latex');
end

function draw_wbvc_cells(ax, pos, radii, gammas, Q)
    n = size(pos,1);
    theta_c = linspace(0, 2*pi, 80);
    for i = 1:n
        poly = compute_wbvc_polygon(i, pos, radii, gammas, Q);
        if size(poly,1) >= 3
            gc = svo_color(gammas(i));
            patch(ax, poly(:,1), poly(:,2), gc, ...
                  'FaceAlpha', 0.22, 'EdgeColor', 'k', ...
                  'LineWidth', 0.7, 'HandleVisibility','off');
        end
    end
end

function draw_agents(ax, pos, radii, gammas, goals)
    n = size(pos,1);
    theta_c = linspace(0, 2*pi, 60);
    for i = 1:n
        gc = svo_color(gammas(i));
        % Goal marker
        if nargin >= 5 && ~isempty(goals)
            plot(ax, goals(i,1), goals(i,2), '*', 'Color', gc, ...
                 'MarkerSize', 9, 'LineWidth', 1.2, ...
                 'HandleVisibility','off');
        end
        % Agent body
        fill(ax, pos(i,1) + radii(i)*cos(theta_c), ...
                 pos(i,2) + radii(i)*sin(theta_c), ...
             gc, 'FaceAlpha', 0.85, 'EdgeColor','k', ...
             'LineWidth', 0.8, 'HandleVisibility','off');
        plot(ax, pos(i,1), pos(i,2), 'k+', 'MarkerSize', 5, ...
             'LineWidth', 1.0, 'HandleVisibility','off');
    end
end

function draw_trajectories(ax, history, gammas, k_end)
    if nargin < 4
        k_end = size(history,1);
    end
    
    if k_end < 2
        return;   % ← nothing to draw yet, exit silently
    end
    
    n = size(history,2);
    for i = 1:n
        xy = squeeze(history(1:k_end, i, :));   % k_end×2 matrix
        
        % Guard: squeeze can collapse to vector if k_end=2 on some MATLAB versions
        if isvector(xy)
            xy = xy(:)';   % force to 1×2, still skip (covered by k_end<2 above)
            continue
        end
        
        rgba = [svo_color(gammas(i)), 0.75];
        plot(ax, xy(:,1), xy(:,2), '-', 'Color', rgba, ...
             'LineWidth', 1.1, 'HandleVisibility','off');
    end
end
function add_legend(ax)
    patch(ax, nan, nan, svo_color(1.0), 'DisplayName','Egoistic (\gamma=1)');
    patch(ax, nan, nan, svo_color(0.5), 'DisplayName','Prosocial (\gamma=0.5)');
    patch(ax, nan, nan, svo_color(0.0), 'DisplayName','Altruistic (\gamma=0)');
    legend(ax, 'Location','northeast','FontSize',8);
end


%% ══════════════════════════════════════════════════════════════════════════
%%  RANDOM COLLISION-FREE PLACEMENT (Scenario 3)
%% ══════════════════════════════════════════════════════════════════════════
function positions = random_collision_free(n, r, Q)
    margin = r * 2.5;
    positions = zeros(n,2);
    placed = 0;
    max_tries = 8000;
    for attempt = 1:max_tries
        if placed == n; break; end
        p = margin + (Q - 2*margin) * rand(1,2);
        ok = true;
        for k = 1:placed
            if norm(p - positions(k,:)) < 2.5*r
                ok = false; break;
            end
        end
        if ok; placed = placed+1; positions(placed,:) = p; end
    end
    if placed < n
        % Grid fallback
        cols = ceil(sqrt(n));
        xs   = linspace(margin, Q-margin, cols);
        for idx = 1:n
            positions(idx,:) = [xs(mod(idx-1,cols)+1), xs(floor((idx-1)/cols)+1)];
        end
    end
end
