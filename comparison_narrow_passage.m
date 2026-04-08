%% COMPARISON SCRIPT: PROPOSED METHOD (STT) vs. SAFETY BARRIER CERTIFICATES
clc; clear; close all;

%% ========================================================================
%% PART 1: PROPOSED METHOD (STT)
%% ========================================================================
n_agents_stt = 2;
dT_stt = 30;
cen_initial_stt = [150 60; 150 322.5]; 
eta_stt = [150 322.5; 150 60]; 
rad_max = 0.9 * dT_stt; 
rad_min = 0.7 * dT_stt;
rad = rad_max * ones(n_agents_stt, 1);
rad_initial = rad;

h_1 = 0.09; h_2 = 1; hat_h_2 = 8.0; b_stt = 1;
delt_stt = 0.009; tf_stt = 150; tc_stt = 0.9 * tf_stt;
t_stt = 0:delt_stt:tf_stt;
robot_stt = cen_initial_stt;
s_a = [0.9 0.1];

rad_o = 112; dO = 3;
obstacles_stt = [0 191.25 rad_o 0 0; 272 191.25 rad_o 0 0];
theta_circle = linspace(0, 2*pi, 600);

snap_times_stt = [0.05 * tf_stt, 0.1 * tf_stt];
snap_iters_stt = round(snap_times_stt / delt_stt) + 1;
snap_data_stt = struct('tube',{},'robot',{},'rad',{},'obs',{},'tube_traj',{},'time',{});
tube_traj = cell(n_agents_stt, 1);
for i = 1:n_agents_stt, tube_traj{i} = cen_initial_stt(i,:); end

% Loop for STT
for iter = 1:length(t_stt)
    for k = 1:n_agents_stt
        if t_stt(iter) < tc_stt
            gamma = -h_1 * tc_stt * (cen_initial_stt(k,:) - eta_stt(k,:)) / (tc_stt - t_stt(iter));
        else
            gamma = [0 0];
        end
        d_cen_2 = [0,0]; d_cen_3 = [0,0]; rad(k) = rad_max; sum_dist = 1e4; num_obs = 0;
        
        for o = 1:size(obstacles_stt,1)
            obs_pos = obstacles_stt(o,1:2); obs_radius = obstacles_stt(o,3);
            d = norm(cen_initial_stt(k,:) - obs_pos) - (obs_radius + rad_min);
            d_1 = norm(cen_initial_stt(k,:) - obs_pos) - obs_radius;
            if d_1 <= rad_max
                m_j = (1/d_1 - 1/rad_max) * h_2 * (1/d^2) * ((cen_initial_stt(k,:) - obs_pos)/d);
                v_j = (k==1) * ([0 -1;1 0]*m_j')' + (k~=1) * ([0 1;-1 0]*m_j')';
                d_cen_2 = d_cen_2 + m_j + v_j; num_obs = num_obs + 1;
                sum_dist = smooth_min(sum_dist, d_1);
            end
        end
        if num_obs > 0, rad(k) = smooth_min(rad_max, sum_dist); end

        other_idx = [1:k-1, k+1:n_agents_stt];
        for l = 1:length(other_idx)
            idx_l = other_idx(l);
            agent_pos = cen_initial_stt(idx_l,:);
            d_agent = norm(cen_initial_stt(k,:) - agent_pos) - 2*rad_min;
            phi = phi_kl(s_a(k), s_a(idx_l), t_stt(iter), tc_stt, b_stt);
            if norm(cen_initial_stt(k,:) - agent_pos) <= 2*rad_max
                hat_m = (1/norm(cen_initial_stt(k,:) - agent_pos) - 1/(2*rad_max)) * hat_h_2 * (1/d_agent^2) * ((cen_initial_stt(k,:) - agent_pos)/d_agent);
                hat_v = (k==1) * ([0 -1;1 0]*hat_m')' + (k~=1) * ([0 1;-1 0]*hat_m')';
                d_cen_3 = d_cen_3 + phi * (hat_m + hat_v);
            end
        end
        tube_next(k,:) = cen_initial_stt(k,:) + delt_stt * (gamma + d_cen_2 + d_cen_3);
    end
    % Radius and Robot Update
    for k = 1:n_agents_stt
        sum_dist_r = rad(k);
        other_idx = [1:k-1, k+1:n_agents_stt];
        for l = 1:length(other_idx)
            idx_l = other_idx(l);
            phi = phi_kl(s_a(k), s_a(idx_l), t_stt(iter), tc_stt, b_stt);
            d_2k = norm(tube_next(k,:) - tube_next(idx_l,:));
            if d_2k <= 2*rad_max
                w = (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi) + (d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                sum_dist_r = smooth_min(sum_dist_r, rad_min+(d_2k-2*rad_min)*w);
            end
        end
        rad(k) = smooth_min(rad_max, sum_dist_r);
        u_robot = control_stt(tube_next(k,:), rad(k), robot_stt(k,:));
        robot_stt(k,:) = robot_stt(k,:) + delt_stt * u_robot;
        tube_traj{k} = [tube_traj{k}; tube_next(k,:)];
    end
    % Snapshot Storage
    for s = 1:length(snap_iters_stt)
        if iter == snap_iters_stt(s)
            snap_data_stt(s).tube = tube_next; snap_data_stt(s).robot = robot_stt;
            snap_data_stt(s).rad = rad; snap_data_stt(s).obs = obstacles_stt;
            snap_data_stt(s).tube_traj = tube_traj; snap_data_stt(s).time = t_stt(iter);
        end
    end
    cen_initial_stt = tube_next; rad_initial = rad;
end

%% ========================================================================
%% PART 2: SAFETY BARRIER CERTIFICATES (CBF-QP)
%% ========================================================================
n_agents_cbf = 2; dt_cbf = 0.02; T_cbf = 15; steps_cbf = T_cbf/dt_cbf;
Ds = 2*13; gamma_cbf = 0.2; k1 = 0.5; k2 = 1.8; alpha_max = 4000.0; beta_max = 2000.0;
pos0_cbf = [150 60; 150 322.5]; goals_cbf = [150 322.5; 150 60];
obs_static = [0 191.25 112; 272 191.25 112]; n_obs = size(obs_static,1);

pos = pos0_cbf; vel = zeros(n_agents_cbf, 2);
pos_hist = zeros(steps_cbf+1, n_agents_cbf, 2); pos_hist(1,:,:) = pos;
snap_times_cbf = [10, 15]; snap_iters_cbf = round(snap_times_cbf / dt_cbf);
snap_data_cbf = struct('pos',{},'time',{});

for k = 1:steps_cbf
    u_nom = zeros(n_agents_cbf, 2);
    for i = 1:n_agents_cbf
        u_nom(i,:) = -k1*(pos(i,:) - goals_cbf(i,:)) - k2*vel(i,:);
    end
    
    Ndv = 2*n_agents_cbf; H_qp = 2*eye(Ndv); A_ineq = []; b_ineq = [];
    % Constraints Agent-Agent and Agent-Obstacle (logic simplified for brevity)
    for i = 1:n_agents_cbf
        % Obstacles
        for o = 1:n_obs
            Dp = pos(i,:) - obs_static(o,1:2); Dv = vel(i,:); d = norm(Dp);
            Ds_obs = Ds + obs_static(o,3); radicand = 2*alpha_max*(d - Ds_obs);
            if radicand <= 0, radicand = 1e-4; end
            h_io = sqrt(radicand) + dot(Dp, Dv)/d;
            rhs = gamma_cbf * h_io * d - dot(Dp, Dv)^2/d^2 + norm(Dv)^2 + alpha_max * dot(Dp, Dv) / sqrt(radicand + 1e-9);
            row = zeros(1, Ndv); row((i-1)*2 + (1:2)) = -Dp;
            A_ineq = [A_ineq; row]; b_ineq = [b_ineq; rhs];
        end
    end
    
    opts = optimoptions('quadprog','Display','off');
    [u_sol, ~, exitflag] = quadprog(H_qp, -2*reshape(u_nom.',[],1), A_ineq, b_ineq, [], [], [], [], [], opts);
    if exitflag <= 0, u_sol = zeros(Ndv,1); end
    u_mat = reshape(u_sol, 2, n_agents_cbf)';
    vel = vel + dt_cbf * u_mat; pos = pos + dt_cbf * vel;
    pos_hist(k+1,:,:) = pos;
    
    for s = 1:length(snap_iters_cbf)
        if k == snap_iters_cbf(s)
            snap_data_cbf(s).pos = pos; snap_data_cbf(s).time = k*dt_cbf;
        end
    end
    if all(vecnorm(pos - goals_cbf, 2, 2) < 0.1), break; end
end

%% ========================================================================
%% VISUALIZATION
%% ========================================================================
stt_colors = [0.3010 0.7450 0.9330; 0.9290 0.6940 0.1250];

% FIGURE 1: PROPOSED METHOD (STT)
figure(1); set(gcf,'Position',[100 100 1200 500], 'Name', 'Proposed Method (STT)');
for s = 1:2
    subplot(1,2,s); hold on; grid on; axis([0 272 0 382.5]);
    for i=1:2, fill(eta_stt(i,1)+dT_stt*cos(theta_circle), eta_stt(i,2)+dT_stt*sin(theta_circle),'g','FaceAlpha',0.15); end
    for o=1:2, fill(obstacles_stt(o,1)+(obstacles_stt(o,3)-dO)*cos(theta_circle), obstacles_stt(o,2)+(obstacles_stt(o,3)-dO)*sin(theta_circle),'r','FaceAlpha',0.4); end
    if s <= length(snap_data_stt)
        sd = snap_data_stt(s);
        for i=1:2
            plot(sd.tube_traj{i}(:,1), sd.tube_traj{i}(:,2), 'Color', stt_colors(i,:), 'LineWidth', 2);
            plot(sd.tube(i,1)+sd.rad(i)*cos(theta_circle), sd.tube(i,2)+sd.rad(i)*sin(theta_circle), '--', 'Color', stt_colors(i,:));
            plot(sd.robot(i,1), sd.robot(i,2), 'ko', 'MarkerFaceColor','k');
        end
        title(sprintf('Proposed Method: t = %.2f s', sd.time));
    end
end

% FIGURE 2: SAFETY BARRIER CERTIFICATES
figure(2); set(gcf,'Position',[150 150 1200 500], 'Name', 'Safety Barrier Certificates (CBF-QP)');
for s = 1:2
    subplot(1,2,s); hold on; grid on; axis([0 272 0 382.5]);
    for o=1:2, fill(obs_static(o,1)+obs_static(o,3)*cos(theta_circle), obs_static(o,2)+obs_static(o,3)*sin(theta_circle),'r','FaceAlpha',0.4); end
    for i=1:2, fill(goals_cbf(i,1)+30*cos(theta_circle), goals_cbf(i,2)+30*sin(theta_circle),'g','FaceAlpha',0.15); end
    if s <= length(snap_data_cbf)
        n_pts = min(snap_iters_cbf(s)+1, size(pos_hist,1));
        for i=1:2
            xy = squeeze(pos_hist(1:n_pts, i, :));
            plot(xy(:,1), xy(:,2), 'Color', stt_colors(i,:), 'LineWidth', 2);
            plot(snap_data_cbf(s).pos(i,1)+21*cos(theta_circle), snap_data_cbf(s).pos(i,2)+21*sin(theta_circle), '--', 'Color', stt_colors(i,:));
            plot(snap_data_cbf(s).pos(i,1), snap_data_cbf(s).pos(i,2), 'ko', 'MarkerFaceColor','k');
        end
        title(sprintf('CBF-QP Method: t = %.1f s', snap_data_cbf(s).time));
    end
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================
function phi = phi_kl(s_a_k, s_a_l, t, t_c_k, b)
    denom = s_a_l + s_a_k;
    phi = (t < t_c_k) * (s_a_k / denom) + (t >= t_c_k) * (s_a_k / denom * exp(-((t-t_c_k)^2)/(b^2)));
end

function y = smooth_min(x1, x2)
    mu = 1; y = -(1/mu) * log(exp(-mu*x1) + exp(-mu*x2));
end

function u = control_stt(tube_pos, rad_val, robot_pos)
    k_gain = 5; u = k_gain * (tube_pos - robot_pos);
end
