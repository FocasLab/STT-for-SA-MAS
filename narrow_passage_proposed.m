%% %% STT for SA-MAS in Narrow Passage with Snapshots

clc; 
clear; 
clf;

n_agents=2;

%% --- STT Parameters ---
dT=30;
cen_initial = [150 60; 150 322.5]; % starting point of Tube
eta = [150 322.5; 150 60]; % goal point in Target set

rad_max =0.9*dT; % maximum allowable tube radius
rad_min=0.7*dT; % minimum allowable tube radius

rad=rad_max*ones(n_agents,1); % initializing radius of tube as rad_max
rad_initial=rad;

%gain used in (12) and (17)
h_1 = 0.09;
h_2 = 1;
hat_h_2=8.0;
b=1;

%% Simulation Parameters
delt = 0.009;
tf = 150;
tc =0.9*tf;
t = 0:delt:tf;

robot = cen_initial; % robots initiali psoition
s_a=[0.9 0.1]; % social awareness value of eac agent

control_agents=zeros(length(t),2,n_agents);

%% Obstacles
rad_o=112;
dO=3;
obstacles = [0   191.25 rad_o 0 0;
             272 191.25 rad_o 0 0];

theta = linspace(0,2*pi,600);

%%  SNAPSHOT SETTINGS
snap_times = [0.05*tf, 0.1*tf];
snap_iters = round(snap_times/delt)+1;

snap_data = struct('tube',{},'robot',{},'rad',{},...
                   'obs',{},'tube_traj',{},'time',{});

% trajectory storage
tube_traj = cell(n_agents,1);
for i=1:n_agents
    tube_traj{i} = cen_initial(i,:);
end

%% --- Simulation Loop ---
for iter = 1:length(t)

    % Update obstacles
    for o = 1:size(obstacles,1)
        obstacles(o,1) = obstacles(o,1) + delt*obstacles(o,4);
        obstacles(o,2) = obstacles(o,2) + delt*obstacles(o,5);
    end

    %% Tube dynamics
    for k=1:n_agents
        % calculation of first term in (12)
        if t(iter)<tc
            gamma = -h_1 * tc * (cen_initial(k,:) - eta(k,:)) / (tc - t(iter));
        else
            gamma = [0 0];
        end
        d_cen_2 = [0,0];
        d_cen_3 = [0,0];

        rad(k)=rad_max;
        sum_dist = 1e4;
        num_obs=0;

        %  % calculation of second term in (12)
        for o=1:size(obstacles,1)
            obs_pos = obstacles(o,1:2);
            obs_radius = obstacles(o,3);
            d = norm(cen_initial(k,:) - obs_pos) - (obs_radius+rad_min);
            d_1= norm(cen_initial(k,:) - obs_pos)-obs_radius;
            if d_1 <= rad_max % consider only those obstacles inside the sensing radius
                m_j= (1/d_1-1/rad_max)*h_2*(1/d^2)*((cen_initial(k,:) - obs_pos)/d);
                if k==1
                    v_j=([0 -1;1 0]*m_j')';
                else
                    v_j=([0 1;-1 0]*m_j')';
                end
                d_cen_2 = d_cen_2 + m_j + v_j;
                num_obs = num_obs + 1;
                sum_dist = smooth_min(sum_dist,d_1);
            end
        end

        if num_obs>0
            rad(k)=smooth_min(rad_max,sum_dist);
        end
         % calculation of secodnd term in (12)
        Ag_a_k=[cen_initial([1:k-1,k+1:end],:) rad_initial([1:k-1,k+1:end],:)];% set of all agents except the agent k
        s_a_values_k=[s_a([1:k-1,k+1:end])];

        for l=1:size(Ag_a_k,1)

            agent_pos = Ag_a_k(l,1:2);
            d_agent = norm(cen_initial(k,:) - agent_pos) - 2*rad_min;
            phi = phi_kl(s_a(k), s_a_values_k(l), t(iter), tc, b);
            if norm(cen_initial(k,:) - agent_pos) <=2*rad_max % consider only those agnet which are inside sensing radius i.e. l \in N_a_k
                hat_m = (1/norm(cen_initial(k,:) - agent_pos)-1/(2*rad_max)) ...
                        *hat_h_2*(1/d_agent^2)*((cen_initial(k,:) - agent_pos)/d_agent);
                if k==1
                    hat_v=([0 -1;1 0]*hat_m')';
                else
                    hat_v=([0 1;-1 0]*hat_m')';
                end

                d_cen_3 = d_cen_3 + phi*(hat_m+hat_v);
            end
        end

        d_cen = gamma + d_cen_2 + d_cen_3;
        % euler update for center of Tube (12)
        tube(k,:) = cen_initial(k,:) + delt*d_cen;

    end

    %% Radius update and control implementation
    for k=1:n_agents

        sum_dist=rad(k);
        Ag_k_l = tube([1:k-1,k+1:end],:);
        s_vals = s_a([1:k-1,k+1:end]);

        for l=1:size(Ag_k_l,1)

            phi = phi_kl(s_a(k), s_vals(l), t(iter), tc, b);
            d_2k = norm(tube(k,:) - Ag_k_l(l,:));

            if d_2k <= 2*rad_max % consider only those agnet which are inside sensing radius i.e. l \in N_a_k
                w= (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi)+(d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                sum_dist=smooth_min(sum_dist,rad_min+(d_2k-2*rad_min)*w);
            end
        end

        rad(k)=smooth_min(rad_max,sum_dist);
        control_agents(iter,:,k)=control(tube(k,:),rad(k),robot(k,:)); % control law in Theorem 11
        robot(k,:) = robot(k,:) + delt*control_agents(iter,:,k);
        % store trajectory
        tube_traj{k} = [tube_traj{k}; tube(k,:)];
    end

    %% STORE SNAPSHOTS
    for s=1:length(snap_iters)
        if iter == snap_iters(s)
            snap_data(s).tube = tube;
            snap_data(s).robot = robot;
            snap_data(s).rad = rad;
            snap_data(s).obs = obstacles;
            snap_data(s).tube_traj = tube_traj;
            snap_data(s).time = t(iter);
        end
    end

    cen_initial = tube;
    rad_initial = rad;

end

%% SUBPLOT VISUALIZATION

figure;
set(gcf,'Position',[100 0 1600 800]);

colors = [0.3010 0.7450 0.9330;
          0.9290 0.6940 0.1250];

for s=1:length(snap_data)

    subplot(1,2,s); hold on; grid on; box on;
    axis([0 272 0 382.5]);

    sd = snap_data(s);

    % goals
    for i=1:n_agents
        xg = eta(i,1) + dT*cos(theta);
        yg = eta(i,2) + dT*sin(theta);
        fill(xg,yg,'g','FaceAlpha',0.2);
    end

    % obstacles
    for o=1:size(sd.obs,1)
        xo = sd.obs(o,1)+(sd.obs(o,3)-dO)*cos(theta);
        yo = sd.obs(o,2)+(sd.obs(o,3)-dO)*sin(theta);
        fill(xo,yo,'r','FaceAlpha',0.5);
    end

    % trajectories
    for i=1:n_agents
        traj = sd.tube_traj{i};
        plot(traj(:,1),traj(:,2),'Color',colors(i,:),'LineWidth',2);
    end

    % tubes
    for i=1:n_agents
        sx = sd.tube(i,1)+sd.rad(i)*cos(theta);
        sy = sd.tube(i,2)+sd.rad(i)*sin(theta);
        plot(sx,sy,'--','Color',colors(i,:),'LineWidth',2);
    end

    % robots
    for i=1:n_agents
        plot(sd.robot(i,1),sd.robot(i,2),'ko','MarkerFaceColor','k');
    end

    title(sprintf('t = %.2f s',sd.time));
end

%% --- Function ---
function phi = phi_kl(s_a_k, s_a_l, t, t_c_k, b)
    denom = s_a_l + s_a_k;
    if t < t_c_k
        phi = s_a_k / denom;
    else
        phi = (s_a_k / denom) * exp(-((t - t_c_k)^2) / (b^2));
    end
end

%% -- smooth minimum ---
function y = smooth_min(x1, x2)
    mu = 1;
    y = -(1/mu) * log( exp(-mu*x1) + exp(-mu*x2) );

end