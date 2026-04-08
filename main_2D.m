%% STT for SA_MAS 2D Mobile Robot Case
clc; clear; clf; close all

n_a=8; % number of agents in systems

%% STT Initialization
[cen_initial,eta] = generateAgentsPolygon(n_a,8);% start s for tube and eta target point for each agents
dT = 1.2; %Target area
h_1 = [0.22 0.15 0.13 0.13 0.12 0.12 0.12 0.12];  % constant for goal directed term in (12), same for all agents
h_2 = 0.04; %constant for obstacle avoidance and inter agent collison avoidance, same for all agents
rad_max =0.9; % maximum allowable tube radius, same for all agents
rad=rad_max*ones(n_a,1); % radius initialization
rad_p=rad;
rad_min=0.6; %minimum allowable tube radius, same for all agents
b = 1; % rate of decrease in SIF function


%% Agent Initialization
robot = cen_initial; % Start position
delt = 2*1e-2; % Step size for movement
tf = 12;
snapTimes = [0 3 5 10.18];  
snapIndex = 1; 
tc =[5; 6; 7; 8; 9; 10; 11; 12];% Prescribed time for each agent
t = 0:delt:tf;
s_a=[0.1 0.99*ones(1,3) 0.1  0.99*ones(1,3)]'; % social awareness index of each agent

%% Store trajectories
robot_history = cell(n_a,1);
tube_history  = cell(n_a,1);
for i=1:n_a
    robot_history{i} = robot(i,:);
    tube_history{i}  = cen_initial(i,:);
end

%% Plot setup
figure(1)
hold on; 
grid on;
box on;
axis([-15 15 -15 15]); 
axis square
ax = gca;
ax.FontSize = 18;
xlabel('$x $ (m)',Interpreter='latex'); ylabel('$y$ (m)',Interpreter='latex');

for i=1:n_a
    % Draw Target
    theta = linspace(0, 2*pi, 100);
    xc = eta(i,1); yc = eta(i,2);
    x = xc + dT*cos(theta);
    y = yc + dT*sin(theta);
    fill(x, y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'g');
    
    if i==1 || i==5
        tube_plot(i) = plot(cen_initial(i,1), cen_initial(i,2), 'o','Color',[0.9290, 0.6940, 0.1250]);
        trajectory(i) = plot(robot(i,1), robot(i,2), '-','Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
        robot_plot(i) = plot(robot(i,1), robot(i,2), 'o', 'Color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
    else
        tube_plot(i) = plot(cen_initial(i,1), cen_initial(i,2), 'bo');
        trajectory(i) = plot(robot(i,1), robot(i,2), 'b-', 'LineWidth', 1.5);
        robot_plot(i) = plot(robot(i,1), robot(i,2), 'ko', 'MarkerFaceColor', 'k');
    end
    if i==1  || i==5
        safety_circle_plot(i) = plot(nan, nan, '--', 'Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
    else
        safety_circle_plot(i) = plot(nan, nan, '--', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.5);
    end
end

%% STT loop
for iter = 1:length(t)
    current_time = t(iter);
    % STT Update
    for k=1:n_a
        if t(iter)<tc(k)
            gamma = -h_1(k) * tc(k) * (cen_initial(k,:) - eta(k,:)) / (tc(k) - t(iter));
        else
            gamma=0;
        end
        d_cen_agent = [0, 0];
        sum_dist = 10^4;
        Ag_a_k=[cen_initial([1:k-1,k+1:end],:) rad_p([1:k-1,k+1:end],:)];% set of all agents except the agent k
        s_a_values_k=[s_a([1:k-1,k+1:end],:)];
        for l = 1:size(Ag_a_k, 1)
            agent_pos = Ag_a_k(l,1:2);
            d_agent= norm(cen_initial(k,:) - agent_pos) - 2*(rad_min); 
            if norm(cen_initial(k,:) - agent_pos) <= 2*rad_max  % consider only those agnet which are inside sensing radius i.e. l \in N_a_k
                hat_m_kl= (1/(norm(cen_initial(k,:)-agent_pos))-1/(2*rad_max)) *h_2 * (1 / d_agent^2) * ((cen_initial(k,:) - agent_pos) / d_agent);        
                hat_v_kl=([0 1;-1 0]*hat_m_kl')';
                phi=phi_kl(s_a(k), s_a_values_k(l), t(iter), tc(k), b);
                d_cen_agent =d_cen_agent + phi*(hat_m_kl+hat_v_kl);
            end
        end
        % Total d_cen (12) and STT update using euler update
        d_cen_total = gamma +d_cen_agent;
        if norm(d_cen_total) > 0
            tube(k,:) = cen_initial(k,:) + delt * (d_cen_total );
        end
    end

    % Update radius of Tube
    for k=1:n_a
        sum_rad=10^6;
        Ag_a_k=[tube([1:k-1,k+1:end],:)];% set of all agents except the agent k
        s_a_values_rad=[s_a([1:k-1,k+1:end],:)];
        for l = 1:size(Ag_a_k, 1)
            phi=phi_kl(s_a(k), s_a_values_rad(l), t(iter), t(k), b);
            d_2k=norm(tube(k,:)-Ag_a_k(l,:));
            if d_2k<=2*rad_max % consider only those agnet which are inside sensing radius i.e. l \in N_a_k
                w= (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi)+(d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
                sum_dist=smooth_min(sum_dist,rad_min+(d_2k-2*rad_min)*w);
            else
                 sum_rad=smooth_min(sum_rad,rad_max);
            end
        end
        rad(k)=min(rad_max,sum_rad);

        %% Control Implementation
        robot(k,:) = robot(k,:) + delt * control(tube(k,:),rad(k),robot(k,:));

        % store histories
        robot_history{k} = [robot_history{k}; robot(k,:)];
        tube_history{k}  = [tube_history{k};  tube(k,:)];
    end

    % ---- Take snapshot at specified times ----
    if snapIndex <= length(snapTimes) && abs(t(iter) - snapTimes(snapIndex)) < 1e-6
        figure(2);
        subplot(2,2,snapIndex); hold on; grid on; axis equal;
        axis([-10 10 -10 10]);
        title(['t = ' num2str(t(iter),'%.1fs')]);
       xlabel('$x $ (m)',Interpreter='latex'); ylabel('$y$ (m)',Interpreter='latex');

        % --- Draw goals ---
        theta_c = linspace(0,2*pi,100);
        for i = 1:n_a
            r_goal = 1.2;
            xc = eta(i,1); yc = eta(i,2);
            fill(xc + r_goal*cos(theta_c), yc + r_goal*sin(theta_c), ...
                 'g','FaceAlpha',0.2,'EdgeColor','g');
        end

        for k_snap = 1:n_a
            % tube trajectory
            plot(tube_history{k_snap}(:,1), tube_history{k_snap}(:,2), ...
                 '--', 'Color', get(tube_plot(k_snap),'Color'), 'LineWidth', 1.2);
            % robot trajectory
            plot(robot_history{k_snap}(:,1), robot_history{k_snap}(:,2), ...
                 '-', 'Color', get(trajectory(k_snap),'Color'), 'LineWidth', 1.5);
            % current robot position
             if k_snap==1 || k_snap==5
                    plot(robot(k_snap,1), robot(k_snap,2), 'o', 'Color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
             else
                 plot(robot(k_snap,1), robot(k_snap,2), 'ko', 'MarkerFaceColor', 'k');

             end
           
            % current tube position
            plot(tube(k_snap,1), tube(k_snap,2), 'o', 'Color', get(tube_plot(k_snap),'Color'));
            % safety circle
            safety_x = tube(k_snap,1) + (rad(k_snap)-0.1)*cos(theta_c);
            safety_y = tube(k_snap,2) + (rad(k_snap)-0.1)*sin(theta_c);
            plot(safety_x, safety_y, '--', 'Color', get(safety_circle_plot(k_snap),'Color'));
        end
        snapIndex = snapIndex + 1;
    end

    % check wheather all agent reached goal
    if norm(tube - eta) < 0.3
        break
    end
    cen_initial=tube;
    rad_p=rad;
end


% initialization of each agent
function [startPos, goalPos] = generateAgentsPolygon(n, R)
    if nargin < 2, R = 12; end
    if n > 10 || n < 3, error('n must be between 3 and 10'); end
    theta = linspace(0, 2*pi, n+1); theta(end) = [];  
    startPos = R * [cos(theta(:)), sin(theta(:))];
    goalPos  = -startPos;
end

%% SIF Function
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
