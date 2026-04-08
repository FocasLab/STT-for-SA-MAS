%%STT for SA_MAS 3D Quadcoptor Case with acceleration input
clc;clear; clf;

n_a=8; % number of agents
%% STT Initialization
cen = [ 0 0 0; % starting point 's' for the tube
         8 8 8;
         0 8 0;
         8 0  8;
         0,0 8;
         8 8 0;
         8 0 0
         0 8 8
          ];
cen_initial=cen;
eta =[8 8 8; % goal position 'eta' for the tube
    0 0 0;
    8 0 8;
    0 8 0;
    8 8 0;
    0 0 8;
    0 8 8
    8 0 0
    ];    
rad_max = 0.9; % maximum allowable tube radius
rad=0.9*ones(n_a,1);% Initialization of tube radius
rad_past=rad; 
rad_min=0.6; %minimum allowable tube radius
h_1 = [0.07 0.06*ones(n_a-1,1)']; % gain for target driven term
hat_h_2 = 0.002;  % gains for inter agent avoidance term
hat_h_3 = 0.002;  % gains for inter agent avoidance term

b = 1;  %rate at which SIF function decreases

%% Agent Initialization
robot = cen;
velo=zeros(n_a,3);
s_a=[0.1 0.99 0.99 0.99 0.99 0.1 0.99 0.99]'; % social awareness index for each agents
tc =[20 25*ones(n_a-1,1)']; % prescribed time of convergence for each agents
%% Simulation Parmeter
delt = 2e-2;             % Time step
tf = 30;
t = 0:delt:tf;



snapTimes = [0 8.5 14 24.3];  
snapIndex = 1; 
robot_history = cell(n_a,1);
tube_history  = cell(n_a,1);
for i=1:n_a
    robot_history{i} = robot(i,:);
    tube_history{i}  = cen(i,:);
end
for i =1:n_a
    % Initialize robot plot
        if i==1 || i==6
        robot_plot(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'o','Color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
        tube_plot(i) = plot3(cen(i,1), cen(i,2), cen(i,3), 'o','Color',[0.9290, 0.6940, 0.1250]);
        trajectory(i) = plot3(robot(i,1), robot(i,2), robot(i,3), '-','Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
        [sx, sy, sz] = sphere(20);   % finer resolution
        hSafe(i) = surf(nan(size(sx)), nan(size(sy)), nan(size(sz)), ...
           'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
           'FaceColor', [0.9290, 0.6940, 0.1250], ...
           'Tag', 'safe_sphere');
    else
        robot_plot(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'ko', 'MarkerFaceColor', 'k');
        tube_plot(i) = plot3(cen(i,1), cen(i,2), cen(i,3), 'bo');
        trajectory(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'b-', 'LineWidth', 1.5);
        [sx, sy, sz] = sphere(20);   % finer resolution
        hSafe(i) = surf(nan(size(sx)), nan(size(sy)), nan(size(sz)), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'FaceColor', [0.3010 0.7450 0.9330], ...
            'Tag', 'safe_sphere');
    end
end
%% Main Simulation Loop
for iter = 1:length(t)
    for k=1:n_a
        % calculation of first term of (12)
        if t(iter)<tc(k)
            gamma = -h_1(k) * tc(k) * (cen_initial(k,:) - eta(k,:)) / (tc(k) - t(iter));
        else
            gamma=0;
        end
        d_cen_3 = [0, 0, 0]; % stores third term of (12)
     
        Ag_k_l= [cen_initial([1:k-1,k+1:end],:) rad_past([1:k-1,k+1:end],:)]; % all agents except k-th agents
        s_a_values_k=[s_a([1:k-1,k+1:end],:)];
        for l=1:size(Ag_k_l,1)
           agent_pos=Ag_k_l(l,1:3);
           agents_radius=Ag_k_l(l,4);
           d_agent=norm(cen_initial(k,:)-agent_pos)-rad_min-rad_min;
           if norm(cen_initial(k,:)-agent_pos)<=2*rad_max
               phi=phi_kl(s_a(k), s_a_values_k(l), t(iter), tc(k), b);
               hat_m=  (1/(norm(cen_initial(k,:)-agent_pos))-1/(2*rad_max))*hat_h_2 * (1 / d_agent^2) * ((cen_initial(k,:) - agent_pos) / d_agent);
               N_agent=null(hat_m);
               hat_v=hat_h_3*N_agent(:,1)';
               d_cen_3=d_cen_3+phi*(hat_m+hat_v);
           end
        end
        % Total d_cen calculation as in (12)
        d_cen = gamma + d_cen_3;
        % euler update to get update center of STT
        if norm(d_cen) > 0
            cen(k,:) = cen_initial(k,:) + delt * (d_cen );
        end
    end

for k=1:n_a
    %calculating d_1(k) for (17)
    sum_dist=10^4;
    if n_a>1
        agent_pos_radius=[cen([1:k-1,k+1:end],:)];
         s_a_values_rad=[s_a([1:k-1,k+1:end],:)];
        for l=1:size(agent_pos_radius,1)
            phi = phi_kl(s_a(k), s_a_values_rad(l), t(iter), tc(k), b);
            d_2k=norm(cen(k,:)-agent_pos_radius(l,:));
            if d_2k< 2*rad_max
                % sum_dist=min(sum_dist,rad_min(k)+(d_2k-(rad_min(k)+rad_min(l)))*(s_a_values_rad(l)/(s_a(k)+s_a_values_rad(l))));
               w= (1-(d_2k-2*rad_min)/(2*rad_max-2*rad_min))*(1-phi)+(d_2k-2*rad_min)/(2*rad_max-2*rad_min)*(1/2);
               sum_dist=smooth_min(sum_dist,rad_min+(d_2k-2*rad_min)*w);
            end
         end
    end
    rad(k)=min(rad_max,sum_dist); % final radius update
 % upper and lower funnel boundary
    gamma_L=-4*exp(-0.06*t(iter))*ones(1,3);
    gamma_U=4*exp(-0.06*t(iter))*ones(1,3);    
    gamma_d=(gamma_U-gamma_L);
    d_gamma_dot=-4*0.06*2*exp(-0.06*t(iter))*ones(1,3);

    %% Control update
    u = control_3d(cen_initial(k,:), rad_past(k,:), robot(k,:),velo(k,:),gamma_L,gamma_U,gamma_d,d_gamma_dot); % External function required
    % updating states
    robot(k,:) = robot(k,:) + delt * velo(k,:);
    velo(k,:)=velo(k,:)+delt*u;
    % storing trajectories
    robot_history{k} = [robot_history{k}; robot(k,:)];
    tube_history{k}  = [tube_history{k};  cen(k,:)];
   
end
cen_initial=cen;
rad_past=rad;
 if snapIndex <= length(snapTimes) && abs(t(iter) - snapTimes(snapIndex)) < 1e-6  
       subplot(2,2,snapIndex);
       hold on; grid on; axis equal;
       view(3)   % set 3D view
       axis([-3 10 -3 10 -3 10]);
       title(['t = ' num2str(t(iter),'%.1fs')]);
       xlabel('$X(m)$'); ylabel('$Y(m)$');zlabel('$Z(m)$');
       for i = 1:n_a
             % Plot goal as translucent sphere
             [sx, sy, sz] = sphere(20);
              surf(eta(i,1) + sx, eta(i,2) + sy, eta(i,3) + sz, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'g');
       end
        % --- Draw obstacles --
        for i = 1:n_a
            k_snap=i;
           if i==1 || i==6
                plot3(robot_history{i}(:,1), robot_history{i}(:,2), robot_history{i}(:,3),'Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
                robot_plot(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'o','Color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
                tube_plot(i) = plot3(cen(i,1), cen(i,2), cen(i,3), 'o','Color',[0.9290, 0.6940, 0.1250]);
                trajectory(i) = plot3(robot(i,1), robot(i,2), robot(i,3), '-','Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
                [sx, sy, sz] = sphere(20);   % finer resolution
                hSafe(i) = surf(nan(size(sx)), nan(size(sy)), nan(size(sz)), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                'FaceColor', [0.9290, 0.6940, 0.1250], ...
                'Tag', 'safe_sphere');
                set(hSafe(i), 'XData', cen(i,1) + rad(i)*sx, ...
                       'YData', cen(i,2) + rad(i)*sy, ...
                       'ZData', cen(i,3) + rad(i)*sz);
                      else
                plot3(robot_history{i}(:,1), robot_history{i}(:,2), robot_history{i}(:,3),'b', 'LineWidth', 1.5);
                robot_plot(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'ko', 'MarkerFaceColor', 'k');
                tube_plot(i) = plot3(cen(i,1), cen(i,2), cen(i,3), 'bo');
                trajectory(i) = plot3(robot(i,1), robot(i,2), robot(i,3), 'b-', 'LineWidth', 1.5);
                [sx, sy, sz] = sphere(20);   % finer resolution
                hSafe(i) = surf(nan(size(sx)), nan(size(sy)), nan(size(sz)), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                'FaceColor', [0.3010 0.7450 0.9330], ...
                'Tag', 'safe_sphere');
                set(hSafe(i), 'XData', cen(i,1) + rad(i)*sx, ...
                       'YData', cen(i,2) + rad(i)*sy, ...
                       'ZData', cen(i,3) + rad(i)*sz);
           end
        end
        snapIndex = snapIndex + 1;
    end

    if norm(cen - eta) < 0.2
        disp(['Goal reached at t = ', num2str(t(iter))]);
        break;
    end
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

