function u = control(tube, rad, robot,velo,rho_L,rho_U,rho_d,rho_d_dot)
k_1 =4; % controller gain in first satge
k_2=12;% controller gain in second satge
e_1 = norm(robot - tube)/rad;
r_2 = (-k_1 * (robot - tube) * log((1+e_1)/(1-e_1))); % calculation of refrence signal r_2
% r_2 = real(r_2);
for i=1:length(robot)
    x=velo(i)-r_2(i);
    e_2=(x-0.5*(rho_U(i)+rho_L(i)))/(0.5*(rho_d(i)));
    eps=log((1+e_2)/(1-e_2));
    zeta=4/(rho_d(i)*(1-e_2^2));
    u(i)=-(k_2*zeta*eps-0.5*rho_d_dot(i)*e_2);% final control ip
end
% u = real(u);
end

