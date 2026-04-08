function u = control(tube,rad,robot)
% u = 10*(tube-robot);

a = 7;

e = norm(robot - tube)/rad;
ub = 20;
%u = -ub*tanh(a*e).*(1-exp(-a*a*e.^2)) * (robot - tube) / norm(robot - tube);
u = -ub * (robot - tube) * log((1+e)/(1-e));
u = real(u);

end

