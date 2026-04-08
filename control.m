function u = control(tube,rad,robot)
e = norm(robot - tube)/rad;
k = 20;
u = -k * (robot - tube) * log((1+e)/(1-e));
u = real(u);
end

