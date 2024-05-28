syms w y z %y_star z_star w_star delta beta eta gamma
vars = [w, y, z]; 

R = 0;
% Try params
while R < 1
delta = rand;
beta = rand;
eta = rand;
gamma = rand;
R = (beta*eta)/((delta + eta)*(delta + gamma))
end

y_dot = -delta*y - beta*y*w + delta;
z_dot = -(delta + eta)*z + beta*y*w;
w_dot = -(delta + gamma)*w + eta*z;

% Find equilibrium
equil = solve([y_dot == 0; z_dot == 0; w_dot == 0], [w,y,z]);
for j = 1:length(equil.y)
    jacob = jacobian ([y_dot;z_dot;w_dot], [y,z,w]);
    jacobsubs = subs(jacob,[y,z,w],[equil.y(j),equil.z(j),equil.w(j)]);
    if all(eig(vpa(jacobsubs,4)) <= 0)
        equil_y = vpa(equil.y(j),4);
        equil_z = vpa(equil.z(j),4);
        equil_w = vpa(equil.w(j),4);
    end
end

% Shift coordinates
%w_dot = subs(w_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
%y_dot = subs(y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
%z_dot = subs(z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);

% SOS
prog = sosprogram(vars);

%[prog,a1] = sossosvar(prog,1);
[prog,a2] = sossosvar(prog,1);
[prog,a3] = sossosvar(prog,1);
a1 = 1;
V = a1*(y-equil_y*(log(y)))+a2*(w-equil_w*(log(w)))+a3*(z-equil_z*(log(z)));
%V = a1*(y*(log(y)-log(equil_y)-1)+equil_y)+a2*(w*(log(w)-log(equil_w)-1)+equil_w)+a3*(z*(log(z)-log(equil_z)-1)+equil_z);
derV1 = -(diff(V,y)*y_dot+diff(V,z)*z_dot+diff(V,w)*w_dot);
derV2 = expand(derV1)*y*z*w;
derV3 = subs(derV2,[w,y,z],[w^2,y^2,z^2]);
derV = derV3*(y^2+z^2+w^2+1)^2;
prog = sosineq(prog,derV);
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL = sosgetsol(prog,V)


