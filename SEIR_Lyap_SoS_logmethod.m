syms w y z W Y Z %y_star z_star w_star delta beta eta gamma
vars = [w, y, z, W, Y, Z]; 

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
Y_dot = (y^-1)*y_dot;
Z_dot = (z^-1)*z_dot;
W_dot = (w^-1)*w_dot;

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

% % Shift coordinates
% w_dot = subs(w_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% y_dot = subs(y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% z_dot = subs(z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% 
% W_dot = subs(W_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% Y_dot = subs(Y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% Z_dot = subs(Z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);


% SOS
prog = sosprogram(vars);

mono = monomials(vars,0:1);
[prog,V] = sospolyvar(prog,mono,'wscoeff');

mono2 = monomials([w,y,z],0:2);
[prog,q] = sospolyvar(prog, mono2,'wscoeff');

% mono2 = monomials(vars,0:2);
% [prog,poly] = sospolyvar(prog, mono2, 'wscoeff') ;

% SoS constraints

% 1) V >= 0 
%prog = sosineq(prog,V - 0.1);

% 2) V_dot
constraint2 = -(diff(V,w)*w_dot + diff(V,y)*y_dot + diff(V,z)*z_dot + ...
                diff(V,W)*W_dot + diff(V,Y)*Y_dot + diff(V,Z)*Z_dot);
constraint2 = expand(constraint2*w*y*z);
constraint2 = subs(constraint2, [w, y, z], [w^2, y^2, z^2]);
constraint2  = constraint2*(y^2+z^2+w^2+1)^2;
%constraint2 = constraint2*q^2;
prog = sosineq(prog,expand(constraint2));

% 3) d^2V/dx^2 >=0 at equilibrium
%constraint3 = jacobian(jacobian(V,[vars]),[vars]);
%prog = sosmatrixineq(prog,constraint3,'quadraticMineq');

prog = soseq(prog,coeff_2 - 1);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

SOL = sosgetsol(prog,V);

sol = vpa(subs(SOL,[W,Y,Z],[log(w), log(y), log(z)]),3)



