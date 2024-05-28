% Start with easy handcraft two system network

syms w1 y1 z1 W1 Y1 Z1 w2 y2 z2 W2 Y2 Z2 
vars = [w1 y1 z1 W1 Y1 Z1 w2 y2 z2 W2 Y2 Z2 ]; 

R = 0;
% Try params, to start make both systems have same values
while R < 1
delta = rand;
beta = rand;
eta = rand;
gamma = rand;
R = (beta*eta)/((delta + eta)*(delta + gamma))
alpha = rand;
end

% System 1
y1_dot = -delta*y1 - beta*y1*w1 + delta - alpha*y1*w2;
z1_dot = -(delta + eta)*z1 + beta*y1*w1 + alpha*y1*w2;
w1_dot = -(delta + gamma)*w1 + eta*z1;
Y1_dot = (y1^-1)*y1_dot;
Z1_dot = (z1^-1)*z1_dot;
W1_dot = (w1^-1)*w1_dot;

% System 2
y2_dot = -delta*y2 - beta*y2*w2 + delta - alpha*y2*w1;
z2_dot = -(delta + eta)*z2 + beta*y2*w2 + alpha*y2*w1;
w2_dot = -(delta + gamma)*w2 + eta*z2;
Y2_dot = (y2^-1)*y2_dot;
Z2_dot = (z2^-1)*z2_dot;
W2_dot = (w2^-1)*w2_dot;

% Find equilibrium
equil = solve([y1_dot == 0; z1_dot == 0; w1_dot == 0; y2_dot == 0; z2_dot == 0; w2_dot == 0],...
    [w1,y1,z1,w2,y2,z2]);
for j = 1:length(equil.y1)
    jacob = jacobian ([y1_dot;z1_dot;w1_dot;y2_dot;z2_dot;w2_dot], [y1,z1,w1,y2,z2,w2]);
    jacobsubs = subs(jacob,[y1,z1,w1,y2,z2,w2], ...
        [equil.y1(j),equil.z1(j),equil.w1(j),equil.y2(j),equil.z2(j),equil.w2(j)]);
    if all(eig(vpa(jacobsubs,4)) <= 0)
        equil_y1 = vpa(equil.y1(j),4);
        equil_z1 = vpa(equil.z1(j),4);
        equil_w1 = vpa(equil.w1(j),4);
        equil_y2 = vpa(equil.y2(j),4);
        equil_z2 = vpa(equil.z2(j),4);
        equil_w2 = vpa(equil.w2(j),4);
    end
end

% SOS
prog = sosprogram(vars);

mono = monomials(vars,0:1);
[prog,V] = sospolyvar(prog,mono,'wscoeff');

% mono2 = monomials([w,y,z],0:2);
% [prog,q] = sospolyvar(prog, mono2,'wscoeff');

% SoS constraints

% 1) V >= 0 
%prog = sosineq(prog,V - 0.1);

% 2) V_dot
constraint2 = -(diff(V,w1)*w1_dot + diff(V,y1)*y1_dot + diff(V,z1)*z1_dot + ...
                diff(V,W1)*W1_dot + diff(V,Y1)*Y1_dot + diff(V,Z1)*Z1_dot) ...
              -(diff(V,w2)*w2_dot + diff(V,y2)*y2_dot + diff(V,z2)*z2_dot + ...
                diff(V,W2)*W2_dot + diff(V,Y2)*Y2_dot + diff(V,Z2)*Z2_dot);
constraint2 = expand(constraint2*w1*y1*z1*w2*y2*z2);
constraint2 = subs(constraint2, [w1, y1, z1, w2, y2, z2], [w1^2, y1^2, z2^2, w2^2, y2^2, z2^2]);
constraint2  = constraint2*(y1^2+z1^2+w1^2+1)^2 * (y2^2+z2^2+w2^2+1)^2;
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

% % Shift coordinates
% w_dot = subs(w_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% y_dot = subs(y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% z_dot = subs(z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% 
% W_dot = subs(W_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% Y_dot = subs(Y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);
% Z_dot = subs(Z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]);



