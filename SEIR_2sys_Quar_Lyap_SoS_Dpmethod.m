clear 
syms w1 y1 z1 w2 y2 z2
vars = [w1, y1, z1, w2, y2, z2]; 

R1 = 22; R2 = 0;
% Try params
while R1 > 0.1 || R2 < 10
delta1 = rand;
beta1 = rand;
eta1 = rand;
gamma1 = rand;
R1 = (beta1*eta1)/((delta1 + eta1)*(delta1 + gamma1))

delta2 = rand;
beta2 = rand;
eta2 = rand;
gamma2 = rand;
R2 = (beta2*eta2)/((delta2 + eta2)*(delta2 + gamma2))

alpha1 = 0.5; %rand;
alpha2 = 0.5; %rand;
end

y1_dot = -delta1*y1 - beta1*y1*w1 + delta1 - alpha1*y1*w2;
z1_dot = -(delta1 + eta1)*z1 + beta1*y1*w1 + alpha1*y1*w2;
w1_dot = -(delta1 + gamma1)*w1 + eta1*z1;

y2_dot = -delta2*y2 - beta2*y2*w2 + delta2 - alpha2*y2*w1;
z2_dot = -(delta2 + eta2)*z2 + beta2*y2*w2 + alpha2*y2*w1;
w2_dot = -(delta2 + gamma2)*w2 + eta2*z2;

% Find equilibrium
equil = solve([y1_dot == 0; z1_dot == 0; w1_dot == 0; y2_dot == 0; z2_dot == 0; w2_dot == 0], [w1,y1,z1,w2,y2,z2]);
for j = 1:length(equil.y1)
    jacob = jacobian ([y1_dot;z1_dot;w1_dot;y2_dot;z2_dot;w2_dot], [y1,z1,w1,y2,z2,w2]);
    jacobsubs = subs(jacob,[y1,z1,w1,y2,z2,w2],[equil.y1(j),equil.z1(j),equil.w1(j),equil.y2(j),equil.z2(j),equil.w2(j)]);
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

orderV = 4;
%[prog,V] = sossosvar(prog,monomials(vars,0:2),'wscoeff'); %,vars.^2]);
[prog,V] = sospolyvar(prog,monomials(vars,0:4),'wscoeff');
phi = 0;
for i = 1:length(vars)
    constr = sym(0.01);
    for j = 1:orderV/2
        [prog,eps(i,j)] = sossosvar(prog,1);
        phi = phi+eps(i,j)*vars(i)^(2*j);
        constr = constr-eps(i,j);
    end
    prog = sosineq(prog,-constr);
end

% V = (x - x*)'Q(x - x*)
V = subs(V,[w1,y1,z1,w2,y2,z2],[w1 - equil_w1,y1 - equil_y1,z1 - equil_z1,w2 - equil_w2,y2 - equil_y2,z2 - equil_z2]);
prog = sosineq(prog,V-phi);

% Region D
%r = (equil.w(1)-equil.w(2))^2 + (equil.y(1)-equil.y(2))^2 + (equil.z(1)-equil.z(2))^2;
%D = w^2 + y^2 + z^2 - (r - 0.05);
bound = 0; 
D1 = -w1 + bound; D2 = -y1 + bound; D3 = -z1 + bound;
D4 = -w2 + bound; D5 = -y2 + bound; D6 = -z2 + bound;
[prog,p1] = sossosvar(prog,monomials(vars,0:1),'wscoeff');
[prog,p2] = sossosvar(prog,monomials(vars,0:1),'wscoeff');
[prog,p3] = sossosvar(prog,monomials(vars,0:1),'wscoeff');
[prog,p4] = sossosvar(prog,monomials(vars,0:1),'wscoeff');
[prog,p5] = sossosvar(prog,monomials(vars,0:1),'wscoeff');
[prog,p6] = sossosvar(prog,monomials(vars,0:1),'wscoeff');

% derV <= 0 for x > 0
derV = -(diff(V,y1)*y1_dot+diff(V,z1)*z1_dot+diff(V,w1)*w1_dot+....
         diff(V,y2)*y2_dot+diff(V,z2)*z2_dot+diff(V,w2)*w2_dot) ...
       + D1*p1 + D2*p2 + D3*p3 + D4*p4 + D5*p5 + D6*p6;
%derV = subs(derV,[w,y,z],[w^2,y^2,z^2]);
prog = sosineq(prog,derV);

prog = soseq(prog,coeff_1 - 1);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL = sosgetsol(prog,V)

SOL2 = SOL;

%SOL = subs(SOL,[w,y,z],[w+equil_w,y+equil_y,z+equil_z])

SOL_wy = subs(SOL,z,equil_z);
SOL_wz = subs(SOL,y,equil_y);
SOL_yz = subs(SOL,w,equil_w);

%rad = vpa(r^(0.5),3)

figure
subplot(2,3,1)
fsurf(SOL_wy,[-0.05 1])
subplot(2,3,2)
fsurf(SOL_wz,[-0.05 1])
subplot(2,3,3)
fsurf(SOL_yz,[-0.05 1])

subplot(2,3,4)
fcontour(SOL_wy,[-0.05 1])
subplot(2,3,5)
fcontour(SOL_wz,[-0.05 1])
subplot(2,3,6)
fcontour(SOL_yz,[-0.05 1])

timeDer = -(diff(SOL,y)*y_dot+diff(SOL,z)*z_dot+diff(SOL,w)*w_dot);
timemin = [diff(timeDer,y); diff(timeDer,z); diff(timeDer,w)];
timeminsol = subs(timemin,[y,z,w],[equil_y,equil_z,equil_w]);
timejacob = jacobian(jacobian(timeDer,[y,z,w]),[y,z,w]);
timejacobsol = subs(timejacob,[y,z,w],[equil_y,equil_z,equil_w]);
hessian_eigenvalues  = eig(timejacobsol)
equilibrium = [equil_w,equil_y,equil_z,1-equil_w-equil_y-equil_z]

timeDer_wy = subs(timeDer,z,equil_z);
timeDer_wz = subs(timeDer,y,equil_y);
timeDer_yz = subs(timeDer,w,equil_w);

%plot time derivative function
figure
subplot(2,3,1)
fsurf(timeDer_wy,[-0.05 1])
subplot(2,3,2)
fsurf(timeDer_wz,[-0.05 1])
subplot(2,3,3)
fsurf(timeDer_yz,[-0.05 1])

subplot(2,3,4)
fcontour(timeDer_wy,[-0.05 1])
subplot(2,3,5)
fcontour(timeDer_wz,[-0.05 1])
subplot(2,3,6)
fcontour(timeDer_yz,[-0.05 1])









%%JUNK

%ineq1 = w+equil_w; %geq
%ineq2 = y+equil_y;
%ineq3 = z+equil_z;
%[prog,s1] = sossosvar(prog,monomials(vars,orderV/2));
%[prog,s2] = sossosvar(prog,monomials(vars,orderV/2));
%[prog,s3] = sossosvar(prog,monomials(vars,orderV/2));
%prog = sosineq(prog,derV-ineq1*s1-ineq2*s2-ineq3*s3);

%return
%a1*(y*(log(y)-log(equil_y)-1)+equil_y)+a2*(w*(log(w)-log(equil_w)-1)+equil_w)+a3*(z*(log(z)-log(equil_z)-1)+equil_z);

%derV2 = expand(derV1)*y*z*w;
%
%derV = derV3*(y^2+z^2+w^2+1)^2;

%[prog,a1] = sossosvar(prog,1);
%[prog,a2] = sossosvar(prog,1);
%[prog,a3] = sossosvar(prog,1);
%a1 = 1;
%V = a1*(y-equil_y*(log(y)))+a2*(w-equil_w*(log(w)))+a3*(z-equil_z*(log(z)));


