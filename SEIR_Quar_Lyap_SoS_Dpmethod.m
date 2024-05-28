clear 
syms w y z 
vars = [w, y, z]; 

R = 0;
% Try params
while R < 1 %|| R > 1.1
delta = rand;
beta = rand;
eta = rand;
gamma = rand;
R = (beta*eta)/((delta + eta)*(delta + gamma));
end

%delta = 0.2;
%beta = 0.8;  
%gamma = 0.05;
%eta = 0.5;

% delta = 0.3;
% beta = 0.8;
% gamma = 0.015;
% eta = 0.3;

n = 1;
delta = 0.3; %0.001:0.01:1 %0.3
beta = 0.8; %0.001:0.01:1 %0.8
gamma = 0.015; %0.0001:0.001:1 %0.015
for eta = 0.3; %0.001:0.01:1 %0.3;

R = (beta*eta)/((delta + eta)*(delta + gamma));
if R < 1
    continue
end

R = 0;
% Try params
while R < 1 %|| R > 1.1
delta = rand;
beta = rand;
eta = rand;
gamma = rand;
R = (beta*eta)/((delta + eta)*(delta + gamma));
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

% % Shift coordinates
w_dot = subs(w_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]); %AP
y_dot = subs(y_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]); %AP
z_dot = subs(z_dot,[w,y,z],[w + equil_w,y + equil_y,z + equil_z]); %AP

% SOS
prog = sosprogram(vars);

orderV = 2;
%[prog,V] = sossosvar(prog,monomials(vars,0:2),'wscoeff'); %,vars.^2]);
%[prog,V] = sospolyvar(prog,monomials(vars,[0,2:orderV]),'wscoeff'); %AP %deffo dont need this
[prog,V] = sospolyvar(prog,monomials(vars,2:orderV),'wscoeff');

% [prog,a1] = sospolyvar(prog,1);
% [prog,a2] = sospolyvar(prog,1);
% [prog,a3] = sospolyvar(prog,1);
% [prog,a4] = sospolyvar(prog,1);
% [prog,a5] = sospolyvar(prog,1); 
% V = a1*w^2 + a2*w*y + a3*w*z + a4*y^2 + a5*y*z + z^2;
phi = 0;
for i = 1:length(vars)
    constr = sym(0.01);
    for j = 1:orderV/2
        [prog,eps(i,j)] = sossosvar(prog,1,'wscoeff');
        phi = phi+eps(i,j)*vars(i)^(2*j);
        constr = constr-eps(i,j);
    end
    prog = sosineq(prog,-constr);
end
% V = (x - x*)'Q(x - x*)
% V = subs(V,[w,y,z],[w - equil_w,y - equil_y,z - equil_z]); %AP
prog = sosineq(prog,V-phi);

% Region D
%r = (equil.w(1)-equil.w(2))^2 + (equil.y(1)-equil.y(2))^2 + (equil.z(1)-equil.z(2))^2; %AP
%D = w^2 + y^2 + z^2 - (r - 0.05);
%bound = -1; 
%D1 = -w + bound*equil_w; D2 = -y + bound*equil_y; D3 = -z + bound*equil_z; %leq 0
%[prog,p1] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP
%[prog,p2] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP
%[prog,p3] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP

% derV <= 0 for x > 0
derV = -(diff(V,y)*y_dot+diff(V,z)*z_dot+diff(V,w)*w_dot);% + D1*p1 + D2*p2 + D3*p3;
%derV = subs(derV,[w,y,z],[(w+equil_w)^2-equil_w,(y+equil_y)^2-equil_y,(z+equil_z)^2-equil_z]);
derV = subs(derV,[w,y,z],[(w)^2-equil_w,(y)^2-equil_y,(z)^2-equil_z]);
prog = sosineq(prog,derV);

% Fix first coefficient
%prog = soseq(prog,coeff_2 - 1);

% Set Lyap to equal zero at minimum
% prog = soseq(prog,coeff_2);
% prog = soseq(prog,coeff_3);
%prog = soseq(prog,coeff_4);

prog = soseq(prog,coeff_6 - 1);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL = sosgetsol(prog,V)

%SOL2 = SOL/sosgetsol(prog,coeff_1)

if prog.solinfo.info.pinf == 0
    temp1 = flip(coeffs(SOL));
    a1(n) = temp1(1);
    a2(n) = temp1(2);
    a3(n) = temp1(3);
    a4(n) = temp1(4);
    a5(n) = temp1(5);
    a6(n) = temp1(6);
    %a7(n) = temp1(7);
    delta_value(n) = delta;
    beta_value(n) = beta;
    gamma_value(n) = gamma;
    eta_value(n) = eta;
    n = n + 1;
end

end
return
% p1 = polyfit(eta_value,eta_value.*a1,1);
% p2 = polyfit(eta_value,eta_value.*a2,1);
% p3 = polyfit(eta_value,eta_value.*a3,1);
% p4 = polyfit(eta_value,eta_value.*a4,1);
% p5 = polyfit(eta_value,eta_value.*a5,1);
% p6 = polyfit(eta_value,eta_value.*a6,1);
% %p7 = polyfit(eta_value,eta_value.*a7,1);

%subs(SOL,[w,y,z],[equil_w,equil_y,equil_z])
%subs(sosgetsol(prog,phi),[w,y,z],[equil_w,equil_y,equil_z])


p1 = polyfit(eta_value,a1,1);
p2 = polyfit(eta_value,a2,1);
p3 = polyfit(eta_value,a3,1);
p4 = polyfit(eta_value,a4,1);
p5 = polyfit(eta_value,a5,1);
p6 = polyfit(eta_value,a6,1);
%p7 = polyfit(eta_value,a7,1);

figure
plot(eta_value,a1)
figure
plot(eta_value,a2)
figure
plot(eta_value,a3)
figure
plot(eta_value,a4)
figure
plot(eta_value,a5)
figure
plot(eta_value,a6)
%figure
%plot(eta_value,a7)

return
%SOL2 = SOL;

SOL = subs(SOL,[w,y,z],[w+equil_w,y+equil_y,z+equil_z])

SOL_wy = subs(SOL,z,equil_z);
SOL_wz = subs(SOL,y,equil_y);
SOL_yz = subs(SOL,w,equil_w);

rad = vpa(r^(0.5),3)

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


