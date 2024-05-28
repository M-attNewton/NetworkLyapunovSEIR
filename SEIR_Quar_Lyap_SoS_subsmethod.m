clear 
syms w y z 
vars = [w, y, z]; 

% delta = 0.3;
% beta = 0.8;
% gamma = 0.015;
% eta = 0.3;

% Caused failure
% delta = 0.0082;
% beta = 0.462;
% gamma = 0.0381;
% eta = 0.6948;
% 
tic
n = 1;
k = 1;
for beta = 0:0.01:1
%R = 20;
% Try params
delta = 0.1;
%beta = 0.5;
eta = 0.5;
gamma = 0.1;

% while R > 1 %|| R > 1.1
%     %delta = rand;
%     beta = rand;
%     eta = rand;
%     gamma = rand;
%     R = (beta*eta)/((delta + eta)*(delta + gamma));
% end
%    delta = 0.1;
%    beta = 0.9825;
%    eta = 0.9467;
%    gamma = 0.1795;
R = (beta*eta)/((delta + eta)*(delta + gamma))

y_dot = -delta*y - beta*y*w + delta;
z_dot = -(delta + eta)*z + beta*y*w;
w_dot = -(delta + gamma)*w + eta*z;

% Find equilibrium
equil = solve([y_dot == 0; z_dot == 0; w_dot == 0], [w,y,z]);
for j = 1:length(equil.y)
    jacob = jacobian ([y_dot;z_dot;w_dot], [y,z,w]);
    jacobsubs = subs(jacob,[y,z,w],[equil.y(j),equil.z(j),equil.w(j)]);
    if all(eig(vpa(jacobsubs,10)) <= 0)
        equil_y = vpa(equil.y(j),10);
        equil_z = vpa(equil.z(j),10);
        equil_w = vpa(equil.w(j),10); % this was just changed from 4 to 10
    end
end

% Shift coordinates
% w_dot = 1/equil_w*subs(w_dot,[w,y,z],[(w+1)*equil_w,(y+1)*equil_y,(z+1)*equil_z]); %AP
% y_dot = 1/equil_y*subs(y_dot,[w,y,z],[(w+1)*equil_w,(y+1)*equil_y,(z+1)*equil_z]); %AP
% z_dot = 1/equil_z*subs(z_dot,[w,y,z],[(w+1)*equil_w,(y+1)*equil_y,(z+1)*equil_z]); %AP
% equil_w = 1;
% equil_y = 1;
% equil_z = 1;

% w_dot = subs(w_dot,[w,y,z],[w+equil_w,y+equil_y,z+equil_z]); %AP
% y_dot = subs(y_dot,[w,y,z],[w+equil_w,y+equil_y,z+equil_z]); %AP
% z_dot = subs(z_dot,[w,y,z],[w+equil_w,y+equil_y,z+equil_z]); %AP

w_dot = expand(w_dot - subs(w_dot,[w,y,z],[0,0,0]));
y_dot = expand(y_dot - subs(y_dot,[w,y,z],[0,0,0]));
z_dot = expand(z_dot - subs(z_dot,[w,y,z],[0,0,0]));
A = subs(jacobian([w_dot;y_dot;z_dot],[w,y,z]),[y,z,w],[0,0,0]);
% SOS
prog = sosprogram(vars);

orderV = 4;
%[prog,V] = sossosvar(prog,monomials(vars,0:2),'wscoeff'); %,vars.^2]);
%[prog,V] = sospolyvar(prog,monomials(vars,0:orderV),'wscoeff'); %AP
% Z = [     w^2
%      w*y
%      y^2
%      w*z
%      y*z
%      z^2
% %    w^3
%    w^2*y
%    w*y^2
%      y^3
% %  w^2*z
% %   w*y*z
%    y^2*z
%  %  w*z^2
% %   y*z^2
% %     z^3
% %     w^4
% %   w^3*y
%  w^2*y^2
% %   w*y^3
%      y^4
% %   w^3*z
% % w^2*y*z
%  w*y^2*z
%   y^3*z
% % w^2*z^2
% % w*y*z^2
%  y^2*z^2
% %   w*z^3
%  %  y*z^3
%      %z^4
%      ];
Z = monomials(vars,2:orderV);
[prog,V] = sospolyvar(prog,Z,'wscoeff'); %AP

%  summag = 0;
%  for i = 1:length(Z);
%      [prog,mag(i)] = sossosvar(prog,1);    
%      prog = sosineq(prog,['coeff_',num2str(i)]+mag(i));
%      prog = sosineq(prog,mag(i)-['coeff_',num2str(i)]);     
%      summag = summag+mag(i);
%  end
%  prog = sossetobj(prog,summag);

phi = 0;
for i = 1:length(vars)
    constr = sym(0.1);
    for j = 1:orderV/2
        [prog,eps(i,j)] = sossosvar(prog,1);
        phi = phi+eps(i,j)*vars(i)^(2*j);
        constr = constr-eps(i,j);
    end
    prog = sosineq(prog,-constr);
end

prog = sosineq(prog,V-phi);

% derV <= 0 for x > 0
derV = -(diff(V,y)*y_dot+diff(V,z)*z_dot+diff(V,w)*w_dot);
derV = subs(derV,[w,y,z],[(w)^2-equil_w,(y)^2-equil_y,(z)^2-equil_z]);
%derV = subs(derV,[w,y,z],[(w)^2-equil_w,(y)^2-equil_y,(z)^2-equil_z]);
%derV = subs(derV,[w,y,z],[(w+equil_w)^2-equil_w,(y+equil_y)^2-equil_y,(z+equil_z)^2-equil_z]);
prog = sosineq(prog,derV,'sparse');

prog = soseq(prog,coeff_31 - 1);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL = sosgetsol(prog,V)
SOL_coeff = coeffs(SOL);
SOL2 = vpa(SOL/max(SOL_coeff),5)
if prog.solinfo.info.pinf == 1
    test_sol(n) = SOL;
    test_delta(n) = delta;
    test_beta(n) = beta;
    test_gamma(n) = gamma;
    test_eta(n) = eta;
    test_R(n) = R;
    test_equil_w(n) = equil_w;
    test_equil_y(n) = equil_y;
    test_equil_z(n) = equil_z;
    n = n + 1;
    %return
end

if prog.solinfo.info.pinf == 0
    temp1(k,:) = flip(coeffs(SOL));
%     a1(k) = temp1(1);
%     a2(k) = temp1(2);
%     a3(k) = temp1(3);
%     a4(k) = temp1(4);
%     a5(k) = temp1(5);
%     a6(k) = temp1(6);
%     %a7(n) = temp1(7);
    delta_value(k) = delta;
    beta_value(k) = beta;
    gamma_value(k) = gamma;
    eta_value(k) = eta;
    k = k + 1;
end


end
time = toc

for i = 1:size(temp1,2)
   
    figure
    plot(beta_value,temp1(:,i))
    
    %p(i,:) polyfit
end

return

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


% Region D
%r = (equil.w(1)-equil.w(2))^2 + (equil.y(1)-equil.y(2))^2 + (equil.z(1)-equil.z(2))^2; %AP
%D = w^2 + y^2 + z^2 - (r - 0.05);
%bound = -1; 
%D1 = -w + bound*equil_w; D2 = -y + bound*equil_y; D3 = -z + bound*equil_z; %leq 0
%[prog,p1] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP
%[prog,p2] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP
%[prog,p3] = sossosvar(prog,monomials(vars,0:orderV/2-1),'wscoeff'); %AP

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


