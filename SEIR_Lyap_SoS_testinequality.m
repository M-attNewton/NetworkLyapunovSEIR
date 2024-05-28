syms u y

vars = [u,y];

p_0 = u^2*y^2 + 1/u^2 + 1/y^2 - 3;

p = u^2*y^2*(p_0);

prog = sosprogram(vars);

mono = monomials(vars,0:2);

solver_opt.solver = 'sedumi';
[prog,q1] = sospolyvar(prog,mono,'wscoeff');

prog = sosineq(prog,q1*p);

prog = sosineq(prog,q1 - 0.01);

prog = soseq(prog,coeff_4 - 1);

prog = sossolve(prog);

SOL = sosgetsol(prog,q1)

findsos(p*SOL)

% Test to see if non-sos can become sos with multiple sos
% syms x
% test1 = x^2 + x;
% vars1 = x;
% prog1 = sosprogram(vars1);
% mono1 = monomials(vars1,0:10);
% solver_opt.solver = 'sedumi';
% [prog1,q11] = sospolyvar(prog1,mono1,'wscoeff');
% prog1 = sosineq(prog1,q11*test1 - 0.1);
% prog1 = sossolve(prog1);
% SOL1 = sosgetsol(prog1,q11)

%q = u^2 + y^2 + 1;

%[Q,Z,f] = findsos(expand(q*p))

%findsos(u*y*(u^2 + y^2 + 1))

%sol = solve([diff(a,u)==0, diff(a,y) == 0, jacobian(jacobian(a,[u,y]),[u,y]) <= 0])

%q = u*y*(u^2 + y^2 + 1);

%p = expand(p_0*q);

%[Q,Z] = findsos(p)
% 
% b = (2*u^2*y^2 - 6*u*y + 2*u + 2*y + u^2 + y^2 + 10);
% 
% a = b/(2*u^2*y^2 - 6*u*y + 2*u + 2*y);
% 
% c = (1/2)*(u^2 + y^2 + 10) 
