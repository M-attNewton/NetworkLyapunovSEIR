clear 

syms y u
p = 1/y + 1/u + y*u - 3;
r = u*y*(u^2 + y^2 + 1);
q = subs(p*r,[u,y],[u^2,y^2]);

[a,b,c] = findsos(q);

for i = 1:size(a,1)
    for j = 1:size(a,2)
        if vpa(a(i,j)) < 0.001
            a(i,j) = 0;
        end
    end
end

test1 = b.'*a*b;
vpa(test1,2)