%Solves the 5th order polynomial to find the lagrangian points
%Uses set of eqns on pg 972 of Vallado

function [L1x,L2x,L3x] = findLagrangianPoints2(mu)

f1 = @(x) x^5 + (3-mu)*x^4 + (3-2*mu)*x^3 - mu*x^2 - 2*mu*x - mu;
f2 = @(x) x^5 - (3-mu)*x^4 + (3-2*mu)*x^3 - mu*x^2 + 2*mu*x - mu;
f3 = @(x) x^5 + (2+mu)*x^4 + (1+2*mu)*x^3 - (1-mu)*x^2 - 2*(1-mu)*x - (1-mu);

L1x = fsolve(f1,0.5);
L2x = fsolve(f2,1.24);
L3x = fsolve(f3,-1.1);

fprintf('At L1, the computed value is %f\n',f1(L1x))
fprintf('At L2, the computed value is %f\n',f2(L2x))
fprintf('At L3, the computed value is %f\n',f3(L3x))

end