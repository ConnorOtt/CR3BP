%Solves the 5th order polynomial to find the lagrangian points
%Uses equations on http://www.geom.uiuc.edu/~megraw/MATH1/lib.html

function [L1x,L2x,L3x] = findLagrangianPoints3(mu)

f1 = @(x) -((1-mu)*(1-mu+x)^2) + mu*(mu+x)^2 + x*(1-mu+x)^2 *(mu+x)^2;
f2 = @(x) -((1-mu)*(1-mu+x)^2) - mu*(mu+x)^2 + x*(1-mu+x)^2 *(mu+x)^2;
f3 = @(x) ((1-mu)*(1-mu+x)^2) + mu*(mu+x)^2 + x*(1-mu+x)^2 *(mu+x)^2;

L1x = fsolve(f1,0.5);
L2x = fsolve(f2,1.24) ;
L3x = fsolve(f3,-1.1);

fprintf('At L1, the computed value is %f\n',f1(L1x))
fprintf('At L2, the computed value is %f\n',f2(L2x))
fprintf('At L3, the computed value is %f\n',f3(L3x))

end