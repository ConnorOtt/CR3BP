%Solves the 5th order polynomial to find the lagrangian points
%Uses Equation 12-18

function [L1x,L2x,L3x] = findLagrangianPoints(mu)

r1cubed = @(x) (x+mu)^3;
r2cubed = @(x) (x+mu-1)^3;
f  = @(x) x - ((1-mu)*(x+mu)/r1cubed(x)) - (mu*(x-1+mu)/r2cubed(x));

L1x = fsolve(f,-mu+0.5);
L2x = fsolve(f,1);
L3x = fsolve(f,-mu-1);

fprintf('At L1, the computed value is %f\n',f(L1x))
fprintf('At L2, the computed value is %f\n',f(L2x))
fprintf('At L3, the computed value is %f\n',f(L3x))

end