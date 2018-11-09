%Solves the 5th order polynomial to find the lagrangian points

function [L1x,L2x,L3x] = findLagrangianPoints(mu)

r1cubed = @(x) (x+mu)^3;
r2cubed = @(x) (x+mu-1)^3;
f  = @(x) x - (1-mu)*(x+mu)/r1cubed(x) - mu*(x-1+mu)/r2cubed(x);

L1x = fsolve(f,-mu+0.5);
L2x = fsolve(f,1);
L3x = fsolve(f,-mu-1);

end