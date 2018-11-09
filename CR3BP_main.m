% -------------------------------------------------------------------------
%
%
%                                 CR3BP       
%
% brendan connor
% -------------------------------------------------------------------------


% Propagation
t0 = 0;         % [s]
t_end = 1e2; % [s]
t_span = [t0, t_end];

M1 = 1.9891e30;
M2 = 1.8988e27;
mu = M2/(M1+M2);
mu = 0.3;

[L1x,L2x,L3x] = findLagrangianPoints3(mu);
x0 = 1-L1x-mu;
y0 = 0;
% r1 = sqrt((x0 + mu)^2 + y0^2);
% r2 = sqrt((x0-1+mu)^2 + y0^2);
% V2 = (x0^2+y0^2) + 2*( (1-mu)/r1 + mu/r2 ) - CJ;
Y0 = [x0,y0,0,0,0,0];


[t_s, F_s] = ode113(@(t, F)cr3bp_eom(t, F, mu), t_span, Y0);

figure; hold on; 
plot(-mu,0,'*k','linewidth',10);
plot(1-mu,0,'*k','linewidth',7);
plot(F_s(:,1),F_s(:,2));
plot(F_s(1,1),F_s(1,2),'*g','linewidth',3);
plot(F_s(end,1),F_s(end,2),'*r','linewidth',3);