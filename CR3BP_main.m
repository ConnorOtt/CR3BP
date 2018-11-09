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

% Set your initial guess (HOW???)
dx = -0.05;
dy = -0.05;
x0 = 0.5-mu + dx;
y0 = sqrt(3)/2 + dy;
u0 = -0.1;
v0 = 0.1;
IG = [x0,y0,0,u0,v0,0];
%Decent initial guess:
%   dx = -0.05; dy = -0.05; u0 = -0.1; v0 = 0.1;

k1 = 0.1;
k2 = 0.1;
TG = 1e2;
[IC,FC,T] = singleShooting(mu,IG,TG,k1,k2);


% [t_s, F_s] = ode113(@(t, F)cr3bp_eom(t, F, mu), t_span, Y0);
% 
% figure; hold on; 
% plot(0.5-mu,sqrt(3)/2,'xk');
% plot(-mu,0,'*k','linewidth',10);
% plot(1-mu,0,'*k','linewidth',7);
% plot(F_s(:,1),F_s(:,2));
% plot(F_s(1,1),F_s(1,2),'*g','linewidth',3);
% plot(F_s(end,1),F_s(end,2),'*r','linewidth',3);