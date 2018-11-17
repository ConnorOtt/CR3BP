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
%   dx = -0.05; dy = -0.05; u0 = -0.1; v0 = 0.1; T = 6;

% k1 = 0.1;
% k2 = 0.1;
TG = 6;
% [IC,FC,T] = singleShooting(mu,IG,TG,k1,k2);

%Run grid shooting
tol = 0.0001;
minDist = 100;
du = abs(u0)/10;
dv = abs(v0)/10;
count=0;
while (minDist>tol && du>0 && dv>0)
    %Do a grid shoot
    [IC, minDist, T] = gridShoot(mu,IG,TG,du,dv);
    
    %Update values to do the grid shoot again
    if (abs(IG(4)-IC(4)) == 0)
        du = du/10;
    else
        du = abs(IG(4)-IC(4))/10;
    end
    if (abs(IG(5)-IC(5)) == 0)
        dv = dv/10;
    else
        dv = abs(IG(5)-IC(5))/10;
    end
    IG = IC;
    TG = T;
    count = count+1;
end

%Plot
[t_s, F_s] = ode113(@(t, F)cr3bp_eom(t, F, mu), [0,3*T], IC);

figure; hold on; 
plot(0.5-mu,sqrt(3)/2,'xk');
plot(-mu,0,'*k','linewidth',10);
plot(1-mu,0,'*k','linewidth',7);
plot(F_s(:,1),F_s(:,2));
plot(F_s(1,1),F_s(1,2),'*g','linewidth',3);
plot(F_s(end,1),F_s(end,2),'*r','linewidth',3);