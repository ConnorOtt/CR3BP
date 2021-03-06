%Brendan Boyd and Connor Ott
%ASEN 5050  Semester Project
%
%This function does a single grid shoot to try to find a periodic orbit 
%Note that only the initial velocities are varied.  The initial position is
%kept constant.
%
%Inputs:
%   mu - mu
%   IG - Initial guess for the state vector
%   TG - Initial guess for orbit period
%Outputs:
%   IC - The resulting initial state vector
%   minDist - The distance scalar at the resulting IC
%   T  - The period of the orbit
%--------------------------------------------------------------------------

function [IC, minDist, T] = gridShoot(mu,IG,TG,du,dv)

%Set the options for the ODE
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

%Pick out the parameters to be varied
u = IG(4); v = IG(5);
% du = u/100;
% dv = v/100;

%Set up the grid
uv = (u-10*du):du:(u+10*du);
vv = (v-10*dv):dv:(v+10*dv);
N = length(uv);
distances = zeros(N,N);
Ts = zeros(N,N);

%Iterate
for i=1:N
    for j=1:N
        %Pick out the associated initial state
        IS = [IG(1:3),uv(i),vv(j),IG(6)];
        
        %Run the ODE
        [t_s,F_s] = ode113(@(t, F)cr3bp_eom(t, F, mu), [0,TG*(1.2)], IS, options);
        
        %Find the closest location and calculate the distance
        dF = zeros(length(t_s),1);
        for k=1:length(t_s)
            dF(k) = norm(F_s(k,:)-IS);
        end
        [distances(i,j),ind] = min(dF(floor(length(dF)*3/4):end));
        ind = ind+floor(length(dF)*3/4)-1;
        Ts(i,j) = t_s(ind);
        
%         % Call these to plot the trajectory with L4
%         figure; hold on; 
%         plot(0.5-mu,sqrt(3)/2,'xk');
%         plot(-mu,0,'*k','linewidth',10);
%         plot(1-mu,0,'*k','linewidth',7);
%         plot(F_s(:,1),F_s(:,2));
    end
end

%Find the smallest value distance and return u and v
[minDist,rowInd] = min(distances);
[minDist,colInd] = min(minDist);
rowInd = rowInd(colInd);
T = Ts(rowInd,colInd);
IC = [IG(1:3),uv(rowInd),vv(colInd),IG(6)];

end