%Brendan Boyd and Connor Ott
%ASEN 5050  Semester Project
%
%This function does the single shooting to try to find a periodic orbit
%
%Inputs:
%   mu - mu
%   IG - Initial guess for the state vector
%   TG - Initial guess for orbit period
%   k1 - Proportional gain for the r vector
%   k2 - Proportional gain for the v vector
%Outputs:
%   IC - The resulting initial state vector
%   FC - The resulting final state vector 
%   T  - The period of the orbit, if applicable (-1 if it didn't converge)
%--------------------------------------------------------------------------

function [IC,FC,T] = singleShooting(mu, IG, TG, k1, k2)

guessDifference = 10;
tol = 0.01;
count = 0;

while (norm(guessDifference)>tol && count<15)
    %Specify the time span to iterate over
    t_span = [0,TG];

    [t_s, F_s] = ode113(@(t, F)cr3bp_eom(t, F, mu), t_span, IG);

    r_s = F_s(:,1:3);
    v_s = F_s(:,4:6);

    %Use the velocity direction to figure out when approximately one orbit has
    %passed
    % vd = v_s-v_s(1,:);   %NOTE: This may need to change to make nv(1,:) a full matrix
    vd = v_s;
    %Find when the signs of x and y velocity change 
    xVelInds1 = find(vd(1:end-1,1)>0 & vd(2:end,1) < 0);  %Goes from - to +
    xVelInds2 = find(vd(1:end-1,1)<0 & vd(2:end,1) > 0);  %Goes from + to -
    yVelInds1 = find(vd(1:end-1,2)>0 & vd(2:end,2) < 0);
    yVelInds2 = find(vd(1:end-1,2)<0 & vd(2:end,2) > 0);
    xVelInds = [xVelInds1;xVelInds2];
    xVelInds = sort(xVelInds);
    yVelInds = [yVelInds1;yVelInds2];
    yVelInds = sort(yVelInds);
    %When the x and y differences have both changed twice, it should be
    %approximately aligned with the initial velocity vector
    ind = max([xVelInds(2),yVelInds(2)]) +1; %-1 to account for length diffs
    %Use this index to figure out which nearby F_s vector is closest to the
    %initial vector
    nNear = 3;
    finalF = F_s(ind-nNear:ind+nNear,:);
    finalFDiff = finalF - IG;
    finalFDiffNorm = vecnorm(finalFDiff');
    [~,indF] = min(finalFDiffNorm);
    indLegit = (ind-nNear) + indF;

    %Now that we have the index of the closest one after ~1 period, retreive
    %vals
    T = t_s(indLegit);
    FC = F_s(indLegit,:);

    %Now we know where the orbit ends and we just gotta do it again with an
    %updated guess

    %Update initial guess
    guessDifference = FC-IG;
    newRG = IG(1:3) + k1*guessDifference(1:3);
    newVg = IG(4:6) + k2*guessDifference(4:6);
    IG = [newRG,newVg];
    TG = 1.3*T; %Want to iterate past current T in case T increases, hence 1.3T
    
    %Increment count
    count = count+1;
    
% Call these to plot the trajectory with L4
figure; hold on; 
plot(0.5-mu,sqrt(3)/2,'xk');
plot(-mu,0,'*k','linewidth',10);
plot(1-mu,0,'*k','linewidth',7);
plot(r_s(:,1),r_s(:,2));
plot(FC(1),FC(2),'*r');
plot(r_s(ind,1),r_s(ind,2),'*y');
plot(r_s(1:indLegit,1),r_s(1:indLegit,2));
end

IC = IG;

% 
% % Call these to plot the trajectory with L4
% figure; hold on; 
% plot(0.5-mu,sqrt(3)/2,'xk');
% plot(-mu,0,'*k','linewidth',10);
% plot(1-mu,0,'*k','linewidth',7);
% plot(r_s(:,1),r_s(:,2));
% plot(FC(1),FC(2),'*r');
% plot(r_s(ind,1),r_s(ind,2),'*y');
% plot(r_s(1:indLegit,1),r_s(1:indLegit,2));

end