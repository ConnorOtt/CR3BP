function dfdt = cr3bp_eom(F, mu)
% -------------------------------------------------------------------------
% Define EOM for 3 body problem
%
% Inputs:
%   F(1:3): [3x1][km]       - Vector containing position of 3rd body in 
%                             reference to barycenter of system in 
%                             rotating p1-p2 frame.
%   F(4:6): [3x1][km/s]     - Vector containing velocity of 3rd body.
%   mu:     [scalar][~]     - Unitless mass ration M2/(M1 + M2).
%
% Outputs:
%   dfdt    [6x1][km, km/s] - Vector containing velocity and
%                             acceleration vectors of 3rd body.
%
% Authors: Brendan Boyd, Connor Ott
% -------------------------------------------------------------------------

r_vec = F(1:3);
v_vec = F(4:6);
drdt = v_vec;

% Pos :)
x = r_vec(1);
y = r_vec(2);
z = r_vec(3);

% Vel :D
x_dot = v_vec(1);
y_dot = v_vec(2);
z_dot = v_vec(3);

% Other values ;)
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt(x^2 + (y + (1-mu))^2 + z^2); % This might be incorrect

x_ddot = 2*y_dot + x - (1-mu)*(x+mu)*r1^3 - mu*(x-1+mu)/r2^3;
y_ddot = -2*x_dot + y - (1-mu)*y/r1^3 - mu*y/r2^3;
z_ddot = -(1-mu)*z/r1^3 - mu*z/r2^3; % Pretty sure this is also wrong


a_vec = [x_ddot, y_ddot, z_ddot]';

dfdt = [v_vec; a_vec];

end