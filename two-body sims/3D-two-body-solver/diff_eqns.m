function dY = diff_eqns(t,Y,s)

% This function contains the differential equations we are solving
% Function Inputs:
%   t : current time
%   Y : vector containing solution at time t
%   s : scalar coupling strength
% Function Outputs:
%   dY : derivative of Y at time t

% Get costates m and n
mu = Y(1:12)';

% Get stability equation solution
M = reshape(Y(13:156),12,12)';
J = reshape(Y(157:300),12,12)';

% Define control input u in terms of mu
u1 = 1/(2*s+1) * ((1+s)*mu(1:3) + s*mu(7:9));
u2 = 1/(2*s+1) * (s*mu(1:3) + (1+s)*mu(7:9));
v = [1 0 0];

% Equilibrium equations for moments and forces
dmu = zeros(1,12);
dmu(1:3) = cross(mu(1:3),u1) + cross(mu(4:6),v);
dmu(4:6) = cross(mu(4:6),u1);
dmu(7:9) = cross(mu(7:9),u2) + cross(mu(10:12),v);
dmu(10:12) = cross(mu(10:12),u2);

% Construct the coefficient matrices in the stability equations
F = 1/(2*s+1) * ...
    [-s*Hat(mu(7:9)) -(2*s+1)*Hat(v) s*Hat(mu(1:3)) zeros(3,3);...
     (1+s)*Hat(mu(4:6)) -(2*s+1)*Hat(u1) s*Hat(mu(4:6)) zeros(3,3);...
     s*Hat(mu(7:9)) zeros(3,3) -s*Hat(mu(1:3)) -(2*s+1)*Hat(v);...
     s*Hat(mu(10:12)) zeros(3,3) (1+s)*Hat(mu(10:12)) -(2*s+1)*Hat(u2)];

G = 1/(2*s+1) * ...
    [diag([1+s 1+s 1+s 0 0 0]) diag([s s s 0 0 0]);...
     diag([s s s 0 0 0]) diag([1+s 1+s 1+s 0 0 0])];

H = blkdiag([-Hat(u1) zeros(3,3); -Hat(v) -Hat(u1)], ...
    [-Hat(u2) zeros(3,3); -Hat(v) -Hat(u2)]);

% Stability equations
dM = F*M;
dJ = G*M+H*J;

% Collect all derivatives into one vector dX
dY = [dmu reshape(dM',1,144) reshape(dJ',1,144)]';

end