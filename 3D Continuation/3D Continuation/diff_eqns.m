function dY = diff_eqns(t,Y,params)

% This function contains the differential equations we are solving
% Function Inputs:
%   t : current time
%   Y : vector containing solution at time t
% Function Outputs:
%   dY : derivative of Y at time t

% Get costates m and n
mu = Y(1:6)';

% Get stability equation solution
M = reshape(Y(7:42),6,6)';
J = reshape(Y(43:78),6,6)';

% Define control input u in terms of mu
k = 1./params.c;
u = mu(1:3).*k;
v = [1 0 0];

% Equilibrium equations for moments and forces
dmu = zeros(1,6);
dmu(1:3) = cross(mu(1:3),u)+cross(mu(4:6),v);
dmu(4:6) = cross(mu(4:6),u);

% Construct the coefficient matrices in the stability equations
F11 = [0                 mu(3)*(k(3)-k(2))  mu(2)*(k(3)-k(2));...
       mu(3)*(k(1)-k(3)) 0                  mu(1)*(k(1)-k(3));...
       mu(2)*(k(2)-k(1)) mu(1)*(k(2)-k(1))  0];
F21 = Hat(mu(4:6))*diag(k); 
F = [F11 -Hat(v); F21 -Hat(u)];

G = diag([k 0 0 0]);

H = [-Hat(u) zeros(3,3); -Hat(v) -Hat(u)];

% Stability equations
dM = F*M;
dJ = G*M+H*J;

% Collect all derivatives into one vector dX
dY = [dmu reshape(dM',1,36) reshape(dJ',1,36)]';

end

function uhat = Hat(u)
% This function computes the angular rotation matrix uhat corresponding to
% a 3-dimensional vector u.
% Function Inputs:
%   u : 1x3 or 3x1 vecto
% Function Outputs:
%   uhat : 3x3 angular rotation matrix corresponding to u

uhat = [0     -u(3) u(2);...
     u(3)  0     -u(1);...
     -u(2) u(1)  0];

end