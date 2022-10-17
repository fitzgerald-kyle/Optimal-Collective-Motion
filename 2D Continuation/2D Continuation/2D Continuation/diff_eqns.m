function dY = diff_eqns(t,Y,params)

% This function contains the differential equations we are solving
% Function Inputs:
%   t : current time
%   Y : vector containing solution at time t
% Function Outputs:
%   dY : derivative of Y at time t

% Get vectors x and p
x = Y(1:3);
p = Y(4:6);

% Get jacobian matrices M and J
M = reshape(Y(7:15),3,3)';
J = reshape(Y(16:24),3,3)';

% Define control input u in terms of x and p
u = p(3);

% Compute derivatives of x and p
dx = [cos(x(3)) sin(x(3)) u];

dp = [0 0 p(1)*sin(x(3))-p(2)*cos(x(3))];

% Compute coefficient matrices for sufficient conditions
Hxx = [0 0 0;...
       0 0 0;...
       0 0 -p(1)*cos(x(3))-p(2)*sin(x(3))];
   
Hpp = [0 0 0;...
       0 0 0;...
       0 0 1];
   
Hxp = [0 0 -sin(x(3));...
       0 0 cos(x(3));...
       0 0 0];
   
Hpx = Hxp';
 
% Compute derivatives of M and J
dM = -Hpx*M-Hxx*J;
dJ = Hpp*M+Hxp*J;

% Store derivatives
dY = [dx dp reshape(dM',1,9) reshape(dJ',1,9)]';

end