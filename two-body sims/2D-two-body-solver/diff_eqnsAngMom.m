function dY = diff_eqnsAngMom(t,Y,mu)

% This function contains the differential equations we are solving
% Function Inputs:
%   t : current time
%   Y : vector containing solution at time t
% Function Outputs:
%   dY : derivative of Y at time t

% Get vectors x and p
x = Y(1:6);
p = Y(7:12);
% Get jacobian matrices M and J
M = reshape(Y(13:48),6,6)';
J = reshape(Y(49:84),6,6)';

% Compute derivatives of x and p
dx1= [cos(x(3)) sin(x(3)) ((1+mu)*p(3) + mu*p(6))/(2*mu+1)];
dx2= [cos(x(6)) sin(x(6)) (mu*p(3) + (1+mu)*p(6))/(2*mu+1)];
dx= [dx1 dx2];

dp1= [0 0 p(1)*sin(x(3))-p(2)*cos(x(3))];
dp2= [0 0 p(4)*sin(x(6))-p(5)*cos(x(6))];
dp = [dp1 dp2];

% Compute coefficient matrices for jacobian equations
Hxx = zeros(6,6);
Hxx(3,3)= -p(1)*cos(x(3))-p(2)*sin(x(3));
Hxx(6,6)= -p(4)*cos(x(6))-p(5)*sin(x(6));

Hpp = zeros(6,6);
Hpp(3,3)= (mu + 1)/(2*mu + 1);
Hpp(6,6)= (mu + 1)/(2*mu + 1);
Hpp(3,6)= mu/(2*mu + 1);
Hpp(6,3)= Hpp(3,6);

Hxp= zeros(6,6);
Hxp(1,3)= -sin(x(3));
Hxp(2,3)= cos(x(3));
Hxp(4,6)= -sin(x(6));
Hxp(5,6)= cos(x(6));
   
Hpx = Hxp';
 
% Compute derivatives of M and J
dM = -Hpx*M-Hxx*J;
dJ = Hpp*M+Hxp*J;

% Store derivatives
dY = [dx dp reshape(dM',1,36) reshape(dJ',1,36)]';

end