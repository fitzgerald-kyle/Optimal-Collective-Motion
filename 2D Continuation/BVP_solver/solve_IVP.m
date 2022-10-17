function output_IVP = solve_IVP(x0,p0,x1,tf,params)

% This function solves the IVP at the current value of p0
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   x1 : desired boundary condition for x(tf)
%   tf : end of time interval
%   params : structure containing parameters
% Function Outputs:
%   output_IVP.t : nx1 vector of time values t along solution
%   output_IVP.x : nx3 vector containing the solution for x(t) 
%              (e.g., output.x(:,2) is the solution for x2(t))
%   output_IVP.p : nx3 vector containing the solution for p(t)
%   output_IVP.M : 3x3xn matrix containting Jacobian of p(t) with respect to p(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.p(10,:)
%              with respect to a change in output.p(1,:))
%   output_IVP.J : Jacobian of x(t) with respect to p(0)
%   output_IVP.eta : error vector between desired value of x(1) and the
%                computed value of x(1)
%   output_IVP.err : norm of output.eta

% Initial conditions for Jacobian matrices M(t) and J(t)
M0 = eye(3);
J0 = zeros(3,3);
% Initial condition vector
Y0 = [x0 p0 reshape(M0',1,9) reshape(J0',1,9)];

% Solve ODEs
[t,sol] = ode45(@diff_eqns,[0 tf],Y0,params.ode_options);

% Store vectors t, x, and p
output_IVP.t = t;
output_IVP.x = sol(:,1:3);
output_IVP.p = sol(:,4:6);
% Store jacobian matrices M and J
output_IVP.M = permute(reshape(sol(:,7:15)',3,3,length(t)),[2 1 3]);
output_IVP.J = permute(reshape(sol(:,16:24)',3,3,length(t)),[2 1 3]);

% Store error between current and desired value of x(tf)
output_IVP.eta = x1-output_IVP.x(end,:);
output_IVP.err = norm(output_IVP.eta);

end