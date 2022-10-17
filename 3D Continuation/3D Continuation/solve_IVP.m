function output_IVP = solve_IVP(q0,mu0,tf,params)

% This function solves the IVP at the current value of p0
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu0 : guess for initial condition mu(0) 
%   tf : end of time interval
%   params : structure containing parameters
% Function Outputs:
%   output_IVP.t : nx1 vector of time values t along solution
%   output_IVP.q : 4x4xn matrix containing the solution for q(t) 
%              (e.g., output.q(:,:,10) is the solution for q(t(10))
%   output_IVP.mu : nx6 vector containing the solution for mu(t)
%   output_IVP.M : 6x6xn matrix containting Jacobian of mu(t) with respect to mu(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.mu(10,:)
%              with respect to a change in output.mu(1,:))
%   output_IVP.J : Jacobian of q(t) with respect to mu(0)
%   output.tconj : vector containing times of conjugate points

% Initial conditions for Jacobian matrices M(t) and J(t)
M0 = eye(6);
J0 = zeros(6,6);

% Initial condition vector
Y0 = [mu0 reshape(M0',1,36) reshape(J0',1,36)];

% Solve ODEs
[t,sol,tconj] = ode45(@(t,Y) diff_eqns(t,Y,params),[0 tf],Y0,params.ode_options);

% Store vectors t, x, and p
output_IVP.t = t;
output_IVP.mu = sol(:,1:6);
output_IVP.q = get_q(q0,output_IVP.mu,t,params);

% Store jacobian matrices M and J
output_IVP.M = permute(reshape(sol(:,7:42)',6,6,length(t)),[2 1 3]);
output_IVP.J = permute(reshape(sol(:,43:78)',6,6,length(t)),[2 1 3]);

% Store times of conjugate points
output_IVP.tconj = tconj;

end