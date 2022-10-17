function output = main_IVP_solver

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

% Initial condition for x and p
q0 = eye(4);
mu0 = [0 0 1 0 0 0];

% Final time
tf = 1;

% Define parameters
params = parameters;

% Solve IVP
output = solve_IVP(q0,mu0,tf,params);

% Plot solution
plot_function(output,output.q(:,:,end))

end