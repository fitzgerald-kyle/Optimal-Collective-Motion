function output = main_BVP_solver

% Function Outputs:
%   output.t : nx1 vector of time values t along solution
%   output.q : 4x4xn matrix containing the solution for q(t) 
%              (e.g., output.q(:,:,10) is the solution for q(t(10))
%   output.mu : nx6 vector containing the solution for mu(t)
%   output.M : 6x6xn matrix containting Jacobian of mu(t) with respect to mu(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.mu(10,:)
%              with respect to a change in output.mu(1,:))
%   output.J : Jacobian of q(t) with respect to mu(0)
%   output.eta : error vector between desired value of q(1) and the
%                computed value of q(1)
%   output.err : norm of output.eta
%   output.mu0 : initial condition mu(0) found by BVP solver
%   output.detJ : 1xn vector containing determinant of J(t)
%   output.tconj : vector containing times of conjugate points


% End of time interval
tf = 1;

% Initial condition for q(0)
q0 = eye(4);

% Desired boundary condition for q(tf)
qf = eye(4);
qf(1,4) = 0.5; % set desired x at tf
qf(3,4) = 0.2; % set desired x at tf

% Structure containing parameters (see parameters function below)
params = parameters; 

% Guess for initial condition mu(0) 
mu0 = [0 0 1 0 0 0];

% Solve the boundary value problem
output = solve_BVP(q0,mu0,qf,tf,params); 

% Compute detJ
output = find_detJ(output);

% Plot solution
plot_function(output,qf)

end

