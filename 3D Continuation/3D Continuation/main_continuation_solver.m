function output = main_continuation_solver

% Function Outputs:
%   (output(i) corresponds to the boundary condition in QF(i,:))
%
%   output(i).t : nx1 vector of time values t along solution
%   output(i).q : 4x4xn matrix containing the solution for q(t) 
%              (e.g., output.q(:,:,10) is the solution for q(t(10))
%   output(i).mu : nx6 vector containing the solution for mu(t)
%   output(i).M : 6x6xn matrix containting Jacobian of mu(t) with respect to mu(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.mu(10,:)
%              with respect to a change in output.mu(1,:))
%   output(i).J : Jacobian of q(t) with respect to mu(0)
%   output(i).eta : error vector between desired value of q(1) and the
%                computed value of q(1)
%   output(i).err : norm of output.eta
%   output(i).mu0 : initial condition mu(0) found by BVP solver
%   output(i).detJ : 1xn vector containing determinant of J(t)
%   output(i).tconj : vector containing times of conjugate points

% End of time interval
tf = 1;

% Initial condition for q(0)
q0 = eye(4);

% List of Desired boundary condition for q(tf)
n = 20;
QF = repmat(eye(4),1,1,n);
QF(1,4,:) = 0.5*ones(1,1,n); % set desired x at tf
QF(3,4,:) = linspace(.2,.5,n); % set desired z at tf

% Structure containing parameters (see parameters function below)
params = parameters; 

% Guess for initial condition mu(0) 
mu0 = [0 0 1 0 0 0];

% Solve the continuation problem
output = solve_continuation(q0,mu0,QF,tf,params);

end