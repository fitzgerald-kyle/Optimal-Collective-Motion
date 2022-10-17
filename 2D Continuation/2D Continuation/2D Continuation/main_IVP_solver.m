function output = main_IVP_solver

% This function uses solves the initial value problem 
% for the following system:
% x1' = sin(x3)   x2' = cos(x3)   x3' = p3
% p1' = 0   p2' = 0   p3' = p1*sin(x3) - p2*cos(x3)
% Function Outputs:
%   output.t : nx1 vector of time values t along solution
%   output.x : nx3 vector containing the solution for x(t) 
%              (e.g., output.x(:,2) is the solution for x2(t))
%   output.p : nx3 vector containing the solution for p(t)
%   output.M : 3x3xn matrix containing Jacobian of p(t) with respect to p(0)
%              (e.g., output.M(:,:,10) is the Jacobian of output.p(10,:)
%              with respect to a change in output.p(1,:))
%   output.J : Jacobian of x(t) with respect to p(0)
%   output.tconj : vector containing times of conjugate points

% Initial condition for x and p
x0 = [0 0 0];
p0 = [0 0 1];

% Final time
tf = 1;

% Define parameters
params = parameters;

% Solve IVP
output = solve_IVP(x0,p0,tf,params);

% Plot solution
plot_function(output,[output.x(end,1) output.x(end,2)])

end