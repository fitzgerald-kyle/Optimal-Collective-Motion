function output = main_BVP_solver

% This function uses a shooting method to solve the boundary value problem
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
%   output.eta : error vector between desired value of x(1) and the
%                computed value of x(1)
%   output.err : norm of output.eta
%   output.p0 : initial condition p(0) found by BVP solver

% End of time interval
tf = 1; 
% Initial condition for x(0)
x0 = [0 0 0]; 
% Desired boundary condition for x(tf)
x1 = [.5 0 0]; 
% Structure containing parameters (see parameters function below)
params = parameters; 
% Guess for initial condition p(0) 
p0 = [119.803450646856,6.31288471415678e-05,-3.14268430602785]; 

figure;
% Solve the boundary value problem
output = solve_BVP(x0,p0,x1,tf,params); 

% Plot solution
plot(output.x(:,1),output.x(:,2),'b-')
hold on
% Plot desired boundary condition
plot(x1(1),x1(2),'ro')
hold off
xlabel('x_1')
ylabel('x_2')
axis([-1 1 -1 1])

end

function params = parameters

% This function contains parameters that are used when solving the BVP

params.nmax = 100; % maximium number of iterations to run BVP loop
params.tol = 1e-6; % error tolerance for solution of BVP
params.maxstep = 1; % maximum step size to make in p0 when solving BVP
params.step_ratio = 0.5; % step size decrease to use in line search
params.decrease_param = 0; % minimum error decrease in line search
params.minstep = 1e-6; % minimum step size in line search
params.ode_options = odeset('AbsTol',1e-4,'RelTol',1e-4); % parameters for solving ODEs

end