function output = main_IVP_solver(mu0,tf,s,showplot)
% Function Inputs:
%   mu0 : 1x12 vector of initial values mu(0)
%   s : scalar coupling strength
%   showplot : boolean that determines whether plot_function is called
% Function Outputs:
%   output_IVP.t : nx1 vector of time values t along solution
%   output_IVP.q : 8x8xn matrix containing the solution for q(t) 
%              (e.g., output.q(:,:,10) is the solution for q(t(10))
%   output_IVP.mu : nx12 vector containing the solution for mu(t)
%   output_IVP.M : 12x12xn matrix containting Jacobian of mu(t) with respect to mu(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.mu(10,:)
%              with respect to a change in output.mu(1,:))
%   output_IVP.J : Jacobian of q(t) with respect to mu(0)
%   output.tconj : vector containing times of conjugate points

%close all

% Initial condition for x and p
q0 = eye(8);

% Define parameters
params = parameters;

% Solve IVP
output = solve_IVP(q0,mu0,tf,params,s);

% Plot solution
if showplot
    figure;
    plot_function(output,output.q(:,:,end));
end

end