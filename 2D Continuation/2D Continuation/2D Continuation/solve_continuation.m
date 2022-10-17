function output = solve_continuation(x0,p0,XF,tf,params)

% This function contains the main loop used to solve the continuation
% problem
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   XF : list of desired boundary condition for x(tf)
%   tf : end of time interval
%   params : structure containing parameters
% Function Outputs:
%   output_BVP.t : nx1 vector of time values t along solution
%   output_BVP.x : nx3 vector containing the solution for x(t) 
%              (e.g., output.x(:,2) is the solution for x2(t))
%   output_BVP.p : nx3 vector containing the solution for p(t)
%   output_BVP.M : 3x3xn matrix containing Jacobian of p(t) with respect to p(0)
%              (e.g., output.M(:,:,10) is the Jacobian of output.p(10,:)
%              with respect to a change in output.p(1,:))
%   output_BVP.J : Jacobian of x(t) with respect to p(0)
%   output_BVP.eta : error vector between desired value of x(1) and the
%                computed value of x(1)
%   output_BVP.err : norm of output.eta
%   output_BVP.p0 : initial condition p(0) found by BVP solver
%   output.tconj : vector containing times of conjugate points

for i=1:length(XF(:,1))
   
    % Set desired boundary conditions
    xf = XF(i,:);
    
    % Solve the boundary value problem
    output_BVP = solve_BVP(x0,p0,xf,tf,params);
    
    % Compute detJ
    output_BVP = find_detJ(output_BVP);
    
    % Store solution of boundary value problem
    output(i) = output_BVP;

    % Update initial guess for p0
    p0 = output_BVP.p0;
    
end

end