function output = solve_continuation(q0,mu0,QF,tf,params)

% This function contains the main loop used to solve the continuation
% problem
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu0 : guess for initial condition mu(0) 
%   QF : list of desired boundary condition for q(tf)
%   tf : end of time interval
%   params : structure containing parameters
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

for i=1:length(QF(1,1,:))
   
    % Set desired boundary conditions
    qf = QF(:,:,i);
    
    % Solve the boundary value problem
    output_BVP = solve_BVP(q0,mu0,qf,tf,params); 
    
    % Compute detJ
    output_BVP = find_detJ(output_BVP);
    
    % Store solution of boundary value problem
    output(i) = output_BVP;

    % Update initial guess for p0
    mu0 = output_BVP.mu0;
    
end

end