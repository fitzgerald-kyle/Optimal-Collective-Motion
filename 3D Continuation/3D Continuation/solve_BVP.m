function output_BVP = solve_BVP(q0,mu0,qf,tf,params)

% This function contains the main loop used to solve the BVP using a
% shooting method
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu0 : guess for initial condition mu(0) 
%   qf : desired boundary condition for q(tf)
%   tf : end of time interval
%   params : structure containing parameters
% Function Outputs:
%   output_BVP.t : nx1 vector of time values t along solution
%   output_BVP.q : 4x4xn matrix containing the solution for q(t) 
%              (e.g., output.q(:,:,10) is the solution for q(t(10))
%   output_BVP.mu : nx6 vector containing the solution for mu(t)
%   output_BVP.M : 6x6xn matrix containting Jacobian of mu(t) with respect to mu(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.mu(10,:)
%              with respect to a change in output.mu(1,:))
%   output_BVP.J : Jacobian of q(t) with respect to mu(0)
%   output_BVP.eta : error vector between desired value of q(1) and the
%                computed value of q(1)
%   output_BVP.err : norm of output.eta
%   output_BVP.mu0 : initial condition mu(0) found by BVP solver
%   output.tconj : vector containing times of conjugate points

%% Initialize counter and output
i = 0;
output_BVP = [];

%% Boundary value problem loop
while i <= params.nmax 
    
    % Iterate counter
    i = i+1;
     
    % Solve IVP
    output_IVP = solve_IVP(q0,mu0,tf,params);
    
%     % Store error between current and desired value of q(tf)
    output_IVP.eta = get_eta(output_IVP.q(:,:,end),qf);
    output_IVP.err = norm(output_IVP.eta);
    
    % Print variables of interest
    fprintf('iteration: %.0f     error: %.10f \n',i,output_IVP.err)
    
    % Plot solution
    plot_function(output_IVP,qf)
    
    % Check if error tolerance is satisfied
    if output_IVP.err <= params.tol
        output_BVP = output_IVP;
        output_BVP.mu0 = mu0;
        break
    end
    
    % Compute search direction for mu0
    jacobian = output_IVP.J(:,:,end);
    dmu0 = -jacobian\output_IVP.eta;
    
    % Enforce maximum step size
    if norm(dmu0) > params.maxstep
        dmu0 = dmu0/norm(dmu0)*params.maxstep;
    end
    
    % Update mu0
%     mu0 = mu0+dmu0';
    
    % Use line search to find step size that ensures error decreases
    [step,dmu0] = line_search(q0,mu0,qf,dmu0,tf,output_IVP.err,params);
    
    % Update mu0
    if ~isempty(step)
        mu0 = mu0+step*dmu0';
    else
        error('line search failed')
    end

end

if output_IVP.err > params.tol
    error('BVP solver failed')
end

end

function [step,dmu0] = line_search(q0,mu0,qf,dmu0,tf,err,params)

% This function uses a line search to ensure that the error decreases when
% updating mu0
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu0 : guess for initial condition mu(0) 
%   qf : desired boundary condition for q(tf)
%   dmu0 : search direction for mu0
%   tf : end of time interval
%   err : err at current value of mu(0)
%   params : structure containing parameters
% Function Outputs:
%   step :  step size that ensures error decreases in the direction dmu0

% Initialize step size and new error
step = 1/params.step_ratio;
newerr = 2*err;
    
% Line search loop
while newerr > (1-params.decrease_param)*err
    
    % Compute step size and p0 after step
    step = params.step_ratio*step;
    newmu0 = mu0+step*dmu0';
    
    % Compute error after step
    output_IVP = solve_IVP(q0,newmu0,tf,params);
    
    % Store error between current and desired value of q(tf)
    neweta = get_eta(output_IVP.q(:,:,end),qf);
    newerr = norm(neweta);
    
    % Check if step is less than minstep
    if norm(step) < params.minstep
        step = 1;
        dmu0 = rand(size(dmu0))-1;
        break
    end 
end
    
end