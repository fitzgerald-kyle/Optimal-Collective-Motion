function output_BVP = solve_BVP(x0,p0,x1,tf,params)

% This function contains the main loop used to solve the BVP using a
% shooting method
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   x1 : desired boundary condition for x(tf)
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

%% Initialize counter and output
i = 0;
output_BVP = [];

%% Boundary value problem loop
while i <= params.nmax 
    
    % Iterate counter
    i = i+1;
     
    % Solve IVP
    output_IVP = solve_IVP(x0,p0,x1,tf,params);
    
    % Print variables of interest
    fprintf('iteration: %.0f     error: %.10f \n',i,output_IVP.err)
    
    % Plot x1 vs x2 trajectory as the BVP is solved
    % (comment this section for faster performance)
    plot(output_IVP.x(:,1),output_IVP.x(:,2),'b-')
    hold on
    % Plot desired boundary condition
    plot(x1(1),x1(2),'ro')
    hold off
    axis([-1 1 -1 1])
    xlabel('x_1')
    ylabel('x_2')
    drawnow
    
    % Check if error tolerance is satisfied
    if output_IVP.err <= params.tol
        output_BVP = output_IVP;
        output_BVP.p0 = p0;
        break
    end
    
    % Compute search direction for p0
    jacobian = output_IVP.J(:,:,end);
    dp0 = jacobian\output_IVP.eta';
    
    % Enforce maximum step size
    if norm(dp0) > params.maxstep
        dp0 = dp0/norm(dp0)*params.maxstep;
    end
    
    % Use line search to find step size that ensures error decreases
    step = line_search(x0,p0,x1,dp0,tf,output_IVP.err,params);
    
    % Update p0
    if ~isempty(step)
        p0 = p0+step*dp0';
    else
        error('line search failed')
    end

end

if output_IVP.err > params.tol
    error('BVP solver failed')
end

end

function step = line_search(x0,p0,x1,dp0,tf,err,params)

% This function uses a line search to ensure that the error decreases when
% updating p0
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   x1 : desired boundary condition for x(tf)
%   dp0 : search direction for p0
%   tf : end of time interval
%   err : err at current value of p(0)
%   params : structure containing parameters
% Function Outputs:
%   step :  step size that ensures error decreases in the direction dp0

% Initialize step size and new error
step = 1/params.step_ratio;
newerr = 2*err;
    
% Line search loop
while newerr > (1-params.decrease_param)*err
    
    % Compute step size and p0 after step
    step = params.step_ratio*step;
    newp0 = p0+step*dp0';
    
    % Compute error after step
    output_IVP = solve_IVP(x0,newp0,x1,tf,params);
    newerr = output_IVP.err;
    
    % Check if step is less than minstep
    if norm(step) < params.minstep
        step = [];
        break
    end 
end
    
end