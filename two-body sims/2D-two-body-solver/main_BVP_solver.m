function [output, nmaxErr] = main_BVP_solver(x0, x1, p0, mu, nmax)

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
%   output.detJ : 1xn vector containing determinant of J(t)
%   output.tconj : vector containing times of conjugate points

% End of time interval
tf = 1; 
% Initial conditions for x1(0), x2(0)
%x0_1= [0 0 0]; x0_2= [0 0 0];
%x0 = [x0_1 x0_2]; 
% Desired boundary conditions for x1(tf), x2(tf)
%x1_1= [.2 0 0]; x1_2= [.4 0 0];
%x1 = [x1_1 x1_2]; 
% Structure containing parameters (see parameters function below)
params = parameters(nmax); 
% Guess for initial conditions p1(0), p2(0)
%p0_1= [1 0 1]; p0_2= [2 0 -1];
%p0 = [p0_1 p0_2];
%p0 = [-131.2566 -0.0015 -26.5430 -120.9796 5.0304e-04 23.5063];

% Solve the boundary value problem
[output, nmaxErr] = solve_BVP(x0,p0,x1,tf,params,mu); 
if nmaxErr
    return
end
%
% Plot solution
%close all
plot_function(output,x1,mu);
%}
output = find_detJ(output);

%figure;
%plot(output.t, output.detJ)

output.inflecPts= inflectionPts(output.x(:,3), output.x(:,6));
output.mu= mu;

output.u1 = 1/(2*mu+1)*((1+mu)*output.p0(3) + mu*output.p0(6));
output.u2 = 1/(2*mu+1)*((1+mu)*output.p0(6) + mu*output.p0(3));
end

function numPts= inflectionPts(theta1, theta2)
    dtheta1= zeros(length(theta1)-1,1);
    dtheta2= zeros(length(theta2)-1,1);
    for i= 2:length(dtheta1)
        dtheta1(i)= theta1(i) - theta1(i-1);
        dtheta2(i)= theta2(i) - theta2(i-1);
    end    

    numPts1= 0; numPts2= 0;
    for i= 3:length(dtheta1)
        if sign(dtheta1(i)) == -sign(dtheta1(i-1)) || ...
                (sign(dtheta1(i)) == -sign(dtheta1(i-2)) && sign(dtheta1(i-1)) == 0)
            numPts1= numPts1 + 1;
        end
        if sign(dtheta2(i)) == -sign(dtheta2(i-1)) || ...
                (sign(dtheta2(i)) == -sign(dtheta2(i-2)) && sign(dtheta2(i-1)) == 0)
            numPts2= numPts2 + 1;
        end
    end    

    numPts= [numPts1 numPts2];
end