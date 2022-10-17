function output_IVP = solve_IVP(x0,p0,tf,params)

% This function solves the IVP at the current value of p0
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   x1 : desired boundary condition for x(tf)
%   tf : end of time interval
%   params : structure containing parameters
% Function Outputs:
%   output_IVP.t : nx1 vector of time values t along solution
%   output_IVP.x : nx3 vector containing the solution for x(t) 
%              (e.g., output.x(:,2) is the solution for x2(t))
%   output_IVP.p : nx3 vector containing the solution for p(t)
%   output_IVP.M : 3x3xn matrix containting Jacobian of p(t) with respect to p(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.p(10,:)
%              with respect to a change in output.p(1,:))
%   output_IVP.J : Jacobian of x(t) with respect to p(0)
%   output.tconj : vector containing times of conjugate points

close all

% Initial conditions for Jacobian matrices M(t) and J(t)
M0 = eye(3);
J0 = zeros(3,3);
% Initial condition vector
Y0 = [x0 p0 reshape(M0',1,9) reshape(J0',1,9)];

% Solve ODEs
[t,sol,tconj,~,ie] = ode45(@(t,Y) diff_eqns(t,Y,params),[0 tf],Y0,params.ode_options);

% Store vectors t, x, and p
output_IVP.t = t;
output_IVP.x = sol(:,1:3);
output_IVP.p = sol(:,4:6);

% Store jacobian matrices M and J
output_IVP.M = permute(reshape(sol(:,7:15)',3,3,length(t)),[2 1 3]);
output_IVP.J = permute(reshape(sol(:,16:24)',3,3,length(t)),[2 1 3]);

% Store times of conjugate points and values of detJ
output_IVP = find_detJ(output_IVP);
tconjpts = is_stable(tconj,ie,output_IVP);
output_IVP.tconj = tconjpts;

x= output_IVP.x;
p= output_IVP.p;

[output_IVP.H_t1t1, output_IVP.H_t1t1_idx] = min( sign(p(:,1)).*sqrt(p(:,1).^2 + p(:,2).^2).* ...
    cos(x(:,3)-atan(p(:,2)./p(:,1))) );
output_IVP.p1rad = sqrt(p(1,1).^2 + p(1,2).^2);
%
figure;
plot(output_IVP.t, output_IVP.detJ);
hold on
plot(0:1, [0,0], '-k'); % solid black line at det(J)=0
hold off
drawnow

figure;
plot_function(output_IVP, output_IVP.x(end,:));
%}

end