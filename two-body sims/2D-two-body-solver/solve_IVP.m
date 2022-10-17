function output_IVP = solve_IVP(x0,p0,tf,params,mu,L)

% This function solves the IVP at the current value of p0
% Function Inputs:
%   x0 : initial condition for x(0)
%   p0 : guess for initial condition p(0) 
%   tf : end of time interval
%   params : structure containing parameters
%   mu : scalar coupling strength
%   L : system scaling coefficient
% Function Outputs:
%   output_IVP.t : nx1 vector of time values t along solution
%   output_IVP.x : nx6 vector containing the solution for x(t) 
%              (e.g., output.x(:,2) is the solution for x2(t))
%   output_IVP.p : nx6 vector containing the solution for p(t)
%   output_IVP.M : 6x6xn matrix containting Jacobian of p(t) with respect to p(0)
%              (e.g., output.M(:,:,10) is the jacobian of output.p(10,:)
%              with respect to a change in output.p(1,:))
%   output_IVP.J : Jacobian of x(t) with respect to p(0)
%   output.tconj : vector containing times of conjugate points

%close all

% Initial conditions for Jacobian matrices M(t) and J(t)
M0 = eye(6);
J0 = zeros(6,6);

% scaling (if applicable)
p0 = L^2*p0;
p0(3) = p0(3)/L;
p0(6) = p0(6)/L;

% Initial condition vector
Y0 = [x0 p0 reshape(M0',1,36) reshape(J0',1,36)];

% Solve ODEs
[t,sol,tconj,~,~] = ode45(@(t,Y) diff_eqnsAngMom(t,Y,mu),[0 tf],Y0,params.ode_options);

% Store vectors t, x, and p
output_IVP.t = t;
output_IVP.x = sol(:,1:6);
output_IVP.p = sol(:,7:12);

% Store jacobian matrices M and J
output_IVP.M = permute(reshape(sol(:,13:48)',6,6,length(t)),[2 1 3]);
output_IVP.J = permute(reshape(sol(:,49:84)',6,6,length(t)),[2 1 3]);

% Store times of conjugate points and values of detJ
output_IVP = find_detJ(output_IVP);
tconjpts = is_stable(tconj,output_IVP);
output_IVP.tconj = tconjpts;
%output_IVP.Ht1t1 = output_IVP.p(:,1).*cos(output_IVP.x(:,3)) + ...
%    output_IVP.p(:,2).*cos(output_IVP.x(:,6));

%{
figure;
plot(output_IVP.t, output_IVP.detJ);
hold on
plot(0:1, [0,0], '-k'); % solid black line at det(J)=0
hold off
drawnow
%}
%{
figure;
plot_function(output_IVP, output_IVP.x(end,:), mu);
%}

end
