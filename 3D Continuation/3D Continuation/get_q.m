function q = get_q(q0,mu,t,params)

% This function finds the state trajectory q given the costate trajectory
% mu
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu : costate trajectory
%   t : vector of time values along solution
%   params : structure containing parameters
% Function Outputs:
%   q : state trajectory

% Compute velocity matrix
X = @(u,v) [0     -u(3) u(2)  v(1);...
            u(3)  0     -u(1) v(2);...
            -u(2) u(1)  0     v(3);...
            0     0     0     0];
       
% Initialize centerline q   
q = zeros(4,4,length(t));
q(:,:,1) = q0;

% Compute step sizes
dt = t(2:end)-t(1:end-1);

% Compute centerline shape q
for i = 1:length(t)-1
    
    u = mu(i,1:3)./params.c;
    v = [1 0 0];
    
    % Compute next pose
    q(:,:,i+1) = q(:,:,i)*expm(dt(i)*X(u,v));
    
end

end