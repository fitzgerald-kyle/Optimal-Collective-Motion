function q = get_q(q0,mu,t,params,s)

% This function finds the state trajectory q given the costate trajectory
% mu
% Function Inputs:
%   q0 : initial condition for q(0)
%   mu : costate trajectory
%   t : vector of time values along solution
%   params : structure containing parameters
%   s : scalar coupling strength
% Function Outputs:
%   q : state trajectory

% Compute velocity matrix
X = @(u1,u2,v) blkdiag([0     -u1(3) u1(2)  v(1);...
            u1(3)  0     -u1(1) v(2);...
            -u1(2) u1(1)  0     v(3);...
            0     0     0     0],...
            [0     -u2(3) u2(2)  v(1);...
            u2(3)  0     -u2(1) v(2);...
            -u2(2) u2(1)  0     v(3);...
            0     0     0     0]);
       
% Initialize centerline q   
q = zeros(8,8,length(t));
q(:,:,1) = q0;

% Compute step sizes
dt = t(2:end)-t(1:end-1);

% Compute centerline shape q
for i = 1:length(t)-1
    
    u1 = 1/(2*s+1) * ((1+s)*mu(i,1:3) + s*mu(i,7:9));
    u2 = 1/(2*s+1) * (s*mu(i,1:3) + (1+s)*mu(i,7:9));
    v = [1 0 0];
    
    % Compute next pose
    q(:,:,i+1) = q(:,:,i)*expm(dt(i)*X(u1,u2,v));
    
end

end