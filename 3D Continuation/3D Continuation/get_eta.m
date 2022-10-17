function eta = get_eta(q1,qf)

% Get current position and orientation at t=tf
r1 = q1(1:3,4);
R1 = q1(1:3,1:3);

% Get desired position and orientation at t=tf
rf = qf(1:3,4);
Rf = qf(1:3,1:3);

% Compute error between current and desired position
dr = R1'*(r1-rf);

% Compute error between current and desired orientation
dR = get_angle_err(Rf'*R1);

% Collect position and orientation errors
eta=[dR; dr];

end

function dR = get_angle_err(Rerr)

% This function computes the error between the current and desired 
% orientation at the end of the strip

dR2 = asin(Rerr(1,3));
dR1 = atan2(-Rerr(2,3),Rerr(3,3));
dR3 = atan2(-Rerr(1,2),Rerr(1,1));
dR = [dR1; dR2; dR3];

end