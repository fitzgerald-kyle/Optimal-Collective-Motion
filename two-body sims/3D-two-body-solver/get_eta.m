function eta = get_eta(q,qf)

% Get current positions and orientations at t=tf
r1 = q(1:3,4);
R1 = q(1:3,1:3);
r2 = q(5:7,8);
R2 = q(5:7,5:7);

% Get desired positions and orientations at t=tf
r1f = qf(1:3,4);
R1f = qf(1:3,1:3);
r2f = qf(5:7,8);
R2f = qf(5:7,5:7);

% Compute errors between current and desired positions
dr1 = R1'*(r1-r1f);
dr2 = R2'*(r2-r2f);

% Compute errors between current and desired orientations
dR1 = get_angle_err(R1f'*R1);
dR2 = get_angle_err(R2f'*R2);

% Collect position and orientation errors
eta=[dR1; dr1; dR2; dr2];

end

function dR = get_angle_err(Rerr)

% This function computes the error between the current and desired 
% orientations at the end of the strip

dR2 = asin(Rerr(1,3));
dR1 = atan2(-Rerr(2,3),Rerr(3,3));
dR3 = atan2(-Rerr(1,2),Rerr(1,1));
dR = [dR1; dR2; dR3];

end