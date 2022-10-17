function uhat = Hat(u)
% This function computes the angular rotation matrix corresponding to
% a 3-dimensional vector u.
% Function Inputs:
%   u : 1x3 or 3x1 vector
% Function Outputs:
%   uhat : 3x3 angular rotation matrix corresponding to u

uhat = [0     -u(3) u(2);...
     u(3)  0     -u(1);...
     -u(2) u(1)  0];

end