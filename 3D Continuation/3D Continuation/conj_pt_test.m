function [detJ,isterminal,direction] = conj_pt_test(t,Y)
% This function computes the determinant of the matrix J(t).  The ODE
% solver within the function solve_IVP finds times when det(J(t))=0.
% INPUTS: t - current arc length.
%         Y - solution of equilibrium and stability equations at 
%             arc length t.
%
% OUTPUTS: detJ - value of det(J(t)).
%          isterminal - equal to 0, indicating that the ODE solver should 
%                       not terminate when detJ=0.
%          direction - equal to [], indicating that we don't care if
%                      det(J(t)) is increasing or decreasing when it 
%                      crosses 0.

isterminal = 0;
direction = [];

% If t=0, then det(J(0))=0, but this isn't a conjugate point, so we only
% check det(J(0)) when t>0
if t == 0
    detJ = 1;
else
    % Get the matrix J
    J = reshape(Y(43:78),6,6)';
    % Compute det(J(t))
    detJ = det(J);
end

end