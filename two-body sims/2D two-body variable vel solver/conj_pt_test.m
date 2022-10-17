function [value,isterminal,direction] = conj_pt_test(t,Y)
% Inputs t = arc length, Y = solution vector at arc length t
% Computes det(J) to see if suff cond's were met. 
% If det(J)=0 at any point, found a non-minimum local extrema == conj pt

isterminal = [0; 0]; % Don't stop solving when detJ = 0
direction = [0; 0]; % We want to know every time detJ = 0, not just pos or 
                    % neg slope instances (if det is incr or decr)

J = reshape(Y(49:84),6,6)';
detJ = det(J);

% Discount any conj points found until det(J(t))>conjtol
% Once det(J)>conjtol, we know it's not mistaking tiny det(J) for 0
conjtol = 1e-50;
value = [detJ; abs(detJ)-conjtol];

end