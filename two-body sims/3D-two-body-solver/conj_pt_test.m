function [detJ,isterminal,direction] = conj_pt_test(t,Y)
% Inputs t = arc length, Y = solution vector at arc length t
% Computes det(J) to see if suff cond's were met. 
% If det(J)=0 at any point, found a non-minimum local extrema == conj pt

isterminal = 0; % Don't stop solving when detJ = 0
direction = 0; % We want to know every time detJ = 0, not just pos or 
                    % neg slope instances (if det is incr or decr)

if t==0
    detJ = 1;
else
    J = reshape(Y(157:300),12,12)';
    detJ = det(J);
end

end