function output = find_detJ(output)

% This function computes the determinant of J and adds it to the output
% variable
% Function Inputs:
%   output : structure containing the current solution
% Function Outputs:
%   tc: time of first conjugate point

% Initialize detJ
output.detJ = zeros(1,length(output.t));

% Compute detJ
for i=1:length(output.t)
    output.detJ(i) = det(output.J(:,:,i));
end

end