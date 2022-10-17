function tconjpts = is_stable(tconj,ie,output)
% Ensure that conj points found are true conj points by requiring det(J)
% to pass a certain threshold before finding det(J)=0
% This prevents Matlab from mistaking small det(J) as det(J) = 0

% Find the index of the first instance where det(J) passes threshold
ithreshold = find(ie==2,1,'first');

% Find regular conjugate points (i.e. det(J) passes through 0)
tconjpts = [];
if (~isempty(ithreshold))
    % Find indices where det(J) is marked as 0
    iconj = find(ie==1);
    % Only keep indices of conj pts past the threshold index
    iconj = iconj(iconj>ithreshold);
    % Get regular conjugate times
    tconjpts = tconj(iconj);
end

% Find indices where det(J) grazes 0 (i.e. abs(det(J)) passes below the
%  threshold without crossing 0)
tgraze = [];
grazeTol= 1e-50;
if (~isempty(ithreshold))
    for i= 2:length(output.detJ)-1
        if i > ithreshold && abs(output.detJ(i)) < grazeTol ...
                && abs(output.detJ(i)) < abs(output.detJ(i-1)) ...
                && abs(output.detJ(i)) < abs(output.detJ(i+1)) ...
                && sign(output.detJ(i)) == sign(output.detJ(i-1)) ...
                && sign(output.detJ(i)) == sign(output.detJ(i+1))
            tgraze(end+1)= output.t(i);
            tgraze= tgraze(:);
        end
    end
end
%{
k=1;
if (~isempty(ithreshold))
    % Find indices where det(J) crosses threshold
    igraze = find(ie==2);
    % Only keep indices past the initial threshold crossing
    igraze = igraze(igraze>ithreshold);
    % Check if any indices are sequential (i.e. det(J) passes
    % through threshold twice without crossing 0)
    igraze = igraze(diff(igraze)==1);
    % Check if det(J) is below threshold between sequential indices
    for j=1:length(igraze)
        d1 = abs(output.detJ(find(output.t<tconj(igraze(j)),1,'last')));
        d2 = abs(output.detJ(find(output.t>tconj(igraze(j)),1,'first')));               
        if d2<d1
            tgraze(k) = tconj(igraze(j));
            k=k+1;
        end
    end
    if (~isempty(tgraze))
        % Find indices of times tgraze within tconj
        igraze = find(ismember(tconj,tgraze));
        % Approximate grazing conjugate points by averaging times when det(J)
        % passes through threshold
        tgraze = 0.5*(tconj(igraze)+tconj(igraze+1));
    end
end
%}
% Combine tgraze and tconj
tconjpts = sort([tgraze; tconjpts]);
    
end