function tconjpts = is_stable(tconj,output)
% Ensure that conj points found are true conj points by requiring det(J)
% to pass a certain threshold before finding det(J)=0
% This prevents Matlab from mistaking small det(J) as det(J) = 0

detJ = output.detJ;

% Find the index of the first instance where det(J) passes threshold
ithreshold = find(abs(detJ) > abs(detJ(2))*1000, 1, 'first');

% Find regular conjugate points (i.e. det(J) passes through 0)
tconjpts = [];
tgraze = [];
if ~isempty(ithreshold)
    % Get regular conjugate times, only keeping those past the threshold index
    tconjpts = tconj(tconj>output.t(ithreshold));

    imin = find(islocalmin(detJ)); % find indices of detJ local minima
    iconjmin = [];
    % keep minima that graze zero from above OR cross zero at one point
    for i= 1:length(imin)
        if detJ(imin(i)) >= 0 || ...
                (detJ(imin(i)+1)>=0 && detJ(imin(i))<0 && detJ(imin(i)-1)>=0)
            iconjmin(end+1) = imin(i);
        end
    end
    
    imax = find(islocalmax(detJ)); % find indices of detJ local maxima
    iconjmax = [];
    % keep maxima that graze zero from below OR cross zero at one point
    for i= 1:length(imax)
        if detJ(imax(i)) <= 0 || ...
                (detJ(imax(i)+1)<=0 && detJ(imax(i))>0 && detJ(imax(i)-1)<=0)
            iconjmax(end+1) = imax(i);
        end
    end
    
    compFactor = 2; % comparison factor for determining e.g. how much smaller 
                    % a min must be than previous max to be labeled as a root
                    
    if ~isempty(imin) && ~isempty(imax)                
        idxOffset = imin(1)<imax(1);

        % keep minima that are past the threshold index and are at least
        % compFactor times smaller than the previous maximum
        for i= 1:length(iconjmin)
            if iconjmin(i) > ithreshold && ...
                    abs( detJ(imax(find(imin==iconjmin(i))-idxOffset)) ) > ...
                    abs(compFactor*detJ(iconjmin(i)))
                tgraze(end+1)= output.t(iconjmin(i));
            end
        end

        % keep maxima that are past the threshold index and are at least
        % compFactor times smaller in absolute value than the previous minimum
        for i= 1:length(iconjmax)
            if iconjmax(i) > ithreshold && ...
                    abs( detJ(imin(find(imax==iconjmax(i))-~idxOffset)) ) > ...
                    compFactor*abs(detJ(iconjmax(i)))
                tgraze(end+1)= output.t(iconjmax(i));
            end
        end
    end
    
    tgraze= tgraze(:);
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
        d1 = abs(detJ(find(output.t<tconj(igraze(j)),1,'last')));
        d2 = abs(detJ(find(output.t>tconj(igraze(j)),1,'first')));               
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
