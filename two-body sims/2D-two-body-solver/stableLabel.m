function appendage= stableLabel(CSVData)
%    if ~contains(fileName, '_UNSTABLE') && ~contains(fileName, '_STABLE=')
        %CSVData= csvread([csvDir, '/', fileName]);
        mu= CSVData(:,1);
        conjPts= CSVData(:,end);
        begMu= [];
        endMu= [];
        for j= 1:length(conjPts)
            if conjPts(j)==0 && (j==1 || (j>1 && conjPts(j-1)>0))
                begMu(end+1)= mu(j);
            end
            if conjPts(j)==0 && (j==length(conjPts) || ...
                    (j<length(conjPts) && conjPts(j+1)>0))
                endMu(end+1)= mu(j);
            end
        end
        appendage= [];
        for j= 1:length(begMu)
            appendage= [appendage, num2str(begMu(j)),'-',num2str(endMu(j)),'_'];
        end
        if isempty(appendage)
            appendage= '_UNSTABLE';
        else
            appendage= ['_STABLE=', appendage(1:end-1)];
        end
    %end
end