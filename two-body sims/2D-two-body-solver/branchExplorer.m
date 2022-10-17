function branchExplorer(x0, x1, folderName, fileName)
    CSVData= csvread([folderName,'/',fileName]);
    
    for i= 2:length(CSVData(:,1))
        if CSVData(i,8)==1 && CSVData(i-1,8)==0
            [soln, nmaxErr]= constrainedAngMomSolver(x0, x1, ...
                CSVData(i,2:7), CSVData(i,1), 200);
            if nmaxErr
                continue
            end
            
            idxOfConj= 0; j= 1;
            while idxOfConj==0
                if soln.t(j) >= soln.tconj
                    idxOfConj= j;
                end
                j= j+1;
            end
            [U,D,~]= svd(soln.J(:,:,idxOfConj));
            disp(D(end,end));
            dp0= U(:,end)';
            
            dp0= dp0/norm(CSVData(i,2:7)-CSVData(i-1,2:7));
            
            mu1= CSVData(i,1) + (CSVData(i,1)-CSVData(i-1,1));
            mu2= CSVData(i,1) - (CSVData(i,1)-CSVData(i-1,1));
            
            pertSolver(x0, x1, CSVData(i,2:7)+dp0, mu1:-.01:0, 0);
            pertSolver(x0, x1, CSVData(i,2:7)-dp0, mu1:-.01:0, 0);
            pertSolver(x0, x1, CSVData(i,2:7)+dp0, mu2:-.01:0, 0);
            pertSolver(x0, x1, CSVData(i,2:7)-dp0, mu2:-.01:0, 0);
        end
    end
end