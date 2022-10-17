function solns= pertSolver(x0, x1, p0, muVec, figs)
% Takes initial momentum p0 (entries 1-3 corresponding to p1, p2, p3 for 
% body 1, and entries 4-6 corresponding to p1, p2, p3 for body 2) along 
% with boundary values x0 and x1 (laid out the same way as p0) and 
% solves the constrained angular momentum optimization problem for each mu 
% in muVec. The momentum solution for each iteration is used as the
% initial momentum guess for the next iteration.
%
% Reads from CSV file 'p0_(______)_mu=______' in the folder 
% 'x0_(______)_x1_(______)' if file already exists. File has nine columns:
% the first is a list of mu values, the next six are the components of p0 
% that solve the optimization problem for those mu values, the eighth is 
% the tolerance needed to solve the BVP, and the ninth is a binary T/F 
% indicating whether or not the solution has a conjugate point. Writes 
% solutions obtained to the same file.
%
% Creates a variety of plots to visualize solutions, inflection points in
% solutions, and conjugate points in solutions. The last plot is a
% bifurcation plot for p0 as a function of changing mu. The input
% 'figs' is a 1 or 0 indicating whether to generate figures.
%
% EXAMPLE INPUT:
% perturbationSolver([0,0,0,0,0,0], [.2,0,0,.4,0,0], [1,0,1,2,0,-1], ...
%     0:.01:1, 1)

    close all

    [folderName, fileName]= createFolderAndFileNames(x0, x1, p0, muVec(1));
    
    [newWrite, altName]= isRepeatedSoln(folderName, p0, muVec(1));
     if ~newWrite
         fileExists= 1;
         fileName= altName;
         CSVData= csvread([folderName, '/', fileName]);
         fprintf('Warning: using %s instead\n', fileName); 
     elseif exist([folderName, '/', fileName], 'file') == 2
        fileExists= 1;
        CSVData= csvread([folderName, '/', fileName]);
    else
        if exist(folderName, 'dir') == 0
            mkdir(folderName);
        end
        fileExists= 0;
        CSVData= [];
    end

    muVec= round(muVec, 10);
    len= length(muVec);

    output= {};
    solns= {}; % will store all relevant data to be plotted
    lastp0= p0;
    for i= 1:length(muVec) % solve BVP for each mu in muVec
        if fileExists && ismember(muVec(i), CSVData(:,1))
            [~, loc]= ismember(muVec(i), CSVData(:,1));

            [output, nmaxErr]= constrainedAngMomSolver(x0, x1, ...
                CSVData(loc, 2:7), CSVData(loc, 1), 200);
            CSVData(loc, 2:7)= output.p0;
            CSVData(loc, 8)= length(output.tconj);
            lastp0= CSVData(loc, 2:7);
        else
            if i==1
                [output, nmaxErr]= constrainedAngMomSolver(x0, x1, p0, ...
                    muVec(i), 1000);
                
                 if ~nmaxErr
                     [newWrite, altName]= isRepeatedSoln(folderName, output.p0, muVec(1));
                     if ~newWrite
                         fileExists= 1;
                         fileName= altName;
                         CSVData= csvread([folderName, '/', fileName]);
                         fprintf('Warning: using %s instead\n', fileName);
                     end
                 end
                 
                 %if ~nmaxErr && (output.p0(2) < -1e-2 || output.p0(3) < -1e-2)
                 %    fprintf('Warning: stopped b/c either p0(2) or p0(3) was negative\n');
                 %    CSVData= [];
                 %    break
                 %end
            else            
               [output, nmaxErr]= constrainedAngMomSolver(x0, x1, ...
                    lastp0, muVec(i), 500);
            end

            if ~fileExists && i==1 && ~nmaxErr
                CSVData= [muVec(i), output.p0, 0];
            elseif ~nmaxErr && ~ismember(muVec(i), CSVData(:,1))
                CSVData= [CSVData; muVec(i), output.p0, 0];
            end

            if ~nmaxErr && ~isempty(output.tconj) % if conj points exist, store in CSV
                CSVData(CSVData(:,1) == muVec(i), 8) = length(output.tconj);
            end
        end

        if nmaxErr
            len= i-1;
            fprintf('Warning: solve stopped at mu=%.4f\n', muVec(i));
            break
        end

        solns{i,1}.t= output.t;
        solns{i,1}.x= output.x;
        solns{i,1}.p= output.p;
        solns{i,1}.detJ= output.detJ;
        solns{i,1}.tconj= output.tconj;
        solns{i,1}.inflecPts= output.inflecPts;
        solns{i,1}.mu= output.mu;

        lastp0= output.p0;

        % print an update after completion of every integer mu
        if mod(muVec(i), 1) == 0
            fprintf('mu solved: %i\n', muVec(i));
        end
    end
    
    if isempty(CSVData)
        return
    end
    
    CSVData= sortByMu(CSVData, length(CSVData(:,1))-len+1);
    muVec= muVec(1:len);
    [muVec, I]= sort(muVec);
    
    [~, fileName]= createFolderAndFileNames(x0, x1, CSVData(1,2:7), CSVData(1,1));
%        if fileExists
%            delete([folderName,'/',fileName]);
%        end
    
    appendage= stableLabel(CSVData);
    fileName= [fileName(1:end-4), appendage, '.csv'];

    tempSolns= cell(len, 1);
    for i= 1:len
        tempSolns{i}= solns{I(i)};
    end
    solns= tempSolns;
    
    % write data to CSV file
    fprintf('writing to CSV...\n');
    writematrix(CSVData, [folderName, '/', fileName]);
    
    fprintf('making images...\n'); 
    %make_ref_images(folderName, fileName, x0, 1, parameters(200));
        
    if ~figs
        return
    end

    minMuArg= muVec(1);
    maxMuArg= muVec(end);
    if maxMuArg == minMuArg
        maxMuArg= minMuArg+1e-10;
    end    

    %% Plot log_10|det(J(t))| to look for conjugate points
    %{
    figure; hold on
    for i= 1:len
        % Plot and change hue of line based on value of mu
        plot(solns{i}.t(2:end), log10(abs(solns{i}.detJ(2:end))), 'Color', ...
            hsv2rgb([(solns{i}.mu-minMuArg)/(maxMuArg-minMuArg)*0.9, 1, 1]));
    end
    xlabel('t'); ylabel('log_{10}|det(J(t))|');
    title(['Red = small \mu, violet = large \mu; ', ...
        num2str(minMuArg), ' \leq \mu \leq ', num2str(maxMuArg)]);
    drawnow
    fprintf('log plot of detJ finished\n');
    %}


    %% Plot det(J(t)) to look for conjugate points
    figure; hold on
    for i= 1:len
        plot(solns{i}.t(2:end), solns{i}.detJ(2:end), 'Color', ...
            hsv2rgb([(solns{i}.mu-minMuArg)/(maxMuArg-minMuArg)*0.9, 1, 1]));
        plot(0:1, [0,0], '-k'); % solid black line at det(J)=0
    end
    xlabel('t'); ylabel('det(J(t))');
    title(['Red = small \mu, violet = large \mu; ', ...
        num2str(minMuArg), ' \leq \mu \leq ', num2str(maxMuArg)]);
    drawnow
    fprintf('plot of detJ finished\n');

%{
    %% plot the times at which conjugate points occur as a function of mu
    figure; hold on
    for i= 1:len
        for j= 1:length(solns{i}.tconj)
            plot(solns{i}.mu, solns{i}.tconj(j), '.b', 'MarkerSize', 10);
        end
    end
    xlabel('\mu'); ylabel('t');
    title('t vs \mu for conjugate points');
    drawnow
    fprintf('plot of conj pts finished\n');
%}
%{
    %% plot inflection points vs mu
    figure; hold on
    for i= 1:len
        plot(solns{i}.mu, solns{i}.inflecPts(1), '.b');
        plot(solns{i}.mu, solns{i}.inflecPts(2), '.g');
    end
    xlabel('\mu'); ylabel('# inflection points');
    title('Number of inflection points as a function of \mu');
    drawnow
    fprintf('plot of inflection pts finished\n');
%}
%
    %% plot solution paths as mu is incremented
    figure; hold on
    xlabel('x_1')
    ylabel('x_2')
    for i= len:len%i= 1:len
        plot(solns{i}.x(:,1), solns{i}.x(:,2), 'b-')
        hold on
        plot(solns{i}.x(:,4), solns{i}.x(:,5), 'g-')
        hold off
        title(['\mu = ', num2str(solns{i}.mu)]);
        axis([-0.5 1 -1 1])
        %pause(5/len);
    end
    drawnow
    fprintf('plot of solution paths finished\n');
%}  
%{   
    %% plot bifurcation diagram
    bifur_plot([folderName '_HEADERS']);
    fprintf('bifurcation plot finished\n');
%}
end


function CSVData = sortByMu(CSVData, beginSortIdx)
% Sorts data matrix, which contains a column for mu and six columns for
% components of p0, in order of increasing mu.

    for i= beginSortIdx:length(CSVData(:,1))
        j= i;
        while j > 1 && CSVData(j-1,1) > CSVData(j,1)
            temp= CSVData(j-1,:);
            CSVData(j-1,:)= CSVData(j,:);
            CSVData(j,:)= temp;
            j= j-1;
        end
    end    
end

function [newWrite, altFile] = isRepeatedSoln(folderName, p0Soln, firstMu)
    allFiles = dir(folderName);
    for i = 1:length(allFiles)
        if ~allFiles(i).isdir
            tempMat= csvread([folderName,'/',allFiles(i).name]);
            if round(tempMat(tempMat(:,1)==firstMu,2:7),3) == round(p0Soln,3)
                altFile= allFiles(i).name;
                newWrite= 0;
                return
            end
        end
    end
    
    altFile= '';
    newWrite= 1;
end