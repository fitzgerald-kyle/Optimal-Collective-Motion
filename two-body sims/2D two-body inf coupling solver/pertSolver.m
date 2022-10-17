function solns= pertSolver(x0, x1, p0, figs)
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

    [folderName, fileName]= createFolderAndFileNames(x0, x1, p0);
    
    [newWrite, altName]= isRepeatedSoln(folderName, p0);
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

    if fileExists
        [solns, nmaxErr]= constrainedAngMomSolver(x0, x1, CSVData(1:6), 200);
        CSVData(1:6)= solns.p0;
        CSVData(7)= length(solns.tconj);
    else
        [solns, nmaxErr]= constrainedAngMomSolver(x0, x1, p0, 1000);

        if ~nmaxErr
            [newWrite, altName]= isRepeatedSoln(folderName, solns.p0);
            if ~newWrite
                fileExists= 1;
                fileName= altName;
                CSVData= csvread([folderName, '/', fileName]);
                fprintf('Warning: using %s instead\n', fileName);
            end
        end

        CSVData= [solns.p0, 0];

        if ~nmaxErr && ~isempty(solns.tconj) % if conj points exist, store in CSV
            CSVData(7) = length(solns.tconj);
        end
    end

    if nmaxErr
        fprintf('Warning: solve stopped\n');
    end
    
    if isempty(CSVData)
        return
    end
    
    [~, fileName]= createFolderAndFileNames(x0, x1, CSVData(1:6));
%        if fileExists
%            delete([folderName,'/',fileName]);
%        end
    
    appendage= stableLabel(CSVData);
    fileName= [fileName(1:end-4), appendage, '.csv'];
    
    % write data to CSV file
    fprintf('writing to CSV...\n');
    writematrix(CSVData, [folderName, '/', fileName]);
    
    fprintf('making images...\n'); 
    make_ref_images(folderName, fileName, x0, 1, parameters(200));
        
    if ~figs
        return
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
    plot(solns.t(2:end), solns.detJ(2:end));
    plot(0:1, [0,0], '-k'); % solid black line at det(J)=0
    xlabel('t'); ylabel('det(J(t))');
    %title(['Red = small \mu, violet = large \mu; ', ...
    %    num2str(minMuArg), ' \leq \mu \leq ', num2str(maxMuArg)]);
    drawnow
    fprintf('plot of detJ finished\n');
    
%
    %% plot solution paths as mu is incremented
    figure; hold on
    xlabel('x_1')
    ylabel('x_2')
    plot(solns.x(:,1), solns.x(:,2), 'b-')
    hold on
    plot(solns.x(:,4), solns.x(:,5), 'g-')
    hold off
    %title(['\mu = ', num2str(solns{i}.mu)]);
    axis([-0.5 1 -1 1])
    drawnow
    fprintf('plot of solution paths finished\n');
%}  
%{   
    %% plot bifurcation diagram
    bifur_plot([folderName '_HEADERS']);
    fprintf('bifurcation plot finished\n');
%}
end



function [newWrite, altFile] = isRepeatedSoln(folderName, p0Soln)
    allFiles = dir(folderName);
    for i = 1:length(allFiles)
        if ~allFiles(i).isdir
            tempMat= csvread([folderName,'/',allFiles(i).name]);
            if round(tempMat(1,1:6),3) == round(p0Soln,3)
                altFile= allFiles(i).name;
                newWrite= 0;
                return
            end
        end
    end
    
    altFile= '';
    newWrite= 1;
end