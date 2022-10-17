function allFilesToMuEquals2(x0, x1, csvDir)
    allFiles= dir(csvDir);
    for i= 1:length(allFiles)
        if ~isdir(allFiles(i).name)
            CSVData= csvread([csvDir, '/', allFiles(i).name]);
            if CSVData(end,1) < 2
                pertSolver(x0, x1, CSVData(end,2:7), CSVData(end,1):.01:2, 0);
            end
        end
    end
end