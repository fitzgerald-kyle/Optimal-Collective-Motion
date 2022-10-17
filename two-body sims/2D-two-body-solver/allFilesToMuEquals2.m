function allFilesToMuEquals2(x0, x1, csvDir)
    allFiles= dir(csvDir);
    for i= 1:length(allFiles)
        if ~isfolder(allFiles(i).name)
            CSVData= csvread([csvDir, '/', allFiles(i).name]);
            pertSolver(x0, x1, CSVData(end,2:7), CSVData(1,1):.01:2, false, true);
        end
    end
end