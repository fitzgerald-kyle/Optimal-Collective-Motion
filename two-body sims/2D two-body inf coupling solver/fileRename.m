function fileRename(directory, appendage)
    allFiles= dir(directory);
    for i= 1:length(allFiles)
        if ~isdir(allFiles(i).name)
            copyfile([directory,'/',allFiles(i).name], ...
                [directory,'/',allFiles(i).name(1:end-13),appendage]);
            delete([directory,'/',allFiles(i).name]);
        end
    end
end