function tempOops(x0, direc)
    allFiles= dir(direc);
    for i= 1:length(allFiles)
        if ~isdir(allFiles(i).name)
            make_ref_images(direc, allFiles(i).name, x0, 1, parameters(200));
            fprintf('done with %s\n',allFiles(i).name);
        end
    end
end