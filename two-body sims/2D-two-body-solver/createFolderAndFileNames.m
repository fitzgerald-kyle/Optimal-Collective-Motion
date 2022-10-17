function [folderName, fileName]= createFolderAndFileNames(x0, x1, p0, muLabel)
% Determines the names of the folder and file corresponding to the
% particular BVP being solved.

    x0str= 'x0_(';
    x1str= '_x1_(';
    fileName= 'p0_(';
    for i= 1:length(x0)
        x0str= [x0str, num2str(x0(i)), ' '];
        x1str= [x1str, num2str(x1(i)), ' '];
        fileName= [fileName, num2str(round(p0(i),2)), ' '];
    end
    
    x0str(end)= ')';
    x1str(end)= ')';
    folderName= [x0str, x1str];
    fileName(end)= ')';
    fileName= [fileName, '_mu=', num2str(muLabel), '.csv'];
end