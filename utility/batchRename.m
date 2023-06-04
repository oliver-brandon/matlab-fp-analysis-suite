clear all;clc;

oldName = '_NA';
newName = '_JZL18';

myDir = uigetdir('/Users/brandon/DA_PRL','Choose the .mat files you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
tic

myFiles = dir(myDir); %gets all mat files in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);

for i = 1:numFiles
    old_filename = myFiles(i).name;

    new_filename = strrep(old_filename,oldName,newName);

    % Rename the file
    oldFilePath = fullfile(myDir, old_filename);
    newFilePath = fullfile(myDir, new_filename);
    movefile(oldFilePath, newFilePath);

end
toc
NERD_STATS(toc,numFiles);