myPath = uigetdir('','Choose the tank(s) you want to save.');
tic
myFiles = dir(myPath);
myFiles = myFiles(~ismember({myFiles.name},{'.','..'}));
numFiles = length(myFiles);
for i = 1:numFiles
    filepath = fullfile(myPath,myFiles(i).name);
    load(filepath);
end
toc
