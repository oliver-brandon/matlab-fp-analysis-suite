clear
myDir = uigetdir('','Choose the mat file(s) you want to convert to json.'); %gets directory%
disp('Choose a folder containing one or more .mat files to convert to .json')
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
savDir = uigetdir('','Choose where you want to save the json file(s).'); %gets directory%
disp('Choose a destination folder to save the .json files')
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all mat files in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);

for i = 1:numFiles
    fprintf('Converting file %d of %d\n',i,numFiles)
    FILEPATH = fullfile(myDir,myFiles(i).name);
    load(FILEPATH);
    [~,name,~] = fileparts(FILEPATH);
    jsonFileName = strcat(savDir, '/', name, '.json');
    if exist(jsonFileName,"file")
        fprintf("%s already exists...skipping\n",name)
        continue
    end
    jsonStr = jsonencode(data,'PrettyPrint',true);
    fid = fopen(jsonFileName, 'w');
    fprintf(fid, '%s', jsonStr);
    fclose(fid);
    
end
toc
fprintf('Files saved to: %s\n',savDir)
NERD_STATS(toc,numFiles);






