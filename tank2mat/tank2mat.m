%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     __                          __         ______                             __      
%    /  |                        /  |       /      \                           /  |     
%   _$$ |_     ______   _______  $$ |   __ /$$$$$$  | _____  ____    ______   _$$ |_   
%  / $$   |   /      \ /       \ $$ |  /  |$$____$$ |/     \/    \  /      \ / $$   |  
%  $$$$$$/    $$$$$$  |$$$$$$$  |$$ |_/$$/  /    $$/ $$$$$$ $$$$  | $$$$$$  |$$$$$$/   
%    $$ | __  /    $$ |$$ |  $$ |$$   $$<  /$$$$$$/  $$ | $$ | $$ | /    $$ |  $$ | __ 
%    $$ |/  |/$$$$$$$ |$$ |  $$ |$$$$$$  \ $$ |_____ $$ | $$ | $$ |/$$$$$$$ |  $$ |/  |
%    $$  $$/ $$    $$ |$$ |  $$ |$$ | $$  |$$       |$$ | $$ | $$ |$$    $$ |  $$  $$/ 
%     $$$$/   $$$$$$$/ $$/   $$/ $$/   $$/ $$$$$$$$/ $$/  $$/  $$/  $$$$$$$/    $$$$/   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                      
% Created by Brandon L. Oliver, M.A.
% Turns TDT tanks into .mat files freeing up storage space and speeding up
% analyses significantly. For instructions, check out the README.

myDir = uigetdir('','Choose the tank(s) you want to save.'); %gets directory%
if myDir == 0
    disp("Select a tank to start")
    return
end
savDir = uigetdir('','Choose where you want to save the .mat(s).'); %gets directory%
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
numFiles = length(myFiles);
totFiles = numFiles; % variable to track how many files actually get saved
LOAD_BAR = waitbar(0,'1 %','Name','Extracting tank files...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(LOAD_BAR,'canceling',0)
for i = 1:numFiles
    if getappdata(LOAD_BAR,'canceling')
        break
    end
    BLOCKPATH = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(BLOCKPATH);
    newfilename = strcat(name,'.mat');
    file_pathname = fullfile(savDir,newfilename);
    if exist(file_pathname,"file") % checks if the file exists in savDir and skips if it does
        fprintf("%s already exists...skipping\n",newfilename)
        totFiles = totFiles - 1;
        continue
    end
    fprintf("Extracting tank %d of %d...\n",i,numFiles)
    data = TDTbin2mat(BLOCKPATH,'TYPE',({'epocs','streams'}));

    disp("Saving...")
    save(file_pathname,"data")
    disp("Done.")
    waitbar(i/numFiles,LOAD_BAR,sprintf('Progress: %d %%',floor(i/numFiles*100)));
    pause(0.1)
end
delete(LOAD_BAR)
disp("Successfully extracted and saved tank data to .mat files")
fprintf("Files saved: %d\n",totFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);