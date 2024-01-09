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
disp('Choose a folder containing one or more tanks that you wish to save.')
if myDir == 0
    disp("Select a tank to start")
    return
end
savDir = uigetdir('','Choose where you want to save the .mat(s).'); %gets directory%
disp('Choose a save location.')
if savDir == 0
    disp("Select a valid save directory")
    return
end

tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
numFiles = length(myFiles);
totFiles = numFiles; % variable to track how many files actually get saved
for i = 1:numFiles
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
    
end

disp("Successfully extracted and saved tank data to .mat files")
fprintf("Files saved: %d\n",totFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);