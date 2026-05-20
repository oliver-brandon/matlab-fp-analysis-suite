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
infusePRL = 0; %  1=infuse all PRL epocs, 0=no epoc infusion
numFibers = 2; % 1=single fiber, 2=dual fiber
swapOnOff = 0; % 1=Swaps TTL onsets with offsets, 0=no swap
myDir = uigetdir(pwd,'Choose the tank(s) you want to save.'); %gets directory%
disp('Choose a folder containing one or more tanks that you wish to save.')
if myDir == 0
    disp("Select a tank to start")
    return
end
savDir = uigetdir(myDir,'Choose where you want to save the .mat(s).'); %gets directory%
disp('Choose a save location.')
if savDir == 0
    disp("Select a valid save directory")
    return
end

tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(~endsWith({myFiles.name}, {'.pdf','.m','.prism','.xlsx'}));
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
    if isfield(data.streams, 'Fi1r')
        data.streams = rmfield(data.streams, {'Fi1r','Fi1d'});
    else
        disp('')
    end
    if isfield(data.streams, 'Fi2r')
        data.streams = rmfield(data.streams, {'Fi2r','Fi2d'});
    else
        disp('')
    end

    if infusePRL == 1
        if numFibers == 1
            disp('Infusing PRL related epocs...')
            data = prl_epocs(data);
            disp('Done.')
        elseif numFibers == 2
            if isfield(data.epocs,'St1_')
                TTLs = 1;
            elseif isfield(data.epocs,'St2_')
                TTLs = 2;
            else
                disp('File is missing TTLs')
            end
            disp('Infusing PRL related epocs...')
            data = prl_df_epocs(data,TTLs);
            disp('Done.')
        else
            disp('Choose a valid number of fibers.')
            break
        end
    else
        disp('')
    
    end
    if swapOnOff == 1
            if isfield(data.epocs, 'IL1_')
                data.epocs.St1_.onset = data.epocs.St1_.offset;
                data.epocs.CL1_.onset = data.epocs.CL1_.offset;
                data.epocs.IL1_.onset = data.epocs.IL1_.offset;
            end
            if isfield(data.epocs, 'aRL_')
                data.epocs.aRw_.onset = data.epocs.aRw_.offset;
                data.epocs.aRL_.onset = data.epocs.aRL_.offset;
                data.epocs.aLL_.onset = data.epocs.aLL_.offset;
            end
            if isfield(data.epocs, 'IL2_')
                data.epocs.St2_.onset = data.epocs.St2_.offset;
                data.epocs.CL2_.onset = data.epocs.CL2_.offset;
                data.epocs.IL2_.onset = data.epocs.IL2_.offset;
            end
            if isfield(data.epocs, 'bRL_')
                data.epocs.bRw_.onset = data.epocs.bRw_.offset;
                data.epocs.bRL_.onset = data.epocs.bRL_.offset;
                data.epocs.bLL_.onset = data.epocs.bLL_.offset;
            end
    else
        disp('')
    
    end

    disp("Saving...")
    save(file_pathname,"data")
    disp("Done.")

    
    
end

disp("Successfully extracted and saved tank data to .mat files")
fprintf("Files saved: %d\n",totFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);