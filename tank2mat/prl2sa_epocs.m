% Created by Brandon L. Oliver, M.A.
% Infuses all mat files containing one animal's stream/epoc data in a
% directory with epocs from a given function (to a user selected directory).
% For instructions, see "README"
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dualFiber = 1; % 1 = dualFiber file, 0 = singleFiber file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('','Choose the mat file(s) you want to infuse.'); %gets directory%
disp('Choose a folder containing one or more mat files that you wish to infuse with prl epochs.')
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
savDir = uigetdir('','Choose where you want to save the infused mat file(s).'); %gets directory%
disp('Choose a folder to save the infused mat files to.')
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);

for i = 1:numFiles

    FILEPATH = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(FILEPATH);
    
    
    load(FILEPATH);
    if dualFiber == 0
        if isfield(data.epocs, 'aRw_') || isfield(data.epocs, 'aRL_')
            fprintf('%s already infused...skipping (%d of %d)\n',name,i,numFiles)
        elseif isfield(data.epocs, 'St1_') || isfield(data.epocs, 'IL1_')
            data.epocs.aRw_ = data.epocs.St1_;
            data.epocs.aRL_ = data.epocs.IL1_;
            data.epocs.aLL_ = data.epocs.CL1_;

            fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
        end
        
    elseif dualFiber == 1
        if isfield(data.epocs, 'aRw_') || isfield(data.epocs, 'aRL_')
            fprintf('%s already infused...skipping (%d of %d)\n',name,i,numFiles)
        elseif isfield(data.epocs, 'St1_') || isfield(data.epocs, 'IL1_')
            data.epocs.aRw_ = data.epocs.St1_;
            data.epocs.aRL_ = data.epocs.IL1_;
            data.epocs.aLL_ = data.epocs.CL1_;

            fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
        end
    save(FILEPATH,"data")
    
    end
end
toc
disp("Successfully infused epocs into mat files")
fprintf("Files infused: %d\n",numFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);