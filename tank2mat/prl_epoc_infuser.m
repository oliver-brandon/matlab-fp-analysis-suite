%                                     __             __           ______                                     
%                                    /  |           /  |         /      \                                    
%  ______   ______   ______   _______$$ |____       $$/ _______ /$$$$$$  __    __  _______  ______   ______  
% /      \ /      \ /      \ /       $$      \      /  /       \$$ |_ $$/  |  /  |/       |/      \ /      \ 
%/$$$$$$  /$$$$$$  /$$$$$$  /$$$$$$$/$$$$$$$  |     $$ $$$$$$$  $$   |  $$ |  $$ /$$$$$$$//$$$$$$  /$$$$$$  |
%$$    $$ $$ |  $$ $$ |  $$ $$ |     $$ |  $$ |     $$ $$ |  $$ $$$$/   $$ |  $$ $$      \$$ |  $$ $$ |  $$/ 
%$$$$$$$$/$$ |__$$ $$ \__$$ $$ \_____$$ |  $$ |     $$ $$ |  $$ $$ |    $$ \__$$ |$$$$$$  $$ \__$$ $$ |      
%$$       $$    $$/$$    $$/$$       $$ |  $$ ______$$ $$ |  $$ $$ |    $$    $$//     $$/$$    $$/$$ |      
% $$$$$$$/$$$$$$$/  $$$$$$/  $$$$$$$/$$/   $$/      $$/$$/   $$/$$/      $$$$$$/ $$$$$$$/  $$$$$$/ $$/       
%         $$ |                               $$$$$$/                                                         
%         $$ |                                                                                               
%         $$/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Brandon L. Oliver, M.A.
% Infuses all mat files containing one animal's stream/epoc data in a
% directory with epocs from a given function (to a user selected directory).
% For instructions, see "README"
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dualFiber = 0; % 1 = dualFiber file, 0 = singleFiber file
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
        if isfield(data.epocs, 'cRewA') || isfield(data.epocs, 'cRewC')...
                || isfield(data.epocs, 'iNoRewA') || isfield(data.epocs, 'iNoRewC')
            fprintf('%s already infused...skipping (%d of %d)\n',name,i,numFiles)
        else
            fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
        end
        
        data = prl_epocs(data);
    elseif dualFiber == 1
        if isfield(data.epocs, 'cRewA') || isfield(data.epocs, 'cRewC') ||...
                isfield(data.epocs, 'iNoRewA') || isfield(data.epocs, 'iNoRewC')
            fprintf('%s already infused...skipping (%d of %d)\n',name,i,numFiles)
        else
            fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
        end

        if isfield(data.epocs,'St1_')
            TTLs = 1;
        elseif isfield(data.epocs,'St2_')
            TTLs = 2;
        else
            disp('File is missing TTLs')
        end
        fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
        data = prl_df_epocs(data,TTLs);
    end
    save(FILEPATH,"data")
    
end

disp("Successfully infused epocs into mat files")
fprintf("Files infused: %d\n",numFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);