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
clear
myDir = uigetdir('','Choose the mat file(s) you want to infuse.'); %gets directory%
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
savDir = uigetdir('','Choose where you want to save the infused mat file(s).'); %gets directory%
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
LOAD_BAR = waitbar(0,'1','Name','Infusing mat files...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(LOAD_BAR,'canceling',0)
for i = 1:numFiles
    if getappdata(LOAD_BAR,'canceling')
        break
    end
    FILEPATH = fullfile(myDir,myFiles(i).name);
    load(FILEPATH);
    data = prl_epocs(data);
    save(FILEPATH,"data")
    waitbar(i/numFiles,LOAD_BAR,sprintf('Progress: %d %%',floor(i/numFiles*100)));
    pause(0.1)
end
delete(LOAD_BAR)
disp("Successfully infused epocs into mat files")
fprintf("Files infused: %d\n",numFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);