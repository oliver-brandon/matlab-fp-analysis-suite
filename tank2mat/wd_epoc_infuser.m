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
myDir = uigetdir('','Choose the mat file(s) you want to infuse.'); %gets directory%
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
tsDir = uigetdir('','Choose where the timetsamp files are located.'); %gets directory%
if tsDir == 0
    disp("Select a valid directory")
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

for i = 1:numFiles

    FILEPATH = fullfile(myDir,myFiles(i).name);
    [~,filename,~] = fileparts(FILEPATH);
    load(FILEPATH);
    data = wd_epocs(filename,tsDir,data);
    save(FILEPATH,"data")
   
end
disp("Successfully infused epocs into mat files")
fprintf("Files infused: %d\n",numFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);