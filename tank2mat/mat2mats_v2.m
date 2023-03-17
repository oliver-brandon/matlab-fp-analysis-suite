%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            _       ______                             __              
%                          /  |     /      \                           /  |             
%  _____  ____    ______   _$$ |_   /$$$$$$  | _____  ____    ______   _$$ |_    _______ 
% /     \/    \  /      \ / $$   |  $$____$$ |/     \/    \  /      \ / $$   |  /       |
% $$$$$$ $$$$  | $$$$$$  |$$$$$$/    /    $$/ $$$$$$ $$$$  | $$$$$$  |$$$$$$/  /$$$$$$$/ 
% $$ | $$ | $$ | /    $$ |  $$ | __ /$$$$$$/  $$ | $$ | $$ | /    $$ |  $$ | __$$      \ 
% $$ | $$ | $$ |/$$$$$$$ |  $$ |/  |$$ |_____ $$ | $$ | $$ |/$$$$$$$ |  $$ |/  |$$$$$$  |
% $$ | $$ | $$ |$$    $$ |  $$  $$/ $$       |$$ | $$ | $$ |$$    $$ |  $$  $$//     $$/ 
% $$/  $$/  $$/  $$$$$$$/    $$$$/  $$$$$$$$/ $$/  $$/  $$/  $$$$$$$/    $$$$/ $$$$$$$/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                      
% Created by Brandon L. Oliver, M.A.
% Separates all mat files containing two animal's stream/epoc data in a
% directory into two separate mat files (to a user selected directory) for 
% easier grouping during analysis. Mat file name must be formatted like the
% following:
% ID1_Task1_Treatment1_ID2_Task2_Treatment2 or
% ID1_Task1_Treatment1_Empty_NA_NA if there is only one animal in the file.
% Fill task/treatment fields with "NA"
clear
VERSION = "2.0";
fprintf("VERSION: %s\n",VERSION)
myDir = uigetdir('','Choose the mat file(s) you want to save.'); %gets directory%
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
savDir = uigetdir('','Choose where you want to save the separated mat file(s).'); %gets directory%
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
totFiles = numFiles*2; % variable to track how many files actually get saved
LOAD_BAR = waitbar(0,'1','Name','Splitting mat files...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(LOAD_BAR,'canceling',0)
for i = 1:numFiles
    FILEPATH = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(FILEPATH);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    animalIDA = char(brokenID{1});
    animalIDC = char(brokenID{4});
    emptylogicA = strcmp(animalIDA,emptyID);
    emptylogicC = strcmp(animalIDC,emptyID);
    load(FILEPATH);
    data.streams = rmfield(data.streams, {'Fi2r','Fi2d','Fi1r','Fi1d'});
    for subjects = 1:2
        if subjects == 1
            if emptylogicA == 1
                disp("Stream A is empty")
                totFiles = totFiles - 1;
                continue
            elseif emptylogicA == 0
                data.streams = rmfield(data.streams, {'x405C','x465C'});
                if isfield(data.epocs,'IL2_') == 1
                    data.epocs = rmfield(data.epocs, {'Pe2_','St2_','IL2_'});
                end
                if isfield(data.epocs,'CL2_') == 1
                    data.epocs = rmfield(data.epocs, {'CL2_'});
                end
                taskA = char(brokenID{2});
                treatmentA = char(brokenID{3});
                
                matchA = regexp(animalIDA, '^[a-zA-Z]{2}','match');
                if ~isempty(matchA)
                    animalIDA = animalIDA(3:end);
                else
                    idxA = isletter(animalIDA);
                    animalIDA(idxA) = regexprep(animalIDA(idxA), '[a-zA-Z]', '');
                end
                
                newfilenameA = strcat(animalIDA,'_',taskA,'_',treatmentA,'.mat');
                file_pathnameA = fullfile(savDir,newfilenameA);
                if exist(file_pathnameA,"file") % checks if the file exists in savDir and skips if it does
                    fprintf("%s already exists...skipping\n",newfilenameA)
                    totFiles = totFiles - 1;
                    continue
                end
                disp("Saving stream A...")
                save(file_pathnameA,"data")
                disp("Done.")
                clear data
            end
        elseif subjects == 2
            load(FILEPATH);
            data.streams = rmfield(data.streams, {'Fi2r','Fi2d','Fi1r','Fi1d'});
            if emptylogicC == 1
                disp("Stream C is empty")
                totFiles = totFiles - 1;
                continue
            elseif emptylogicC == 0
                data.streams = rmfield(data.streams, {'x405A','x465A'});
                if isfield(data.epocs,'IL1_') == 1
                    data.epocs = rmfield(data.epocs, {'Pe1_','St1_','IL1_'});
                end
                if isfield(data.epocs,'CL1_') == 1
                    data.epocs = rmfield(data.epocs, {'CL1_'});
                end
                taskC = char(brokenID{5});
                treatmentC = char(brokenID{6});
                matchC = regexp(animalIDC, '^[a-zA-Z]{2}','match');
                if ~isempty(matchC)
                    animalIDC = animalIDC(3:end);
                else
                    idxC = isletter(animalIDC);
                    animalIDA(idxC) = regexprep(animalIDC(idxC), '[a-zA-Z]', '');
                end
                newfilenameC = strcat(animalIDC,'_',taskC,'_',treatmentC,'.mat');
                file_pathnameC = fullfile(savDir,newfilenameC);
                if exist(file_pathnameC,"file") % checks if the file exists in savDir and skips if it does
                    fprintf("%s already exists...skipping\n",newfilenameC)
                    totFiles = totFiles - 1;
                    continue
                end
                disp("Saving stream C...")
                save(file_pathnameC,"data")
                disp("Done.")
                clear data
            end
        end
    end
    waitbar(i/numFiles,LOAD_BAR,sprintf('Progress: %d %%',floor(i/totFiles*100)));
    pause(0.1)
end
delete(LOAD_BAR)
disp("Successfully separated and saved streams to individual .mat files")
fprintf("Files saved: %d\n",totFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,totFiles);