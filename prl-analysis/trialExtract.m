% clearvars -except data
clear; clc; close all;
withinSession = 1; % 1 = yes, 0 = no
myDir = uigetdir('/Users/brandon/personal-drive/prl/GrabDA',"Select a folder containing one or more files"); 

if myDir == 0
    disp("Select a folder containing one or more files")
    return
end
myFiles = dir(myDir);
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
myFiles = myFiles(endsWith({myFiles.name},{'.mat'}));
numFiles = length(myFiles);
allTrials = [];
for x = 1:numFiles
    fprintf('Loading file %d of %d...\n',x,length(myFiles))
    filename = fullfile(myDir, myFiles(x).name);
    load(filename)
    [~,name,~] = fileparts(filename);
    brokenID = strsplit(name,'_');
    animalID = char(brokenID{1});
    % get animal ID and remove any characters from the ID number
    animalID = str2double(regexp(animalID, '\d+', 'match'));
    % task = char(brokenID{2});
    % treatment = char(brokenID{3});
    if withinSession == 0
        if isfield(data.epocs, 'St1_')
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
        elseif isfield(data.epocs, 'St2_')
            cue = data.epocs.St2_.onset;
            cRew = data.epocs.cRewC.onset;
            cNoRew = data.epocs.cNoRewC.onset;
            iRew = data.epocs.iRewC.onset;
            iNoRew = data.epocs.iNoRewC.onset;
        end
    elseif withinSession == 1
        cue = data.epocs.St1_.onset;
        cRew = data.epocs.cRewA.onset;
        cNoRew = data.epocs.cNoRewA.onset;
        iRew = data.epocs.iRewA.onset;
        iNoRew = data.epocs.iNoRewA.onset;
    else
        disp('Indicate session type')
        return
    end

    for i = 1:length(cRew)
        if cRew(i,1) == 0
            cRew(1,1:3) = nan;
            continue
        end
        cRew(i,2) = 1;
        cRew(i,3) = 1;
    end

    for ii = 1:length(cNoRew)
        if cNoRew(ii,1) == 0
            cNoRew(1,1:3) = nan;
            continue
        end
        cNoRew(ii,2) = 1;
        cNoRew(ii,3) = 0;
    end

    for iii = 1:length(iRew)
        if iRew(iii,1) == 0
            iRew(1,1:3) = nan;
            continue
        end
        iRew(iii,2) = 0;
        iRew(iii,3) = 1;
    end

    for iiii = 1:length(iNoRew)
        if iNoRew(iiii,1) == 0
            iNoRew(1,1:3) = nan;
            continue
        end
        iNoRew(iiii,2) = 0;
        iNoRew(iiii,3) = 0;
    end
    
    trials = vertcat(cRew, cNoRew, iRew, iNoRew);
    trials = sortrows(trials, 1);

    % if length(trials) > 30
    %     trials = trials(1:30,:);
    % end
    nanRows = any(isnan(trials),2);
    trials = trials(~nanRows,:);

    % replace the onset times with the animalID
    for jj = 1:length(trials)
        trials(jj,1) = animalID;
    end
    
    

    % insert column between 1 and 2 with trial numbers
    trials = horzcat(trials(:,1), (1:length(trials))', trials(:,2:3));
    allTrials = vertcat(allTrials, trials);
end

% convert allTrials to a table and add column headers
allTrials = array2table(allTrials, 'VariableNames', {'id', 'trial', 'choice', 'reward'});
% save the table to a .csv file
% writetable(allTrials, 'allTrialsSfRev1VEH.csv')