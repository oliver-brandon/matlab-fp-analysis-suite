clear;close all;

toSave = 0; % set == 1 if want to overwrite the mat files
numTrials = 30;

dataDir = uigetdir('/Users/brandon/_Lab');
if dataDir == 0
    disp('Select a folder containing one or more .mat files to start');
    return
end

tic

dataFiles = dir(dataDir);
dataFiles = dataFiles(~startsWith({dataFiles.name},{'.','..','._'}));
dataFiles = dataFiles(endsWith({dataFiles.name},{'.mat'}));
numFiles = numel(dataFiles);

IDs = cell(numFiles, 1);
phaseList = cell(numFiles, 1);
treatList = cell(numFiles, 1);
percentCorrect = zeros(numFiles, 1);
perseveration = [];
errorProbability = [];
for i = 1:numFiles
    filepath = fullfile(dataDir, dataFiles(i).name);
    try
        DATA = load(filepath);  % Load the data into a separate variable
    catch exception
        disp(['Error loading file: ' filepath]);
        disp(exception.message);
        continue;  % Skip to the next file if loading fails
    end
    data = DATA.data; % sets the top structure back to 'data'
    
    [~, name, ~] = fileparts(filepath);
    brokenID = strsplit(name, '_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    phaseList{i} = cellstr(strtrim(brokenID{2}));
    treatList{i} = cellstr(strtrim(brokenID{3}));


    if isfield(data.streams, 'x405A')
        cue = data.epocs.St1_.onset;
        cRew = data.epocs.cRewA.onset;
        cNoRew = data.epocs.cNoRewA.onset;
        iRew = data.epocs.iRewA.onset;
        iNoRew = data.epocs.iNoRewA.onset;
        [session_identifiers,lever_session_ts,trial_number,trial_name] = ...
            sessionArraySort(cue,cRew,cNoRew,iRew,iNoRew);
        if ~isfield(data.epocs, 'CL1_')
            correct = 0;
        else
            correct = height(data.epocs.CL1_.onset);
        end
        if ~isfield(data.epocs, 'IL1_')
            incorrect = 0;
        else
            incorrect = height(data.epocs.IL1_.onset);
        end
    elseif isfield(data.streams, 'x405C')
        cue = data.epocs.St2_.onset;
        cRew = data.epocs.cRewC.onset;
        cNoRew = data.epocs.cNoRewC.onset;
        iRew = data.epocs.iRewC.onset;
        iNoRew = data.epocs.iNoRewC.onset;
        [session_identifiers,lever_session_ts,trial_number,trial_name] = ...
            sessionArraySort(cue,cRew,cNoRew,iRew,iNoRew);
        if session_identifiers(1,2) == 0 && session_identifiers(2,2) == 0
            session_identifiers = session_identifiers(2:end,:);
        end
        if ~isfield(data.epocs, 'CL2_')
            correct = 0;
        else
            correct = height(data.epocs.CL2_.onset);
        end
        if ~isfield(data.epocs, 'IL2_')
            incorrect = 0;
        else
            incorrect = height(data.epocs.IL2_.onset);
        end

    end
    if correct + incorrect < 30
        percentCorrect(i,1) = nan;
    else
        percentCorrect(i,1) = correct/numTrials;
    end
    data.behavior.perCor = percentCorrect(i,1);
    %% Perseveration %%
    [totPersev] = prlPerseveration(session_identifiers);
    perseveration(end+1,1) = totPersev(1,1);
    %% Error Probability %%
    [data, errorProbLeverTS, errorProbCueTS] = errorProbExtract(...
    data,...
    session_identifiers, ...
    1, ...
    1 ...
    );
    errorProbability(end+1,1) = height(errorProbLeverTS) / sum(session_identifiers(:,2) == 1);
    %%



end
IDs = cell2table(IDs,'VariableNames',{'ID/Sex'});
phaseList = cell2table(phaseList,'VariableNames',{'Phase'});
treatList = cell2table(treatList,'VariableNames',{'Treatment'});
percentCorrect = array2table(percentCorrect,'VariableNames',{'Percent Correct'});
percentCorrect = horzcat(IDs,phaseList,treatList,percentCorrect);
percentCorrect = sortrows(percentCorrect,{'Phase','Treatment'},{'ascend','descend'});
perseveration = array2table(perseveration,'VariableNames',{'Total Perseveration'});
perseveration = horzcat(IDs,phaseList,treatList,perseveration);
perseveration = sortrows(perseveration,{'Phase','Treatment','ID/Sex'},{'ascend','descend','ascend'});
errorProbability = array2table(errorProbability,'VariableNames',{'p(error)'});
errorProbability = horzcat(IDs,phaseList,treatList,errorProbability);
errorProbability = sortrows(errorProbability,{'Phase','Treatment','ID/Sex'},{'ascend','descend','ascend'});

prl_behavior.acq1 = percentCorrect(strcmp(percentCorrect.Phase,'Acq1'), :);
prl_behavior.acq2 = percentCorrect(strcmp(percentCorrect.Phase,'Acq2'), :);
prl_behavior.rev1 = percentCorrect(strcmp(percentCorrect.Phase,'Rev1'), :);
prl_behavior.rev2 = percentCorrect(strcmp(percentCorrect.Phase,'Rev2'), :);
prl_behavior.rev3 = percentCorrect(strcmp(percentCorrect.Phase,'Rev3'), :);
if toSave == 1
    save(filepath,"data")
else
    disp('')
end

toc
NERD_STATS(toc,numFiles);
