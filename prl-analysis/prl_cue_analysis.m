clear; close all;
warning off

%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % baseline signal to include before TTL 
baseline = [-2 -1]; % baseline signal for dFF/zscore
amp_window = [0.5 1]; % time window to grab amplitude from
auc_window = [0 5];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
adjustBase = 1; % 1 = yes, 2 = no
baseAdjust = -2;
removeOutlierTrials = 0; % 1 = remove
plotOmitted = 0; % 1 = plot
removeFirstTrial = 1; % 1 = keep, 2 = remove
dualFiber = 0; % 1 = dual, 0 = single; Change ISOS/SIGNAL manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/My Drive/prl/PRL_GRABDA','Choose the .mat files you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
tic

myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
IDs = {};
treatList = {};
prl_phase = {};
omitted = struct('file', {}, 'timestamp', {}, 'trial', {}, 'signal', {});
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    fprintf('Analyzing %s (%d of %d)\n', name, i, numFiles)
    brokenID = strsplit(name,'_');
    tempID = cellstr(brokenID(1));
    tempPhase = cellstr(brokenID(2));
    tempTreat = cellstr(brokenID(3));
    IDs = vertcat(IDs,tempID);
    treatList = vertcat(treatList,tempTreat);
    prl_phase = vertcat(prl_phase,tempPhase);
    load(filename)
    if dualFiber == 1
        ISOS = 'x405C';
        SIGNAL = 'x465C';
        if isfield(data.epocs,'St1_')
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
        elseif isfield(data.epocs,'St2_')
            cue = data.epocs.St2_.onset;
            cRew = data.epocs.cRewC.onset;
            cNoRew = data.epocs.cNoRewC.onset;
            iRew = data.epocs.iRewC.onset;
            iNoRew = data.epocs.iNoRewC.onset;
        else
            disp('Cannot find epocs')
            break
        end
        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
        data.analysis.sessionID = session_identifiers;
    else
        if isfield(data.streams, 'x405A')
            ISOS = 'x405A';
            SIGNAL = 'x465A';
    
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
            [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
                cNoRew,iRew,iNoRew);
            data.analysis.sessionID = session_identifiers;
            
        elseif isfield(data.streams, 'x405C')
            ISOS = 'x405C';
            SIGNAL = 'x465C';
    
            cue = data.epocs.St2_.onset;
            cRew = data.epocs.cRewC.onset;
            cNoRew = data.epocs.cNoRewC.onset;
            iRew = data.epocs.iRewC.onset;
            iNoRew = data.epocs.iNoRewC.onset;
            [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
                cNoRew,iRew,iNoRew);
            data.analysis.sessionID = session_identifiers;
        end
    end
    %time array used for all streams%
    session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    ind = find(session_time>t,1);% find first index of when time crosses threshold
    session_time = session_time(ind:end); % reformat vector to only include allowed time
    SIGNAL_raw = data.streams.(SIGNAL).data(ind:end);
    ISOS_raw = data.streams.(ISOS).data(ind:end);
   
    %downsample streams and time array by N times%
    ISOS_raw = downsample(ISOS_raw, N);
    SIGNAL_raw = downsample(SIGNAL_raw, N);
    minStreamLength = min(length(ISOS_raw),length(SIGNAL_raw));
    ISOS_raw = ISOS_raw(1:minStreamLength);
    SIGNAL_raw = SIGNAL_raw(1:minStreamLength);
    SIGNAL_dFF = deltaFF(SIGNAL_raw);
    SIGNAL_z = zScore(SIGNAL_dFF);
    
    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;
    
    idx = find(ts1>baseAdjust,1);
    [~,baseSt] = min(abs(ts1 - baseline(1)));
    [~,baseEn] = min(abs(ts1 - baseline(2)));
    [~,ampSt] = min(abs(ts1 - amp_window(1)));
    [~,ampEn] = min(abs(ts1 - amp_window(2)));
    [~,aucSt] = min(abs(ts1 - auc_window(1)));
    [~,aucEn] = min(abs(ts1 - auc_window(2)));

    cueArray = session_identifiers(1:2:end-1,:);
    leverArray = session_identifiers(2:2:end,:);
    % if cueArray(1,1) > leverArray(1,1)
    %     leverArray = leverArray(2:end,:);
    % end
    cuelever_cRew = [];
    cuelever_cNoRew = [];
    cuelever_iRew = [];
    cuelever_iNoRew = [];
    levercue_cRew = [];
    levercue_cNoRew = [];
    levercue_iRew = [];
    levercue_iNoRew = [];
    cueTT_dFF = zeros(height(cueArray),epocArrayLen);
    cueTT_z = zeros(height(cueArray),epocArrayLen);
    cueTT_raw = zeros(height(cueArray),epocArrayLen);

    cuelever_Correct = [];
    cuelever_Incorrect = [];
    cue_amplitude = [];
    % some files end with a cue (0) instead of lever press %
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(3:end-1,:);
        cueArray = cueArray(1:end-1,1);
    end

    for m = removeFirstTrial:height(cueArray)
        cueTTStart = cueArray(m,1) - baseWindow;
        if cueTTStart < t
            continue
        end
        cueTTEnd = cueTTStart + timeWindow + baseWindow;
        [~,cTTSt] = min(abs(session_time - cueTTStart));
        [~,cTTEn] = min(abs(session_time - cueTTEnd));
        cueTTSigRaw = SIGNAL_raw(1,cTTSt:cTTEn);
        if length(cueTTSigRaw) < epocArrayLen
            mn = mean(cueTTSigRaw(1,end-10:end));
            cueTTSigRaw(1,end:epocArrayLen) = mn;
        elseif length(cueTTSigRaw) > epocArrayLen
            op = length(cueTTSigRaw);
            arrayDif = op - epocArrayLen;
            cueTTSigRaw = cueTTSigRaw(1,1:end-arrayDif);
        end
        cueTT_raw(m,:) = cueTTSigRaw;
                
    end

    %% loop for cue BEFORE trial %%
    for n = 1:height(cueTT_raw)
        % dF/F
        meanBase = mean(cueTT_raw(n,baseSt:baseEn));
        stdBase = std(cueTT_raw(n,baseSt:baseEn));
        
        cueTT_dFF(n,1:epocArrayLen) = cueTT_raw(n, 1:epocArrayLen) - meanBase;
        cueTT_dFF(n,1:epocArrayLen) = 100*(cueTT_dFF(n,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(cueTT_dFF(n,baseSt:baseEn));
        stdBase_dFF = std(cueTT_dFF(n,baseSt:baseEn));
        cueTT_z(n,1:epocArrayLen) = (cueTT_dFF(n,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        if adjustBase == 1
            if cueTT_z(n,idx) < 0
                val = cueTT_z(n,idx);
                diff = 0 - val;
                cueTT_z(n,1:epocArrayLen) = cueTT_z(n,1:epocArrayLen) + abs(diff);
            elseif cueTT_z(n,idx) > 0
                val = cueTT_z(n,idx);
                diff = 0 - val;
                cueTT_z(n,1:epocArrayLen) = cueTT_z(n,1:epocArrayLen) - abs(diff);
            end
        else
            disp("")
        end
        if removeOutlierTrials == 1
            tempAMP = max(cueTT_z(n,ampSt:ampEn));
            if isnan(tempAMP)
                continue
            elseif tempAMP < 1 || tempAMP > 20
                omitted(end+1).file = name;
                omitted(end).signal = cueTT_z(n,:);
                omitted(end).trial = n;
                omitted(end).timestamp = cueArray(n,1);
                cueTT_z(n,1:epocArrayLen) = nan;
            end
        end
        
        % sorts the cue psth by the next trial type %
        if leverArray(n,2) == 1
            cuelever_cRew = [cuelever_cRew;cueTT_z(n,:)];
            cuelever_Correct = [cuelever_Correct;cueTT_z(n,:)];
            % cue_amplitude(end+1,1) = max(cuelever_Correct(end,ampSt:ampEn));
        elseif leverArray(n,2) == 2
            cuelever_cNoRew = [cuelever_cNoRew;cueTT_z(n,:)];
            cuelever_Correct = [cuelever_Correct;cueTT_z(n,:)];
            % cue_amplitude(end+1,1) = max(cuelever_Correct(end,ampSt:ampEn));
        elseif leverArray(n,2) == 3
            cuelever_iRew = [cuelever_iRew;cueTT_z(n,:)];
            cuelever_Incorrect = [cuelever_Incorrect;cueTT_z(n,:)];
            % cue_amplitude(end+1,2) = max(cuelever_Incorrect(end,ampSt:ampEn));
        elseif leverArray(n,2) == 4
            cuelever_iNoRew = [cuelever_iNoRew;cueTT_z(n,:)];
            cuelever_Incorrect = [cuelever_Incorrect;cueTT_z(n,:)];
            % cue_amplitude(end+1,2) = max(cuelever_Incorrect(end,ampSt:ampEn));
        end
    end

    if isempty(cuelever_cRew)
        cuelever_cRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_cRewAMP(i,1) = nan;
        cuelever_cRewAUC(i,1) = nan;
    else
        cuelever_cRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_cRew,'omitnan');
        cuelever_cRewAMP(i,1) = calculateAMP(cuelever_cRewSTREAMSz(i,ampSt:ampEn));
        % % Calculate AUC above x=0 %
        % positive_indices = cuelever_cRewSTREAMSz(i,:) > 0;
        % y_pos = cuelever_cRewSTREAMSz(i,positive_indices);
        % x_pos = ts1(1,positive_indices);
        cuelever_cRewAUC(i,1) = calculateAUC(cuelever_cRewSTREAMSz(i,aucSt:aucEn),ts1(1, aucSt:aucEn));
    end
    if isempty(cuelever_cNoRew)
        cuelever_cNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_cNoRewAMP(i,1) = nan;
        cuelever_cNoRewAUC(i,1) = nan;
    else
        cuelever_cNoRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_cNoRew,'omitnan');
        cuelever_cNoRewAMP(i,1) = calculateAMP(cuelever_cNoRewSTREAMSz(i,ampSt:ampEn));
        % Calculate AUC above x=0 %
        % positive_indices = cuelever_cNoRewSTREAMSz(i,:) > 0;
        % y_pos = cuelever_cNoRewSTREAMSz(i,positive_indices);
        % x_pos = ts1(1,positive_indices);
        cuelever_cNoRewAUC(i,1) = calculateAUC(cuelever_cNoRewSTREAMSz(i,aucSt:aucEn),ts1(1, aucSt:aucEn));;
    end
    if isempty(cuelever_iRew)
        cuelever_iRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_iRewAMP(i,1) = nan;
        cuelever_iRewAUC(i,1) = nan;
    else
        cuelever_iRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_iRew,'omitnan');
        cuelever_iRewAMP(i,1) = calculateAMP(cuelever_iRewSTREAMSz(i,ampSt:ampEn));
        % % Calculate AUC above x=0 %
        % positive_indices = cuelever_iRewSTREAMSz(i,:) > 0;
        % y_pos = cuelever_iRewSTREAMSz(i,positive_indices);
        % x_pos = ts1(1,positive_indices);
        cuelever_iRewAUC(i,1) = calculateAUC(cuelever_iRewSTREAMSz(i,aucSt:aucEn),ts1(1, aucSt:aucEn));
    end
    if isempty(cuelever_iNoRew)
        cuelever_iNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_iNoRewAMP(i,1) = nan;
        cuelever_iNoRewAUC(i,1) = nan;
    else
        cuelever_iNoRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_iNoRew,'omitnan');
        cuelever_iNoRewAMP(i,1) = calculateAMP(cuelever_iNoRewSTREAMSz(i,ampSt:ampEn));
        % % Calculate AUC above x=0 %
        % positive_indices = cuelever_iNoRewSTREAMSz(i,:) > 0;
        % y_pos = cuelever_iNoRewSTREAMSz(i,positive_indices);
        % x_pos = ts1(1,positive_indices);
        cuelever_iNoRewAUC(i,1) = calculateAUC(cuelever_iNoRewSTREAMSz(i,aucSt:aucEn),ts1(1, aucSt:aucEn));;
    end
    if isempty(cuelever_Correct)
        cuelever_correctSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_correctAMP(i,1) = nan;
        cuelever_correctAUC(i,1) = nan;
    else
        for x = 1:height(cuelever_Correct)
            cuelever_Correct(sum(cuelever_Correct(x,:)) == 0) = [];
        end
        if height(cuelever_Correct) > 1   
            cuelever_correctSTREAMSz(i,1:epocArrayLen) = mean(cuelever_Correct,'omitnan');
        elseif height(cuelever_Correct) == 1
            cuelever_correctSTREAMSz(i,1:epocArrayLen) = cuelever_Correct;
        else
            disp('Error calculating mean for cuelever_Correct')
        end
        % if cuelever_correctSTREAMSz(i,idx) < 0
        %     val = cuelever_correctSTREAMSz(i,idx);
        %     diff = 0 - val;
        %     cuelever_correctSTREAMSz(i,1:epocArrayLen) = cuelever_correctSTREAMSz(i,1:epocArrayLen) + abs(diff);
        % elseif cuelever_correctSTREAMSz(i,idx) > 0
        %     val = cuelever_correctSTREAMSz(i,idx);
        %     diff = 0 - val;
        %     cuelever_correctSTREAMSz(i,1:epocArrayLen) = cuelever_correctSTREAMSz(i,1:epocArrayLen) - abs(diff);
        % end
        cuelever_correctAMP(i,1) = calculateAMP(cuelever_correctSTREAMSz(i,ampSt:ampEn));
        cuelever_correctAUC(i,1) = calculateAUC(cuelever_correctSTREAMSz(i,aucSt:aucEn),ts1(1,aucSt:aucEn));
    end
    if isempty(cuelever_Incorrect)
        cuelever_incorrectSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_incorrectAMP(i,1) = nan;
        cuelever_incorrectAUC(i,1) = nan;
    else
        for x = 1:height(cuelever_Incorrect)
            cuelever_Incorrect(sum(cuelever_Incorrect(x,:)) == 0) = [];
        end
        if height(cuelever_Incorrect) > 1   
            cuelever_incorrectSTREAMSz(i,1:epocArrayLen) = mean(cuelever_Incorrect,'omitnan');
        elseif height(cuelever_Incorrect) == 1
            cuelever_incorrectSTREAMSz(i,1:epocArrayLen) = cuelever_Incorrect;
        else
            disp('Error calculating mean for cuelever_Incorrect')
        end
        % if cuelever_incorrectSTREAMSz(i,idx) < 0
        %     val = cuelever_incorrectSTREAMSz(i,idx);
        %     diff = 0 - val;
        %     cuelever_incorrectSTREAMSz(i,1:epocArrayLen) = cuelever_incorrectSTREAMSz(i,1:epocArrayLen) + abs(diff);
        % elseif cuelever_incorrectSTREAMSz(i,idx) > 0
        %     val = cuelever_incorrectSTREAMSz(i,idx);
        %     diff = 0 - val;
        %     cuelever_incorrectSTREAMSz(i,1:epocArrayLen) = cuelever_incorrectSTREAMSz(i,1:epocArrayLen) - abs(diff);
        % end
        cuelever_incorrectAMP(i,1) = calculateAMP(cuelever_incorrectSTREAMSz(i,ampSt:ampEn));
        cuelever_incorrectAUC(i,1) = calculateAUC(cuelever_incorrectSTREAMSz(i,aucSt:aucEn),ts1(1,aucSt:aucEn));
    end

    %% loop for cue AFTER trial %%
    for n = 1:height(cueTT_raw)-1
        % dF/F
        meanBase = mean(cueTT_raw(n+1,baseSt:baseEn));
        stdBase = std(cueTT_raw(n+1,baseSt:baseEn));
        
        cueTT_dFF(n,1:epocArrayLen) = cueTT_raw(n, 1:epocArrayLen) - meanBase;
        cueTT_dFF(n,1:epocArrayLen) = 100*(cueTT_dFF(n,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(cueTT_dFF(n,baseSt:baseEn));
        stdBase_dFF = std(cueTT_dFF(n,baseSt:baseEn));
        cueTT_z(n,1:epocArrayLen) = (cueTT_dFF(n,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;

        % sorts the cue psth by the preceding trial type %
        if leverArray(n,2) == 1
            levercue_cRew = [levercue_cRew;cueTT_z(n+1,:)];
        elseif leverArray(n,2) == 2
            levercue_cNoRew = [levercue_cNoRew;cueTT_z(n+1,:)];
        elseif leverArray(n,2) == 3
            levercue_iRew = [levercue_iRew;cueTT_z(n+1,:)];
        elseif leverArray(n,2) == 4
            levercue_iNoRew = [levercue_iNoRew;cueTT_z(n+1,:)];
        end
    end

    if isempty(levercue_cRew)
        levercue_cRewSTREAMSz(i,1:epocArrayLen) = nan;
        levercue_cRewAMP(i,1) = nan;
        levercue_cRewAUC(i,1) = nan;
    else
        levercue_cRewSTREAMSz(i,1:epocArrayLen) = mean(levercue_cRew);
        levercue_cRewAMP(i,1) = max(levercue_cRewSTREAMSz(i,ampSt:ampEn));
        levercue_cRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),levercue_cRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(levercue_cNoRew)
        levercue_cNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        levercue_cNoRewAMP(i,1) = nan;
        levercue_cNoRewAUC(i,1) = nan;
    else
        levercue_cNoRewSTREAMSz(i,1:epocArrayLen) = mean(levercue_cNoRew);
        levercue_cNoRewAMP(i,1) = max(levercue_cNoRewSTREAMSz(i,ampSt:ampEn));
        levercue_cNoRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),levercue_cNoRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(levercue_iRew)
        levercue_iRewSTREAMSz(i,1:epocArrayLen) = nan;
        levercue_iRewAMP(i,1) = nan;
        levercue_iRewAUC(i,1) = nan;
    else
        levercue_iRewSTREAMSz(i,1:epocArrayLen) = mean(levercue_iRew);
        levercue_iRewAMP(i,1) = max(levercue_iRewSTREAMSz(i,ampSt:ampEn));
        levercue_iRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),levercue_iRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(levercue_iNoRew)
        levercue_iNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        levercue_iNoRewAMP(i,1) = nan;
        levercue_iNoRewAUC(i,1) = nan;
    else
        levercue_iNoRewSTREAMSz(i,1:epocArrayLen) = mean(levercue_iNoRew);
        levercue_iNoRewAMP(i,1) = max(levercue_iNoRewSTREAMSz(i,ampSt:ampEn));
        levercue_iNoRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),levercue_iNoRewSTREAMSz(i,aucSt:aucEn));
    end
end
IDs = cell2table(IDs,'VariableNames',{'ID'});
treatList = cell2table(treatList,'VariableNames',{'Treatment'});
phase = cell2table(prl_phase, 'VariableNames',{'Phase'});
%% cue before lever table %%
cuelever_cRewTable = array2table(cuelever_cRewSTREAMSz);
cuelever_cRewTable = horzcat(IDs,treatList,phase,cuelever_cRewTable);
cuelever_cNoRewTable = array2table(cuelever_cNoRewSTREAMSz);
cuelever_cNoRewTable = horzcat(IDs,treatList,phase,cuelever_cNoRewTable);
cuelever_iRewTable = array2table(cuelever_iRewSTREAMSz);
cuelever_iRewTable = horzcat(IDs,treatList,phase,cuelever_iRewTable);
cuelever_iNoRewTable = array2table(cuelever_iNoRewSTREAMSz);
cuelever_iNoRewTable = horzcat(IDs,treatList,phase,cuelever_iNoRewTable);

cuelever_correctTable = array2table(cuelever_correctSTREAMSz);
cuelever_correctTable = horzcat(IDs,treatList,phase,cuelever_correctTable);
cuelever_incorrectTable = array2table(cuelever_incorrectSTREAMSz);
cuelever_incorrectTable = horzcat(IDs,treatList,phase,cuelever_incorrectTable);



cueleverTTAMP = [cuelever_cRewAMP cuelever_cNoRewAMP cuelever_iRewAMP cuelever_iNoRewAMP];
cueleverTTAMP = array2table(cueleverTTAMP,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
cueleverTTAMP = horzcat(IDs,treatList,phase,cueleverTTAMP);
cueleverTTAUC = [cuelever_cRewAUC cuelever_cNoRewAUC cuelever_iRewAUC cuelever_iNoRewAUC];
cueleverTTAUC = array2table(cueleverTTAUC,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
cueleverTTAUC = horzcat(IDs,treatList,phase,cueleverTTAUC);

cueleverAMP = [cuelever_correctAMP cuelever_incorrectAMP];
cueleverAMP = array2table(cueleverAMP,'VariableNames',{'cue_correct','cue_incorrect'});
cueleverAMP = horzcat(IDs,treatList,phase,cueleverAMP);
cueleverAUC = [cuelever_correctAUC cuelever_incorrectAUC];
cueleverAUC = array2table(cueleverAUC,'VariableNames',{'cue_correct','cue_incorrect'});
cueleverAUC = horzcat(IDs,treatList,phase,cueleverAUC);

%% cue after lever table %%
levercue_cRewTable = array2table(levercue_cRewSTREAMSz);
levercue_cRewTable = horzcat(IDs,levercue_cRewTable);
levercue_cNoRewTable = array2table(levercue_cNoRewSTREAMSz);
levercue_cNoRewTable = horzcat(IDs,levercue_cNoRewTable);
levercue_iRewTable = array2table(levercue_iRewSTREAMSz);
levercue_iRewTable = horzcat(IDs,levercue_iRewTable);
levercue_iNoRewTable = array2table(levercue_iNoRewSTREAMSz);
levercue_iNoRewTable = horzcat(IDs,levercue_iNoRewTable);

levercueTTAMP = [levercue_cRewAMP levercue_cNoRewAMP levercue_iRewAMP levercue_iNoRewAMP];
levercueTTAMP = array2table(levercueTTAMP,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
levercueTTAMP = horzcat(IDs,levercueTTAMP);
levercueTTAUC = [levercue_cRewAUC levercue_cNoRewAUC levercue_iRewAUC levercue_iNoRewAUC];
levercueTTAUC = array2table(levercueTTAUC,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
levercueTTAUC = horzcat(IDs,levercueTTAUC);

cuelever_correctTable = sortrows(cuelever_correctTable,{'Phase','Treatment'},{'ascend','descend'});
cuelever_incorrectTable = sortrows(cuelever_incorrectTable,{'Phase','Treatment'},{'ascend','descend'});

cueleverAMP = sortrows(cueleverAMP,{'Phase','Treatment'},{'ascend','descend'});
cueleverAUC = sortrows(cueleverAUC,{'Phase','Treatment'},{'ascend','descend'});

cuelever_correctTable = rows2vars(cuelever_correctTable);
cuelever_incorrectTable = rows2vars(cuelever_incorrectTable);

prl_cue.cuebeforecorrect = cuelever_correctTable;
prl_cue.cuebeforeincorrect = cuelever_incorrectTable;
prl_cue.cueamp = cueleverAMP;
prl_cue.cueauc = cueleverAUC;
%% cue before lever fig %%
f1 = figure;
subplot(4,1,1)
plot(ts1,mean(cuelever_cRewSTREAMSz,1,'omitnan'),'b')
hold on
f1std = std(cuelever_cRewSTREAMSz,'omitnan')/sqrt(height(cuelever_cRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_cRewSTREAMSz,1,'omitnan') + f1std, fliplr(mean(cuelever_cRewSTREAMSz,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('C+ Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,2)
plot(ts1,mean(cuelever_cNoRewSTREAMSz,1,'omitnan'),'b')
hold on
f1std = std(cuelever_cNoRewSTREAMSz,'omitnan')/sqrt(height(cuelever_cNoRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_cNoRewSTREAMSz,1,'omitnan') + f1std, fliplr(mean(cuelever_cNoRewSTREAMSz,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('C- Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,3)
plot(ts1,mean(cuelever_iRewSTREAMSz,1,'omitnan'),'b')
hold on
f1std = std(cuelever_iRewSTREAMSz,'omitnan')/sqrt(height(cuelever_iRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_iRewSTREAMSz,1,'omitnan') + f1std, fliplr(mean(cuelever_iRewSTREAMSz,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('I+ Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,4)
plot(ts1,mean(cuelever_iNoRewSTREAMSz,1,'omitnan'),'b')
hold on
f1std = std(cuelever_iNoRewSTREAMSz,'omitnan')/sqrt(height(cuelever_iNoRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_iNoRewSTREAMSz,1,'omitnan') + f1std, fliplr(mean(cuelever_iNoRewSTREAMSz,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('I- Cue')
ylabel('GrabDA dF/F (z-Score)')
sgtitle('Cue Before Trial')

%% cue after lever fig %%
f2 = figure;
subplot(4,1,1)
plot(ts1,mean(levercue_cRewSTREAMSz,1,'omitnan'),'b')
hold on
f2std = std(levercue_cRewSTREAMSz,'omitnan')/sqrt(height(levercue_cRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(levercue_cRewSTREAMSz,1,'omitnan') + f2std, fliplr(mean(levercue_cRewSTREAMSz,1,'omitnan') - f2std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('C+ Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,2)
plot(ts1,mean(levercue_cNoRewSTREAMSz,1,'omitnan'),'b')
hold on
f2std = std(levercue_cNoRewSTREAMSz,'omitnan')/sqrt(height(levercue_cNoRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(levercue_cNoRewSTREAMSz,1,'omitnan') + f2std, fliplr(mean(levercue_cNoRewSTREAMSz,1,'omitnan') - f2std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('C- Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,3)
plot(ts1,mean(levercue_iRewSTREAMSz,1,'omitnan'),'b')
hold on
f2std = std(levercue_iRewSTREAMSz,'omitnan')/sqrt(height(levercue_iRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(levercue_iRewSTREAMSz,1,'omitnan') + f2std, fliplr(mean(levercue_iRewSTREAMSz,1,'omitnan') - f2std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('I+ Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(4,1,4)
plot(ts1,mean(levercue_iNoRewSTREAMSz,1,'omitnan'),'b')
hold on
f2std = std(levercue_iNoRewSTREAMSz,'omitnan')/sqrt(height(levercue_iNoRewSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(levercue_iNoRewSTREAMSz,1,'omitnan') + f2std, fliplr(mean(levercue_iNoRewSTREAMSz,1,'omitnan') - f2std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('I- Cue')
ylabel('GrabDA dF/F (z-Score)')
sgtitle('Cue After Trial')

%% cue before lever correct vs incorrect %%
f3 = figure;
subplot(2,1,1)
plot(ts1,mean(cuelever_correctSTREAMSz,1,'omitnan'),'b')
hold on
f3std = std(cuelever_correctSTREAMSz,'omitnan')/sqrt(height(cuelever_correctSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_correctSTREAMSz,1,'omitnan') + f3std, fliplr(mean(cuelever_correctSTREAMSz,1,'omitnan') - f3std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Cue -> Correct')
ylabel('GrabDA dF/F (z-Score)')

subplot(2,1,2)
plot(ts1,mean(cuelever_incorrectSTREAMSz,1,'omitnan'),'b')
hold on
f3std = std(cuelever_incorrectSTREAMSz,'omitnan')/sqrt(height(cuelever_incorrectSTREAMSz));
XX = [ts1, fliplr(ts1)];
YY = [mean(cuelever_incorrectSTREAMSz,1,'omitnan') + f3std, fliplr(mean(cuelever_incorrectSTREAMSz,1,'omitnan') - f3std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Cue -> Incorrect')
ylabel('GrabDA dF/F (z-Score)')

if plotOmitted == 1
    % Plot omitted signals on the same graph
    f4 = figure;
    hold on
    for ii = 1:length(omitted)
        plot(ts1,omitted(ii).signal);
        
    end
    hold off
    xlabel('Time');
    ylabel('zScore');
    title('Omitted Signals');
    % Create a cell array of filenames for the legend
    legendLabels = {};
    
    % Convert filename values to strings and add them to legendLabels
    for ii = 1:length(omitted)
        legendLabels{end+1} = sprintf('Trial %d - %s', omitted(ii).trial, char(omitted(ii).file));
    end
    
    % Display legend if there are valid labels
    if ~isempty(legendLabels)
        legend(legendLabels);
    end
    f5 = figure;
    plot(session_time,SIGNAL_z)
    hold on
    for ii = 1:length(omitted)
        xline(omitted(ii).timestamp, 'LineWidth',0.5,'Color','red')
    end
    hold off
end

f6 = figure;
bar1 = mean(cuelever_correctAMP,'omitnan');
bar2 = mean(cuelever_incorrectAMP,'omitnan');
barData = [bar1,bar2];
bar(barData);
xticks(1:2);
xticklabels({'Correct','Incorrect'});

toc

NERD_STATS(toc,numFiles);
clearvars -except prl_cue