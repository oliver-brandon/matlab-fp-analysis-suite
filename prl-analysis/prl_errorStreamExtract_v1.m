clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0 2]; % time window to grab amplitude from
auc_window = [0 timeWindow];
t = 5; % seconds to clip from start of streams
N = 10; % Downsample N times
adjust = -1; % Locks the streams vertically to 'adjust' seconds
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
errorType = 4; % 1 = winStay, 2 = winShift, 3 = loseStay, 4 = loseShift
lever = 1; % errors based on correct lever = 1, incorrect lever = 2
errorPos = 2; % 1 = WIN/LOSE, 2 = STAY/SHIFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(...
    '/Users/brandon/My Drive/prl/PRL_GRABDA/goodSignals','Choose the .mat files you want to analyze.'...
    ); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
tic

if lever == 1
    levStr = 'Lever: Correct';
else
    levStr = 'Lever: Incorrect';
end
    
if errorType == 1 && errorPos == 1
    errStr = 'Trial Outcome: WIN-stay';
elseif errorType == 1 && errorPos == 2
    errStr = 'Trial Outcome: win-STAY';
elseif errorType == 2 && errorPos == 1
    errStr = 'Trial Outcome: WIN-shift';
elseif errorType == 2 && errorPos == 2
    errStr = 'Trial Outcome: win-SHIFT';
elseif errorType == 3 && errorPos == 1
    errStr = 'Trial Outcome: LOSE-stay';
elseif errorType == 3 && errorPos == 2
    errStr = 'Trial Outcome: lose-STAY';
elseif errorType == 4 && errorPos == 1
    errStr = 'Trial Outcome: LOSE-shift';
elseif errorType == 4 && errorPos == 2
    errStr = 'Trial Outcome: lose-SHIFT';
end
TITLE_Lev = strcat('Epoch: Lever Press',{' -- '},levStr,{' -- '},errStr);
TITLE_Cue = strcat('Epoch: Cue Light',{' -- '},levStr,{' -- '},errStr);
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
IDs = {};
phaseList = {};
treatList = {};
master_LevStreams_z = [];
master_LevAmp_z = [];
master_LevAuc_z = [];
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    brokenID = strsplit(name,'_');
    IDs{i,1} = cellstr(strtrim(brokenID{1}));
    phaseList{i,1} = cellstr(strtrim(brokenID{2}));
    treatList{i,1} = cellstr(strtrim(brokenID{3}));
    fprintf('Analyzing file %d of %d\n',i,numFiles);
    load(filename)

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
        [data, errorProbLeverTS, errorProbCueTS] = errorProbExtract(data, session_identifiers, errorType, lever);
        
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
        [data, errorProbLeverTS, errorProbCueTS] = errorProbExtract(data, session_identifiers, errorType, lever);
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

    % [SIGNAL_raw] = correctSignal(ISOS_raw,SIGNAL_raw);
    
    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;
    % checks if the lever array is empty %
    if isempty(errorProbLeverTS)
        master_LevStreams_z(i,1:epocArrayLen) = NaN;
        master_LevAmp_z(i,1) = NaN;
        master_LevAuc_z(i,1) = NaN;
        master_CueStreams_z(i,1:epocArrayLen) = NaN;
        master_CueAmp_z(i,1) = NaN;
        master_CueAuc_z(i,1) = NaN;
        continue
    end
    %% extracts lever streams, amp, and auc %%
    idx = find(ts1>adjust,1);
    for ii = 1:height(errorProbLeverTS(:,errorPos))  
        windowStart = errorProbLeverTS(ii,errorPos)-baseWindow;
        windowEnd = windowStart+timeWindow+baseWindow;
        [~,windSt] = min(abs(session_time - windowStart));
        [~,windEn] = min(abs(session_time - windowEnd));
        epocSigRaw_lev = SIGNAL_raw(1,windSt:windEn);

        if length(epocSigRaw_lev) < epocArrayLen
            mn = mean(epocSigRaw_lev(1,end-10:end));
            epocSigRaw_lev(1,end:epocArrayLen) = mn;
        elseif length(epocSigRaw_lev) > epocArrayLen
            op = length(epocSigRaw_lev);
            arrayDif = op - epocArrayLen;
            epocSigRaw_lev = epocSigRaw_lev(1,1:end-arrayDif);
        end
        streams_raw_lev(ii,1:epocArrayLen) = epocSigRaw_lev;

    end

    [~,baseSt] = min(abs(ts1 - (baseline(1))));
    [~,baseEn] = min(abs(ts1 - (baseline(2))));
    [~,ampSt] = min(abs(ts1 - (amp_window(1))));
    [~,ampEn] = min(abs(ts1 - (amp_window(2))));
    [~,aucSt] = min(abs(ts1 - (auc_window(1))));
    [~,aucEn] = min(abs(ts1 - (auc_window(2))));

    positive_indices = [];
    for j = 1:height(streams_raw_lev)
        % dF/F
        meanBase = mean(streams_raw_lev(j,baseSt:baseEn));
        stdBase = std(streams_raw_lev(j,baseSt:baseEn));
        streams_dFF(j,1:epocArrayLen) = streams_raw_lev(j,1:epocArrayLen) - meanBase;
        streams_dFF(j,1:epocArrayLen) = 100*(streams_dFF(j,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
        stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
        streams_z_lev(j,1:epocArrayLen) = (streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % adjusts streams to baseline of zero at -0.5s %
        if streams_z_lev(j,idx) < 0
            val = streams_z_lev(j,idx);
            diff = 0 - val;
            streams_z_lev(j,1:epocArrayLen) = streams_z_lev(j,1:epocArrayLen) + abs(diff);
        elseif streams_z_lev(j,idx) > 0
            val = streams_z_lev(j,idx);
            diff = 0 - val;
            streams_z_lev(j,1:epocArrayLen) = streams_z_lev(j,1:epocArrayLen) - abs(diff);
        end

    end
    master_LevStreams_z(i,:) = mean(streams_z_lev,'omitnan');
    master_LevAmp_z(i,1) = calculateAMP(master_LevStreams_z(i,ampSt:ampEn),ts1(1,ampSt:ampEn));
    master_LevAuc_z(i,1) = calculateAUC(master_LevStreams_z(i,aucSt:aucEn),ts1(1,aucSt:aucEn));
    
    %% extracts cue streams, amp, and auc %%
    idx = find(ts1>adjust,1);
    for ii = 1:height(errorProbCueTS(:,errorPos))  
        windowStart = errorProbCueTS(ii,errorPos)-baseWindow;
        windowEnd = windowStart+timeWindow+baseWindow;
        [~,windSt] = min(abs(session_time - windowStart));
        [~,windEn] = min(abs(session_time - windowEnd));
        epocSigRaw_cue = SIGNAL_raw(1,windSt:windEn);
    
        if length(epocSigRaw_cue) < epocArrayLen
            mn = mean(epocSigRaw_cue(1,end-10:end));
            epocSigRaw_cue(1,end:epocArrayLen) = mn;
        elseif length(epocSigRaw_cue) > epocArrayLen
            op = length(epocSigRaw_cue);
            arrayDif = op - epocArrayLen;
            epocSigRaw_cue = epocSigRaw_cue(1,1:end-arrayDif);
        end
        streams_raw_cue(ii,1:epocArrayLen) = epocSigRaw_cue;
    
    end
    
    positive_indices = [];
    for j = 1:height(streams_raw_cue)
        % dF/F
        meanBase = mean(streams_raw_cue(j,baseSt:baseEn));
        stdBase = std(streams_raw_cue(j,baseSt:baseEn));
        streams_dFF(j,1:epocArrayLen) = streams_raw_cue(j,1:epocArrayLen) - meanBase;
        streams_dFF(j,1:epocArrayLen) = 100*(streams_dFF(j,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
        stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
        streams_z_cue(j,1:epocArrayLen) = (streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % adjusts streams to baseline of zero at -0.5s %
        if streams_z_cue(j,idx) < 0
            val = streams_z_cue(j,idx);
            diff = 0 - val;
            streams_z_cue(j,1:epocArrayLen) = streams_z_cue(j,1:epocArrayLen) + abs(diff);
        elseif streams_z_cue(j,idx) > 0
            val = streams_z_cue(j,idx);
            diff = 0 - val;
            streams_z_cue(j,1:epocArrayLen) = streams_z_cue(j,1:epocArrayLen) - abs(diff);
        end


    end

    master_CueStreams_z(i,:) = mean(streams_z_cue,'omitnan');
    master_CueAmp_z(i,1) = calculateAMP(master_CueStreams_z(i,ampSt:ampEn),ts1(1,ampSt:ampEn));
    master_CueAuc_z(i,1) = calculateAUC(master_CueStreams_z(i,aucSt:aucEn),ts1(1,aucSt:aucEn));
end
IDs = cell2table(IDs,'VariableNames',{'ID'});
treatList = cell2table(treatList,'VariableNames',{'Treatment'});
phaseList = cell2table(phaseList, 'VariableNames',{'Phase'});

errorProbStreams_lev = array2table(master_LevStreams_z);
errorProbStreams_lev = horzcat(IDs,treatList,phaseList,errorProbStreams_lev);
errorProbStreams_lev = sortrows(errorProbStreams_lev,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});
errorProbAmp_lev = array2table(master_LevAmp_z,'VariableNames',{'Amplitude'});
errorProbAmp_lev = horzcat(IDs,treatList,phaseList,errorProbAmp_lev);
errorProbAmp_lev = sortrows(errorProbAmp_lev,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});
errorProbAUC_lev = array2table(master_LevAuc_z,'VariableNames',{'AUC'});
errorProbAUC_lev = horzcat(IDs,treatList,phaseList,errorProbAUC_lev);
errorProbAUC_lev = sortrows(errorProbAUC_lev,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});

prl_error_analysis.acq1.lev_amplitude = errorProbAmp_lev(strcmp(errorProbAmp_lev.Phase,'Acq1'),:);
prl_error_analysis.acq2.lev_amplitude = errorProbAmp_lev(strcmp(errorProbAmp_lev.Phase,'Acq2'),:);
prl_error_analysis.rev1.lev_amplitude = errorProbAmp_lev(strcmp(errorProbAmp_lev.Phase,'Rev1'),:);
prl_error_analysis.rev2.lev_amplitude = errorProbAmp_lev(strcmp(errorProbAmp_lev.Phase,'Rev2'),:);
prl_error_analysis.rev3.lev_amplitude = errorProbAmp_lev(strcmp(errorProbAmp_lev.Phase,'Rev3'),:);
prl_error_analysis.acq1.lev_auc = errorProbAUC_lev(strcmp(errorProbAUC_lev.Phase,'Acq1'),:);
prl_error_analysis.acq2.lev_auc = errorProbAUC_lev(strcmp(errorProbAUC_lev.Phase,'Acq2'),:);
prl_error_analysis.rev1.lev_auc = errorProbAUC_lev(strcmp(errorProbAUC_lev.Phase,'Rev1'),:);
prl_error_analysis.rev2.lev_auc = errorProbAUC_lev(strcmp(errorProbAUC_lev.Phase,'Rev2'),:);
prl_error_analysis.rev3.lev_auc = errorProbAUC_lev(strcmp(errorProbAUC_lev.Phase,'Rev3'),:);
prl_error_analysis.acq1.lev_streams = errorProbStreams_lev(strcmp(errorProbStreams_lev.Phase,'Acq1'),:);
prl_error_analysis.acq2.lev_streams = errorProbStreams_lev(strcmp(errorProbStreams_lev.Phase,'Acq2'),:);
prl_error_analysis.rev1.lev_streams = errorProbStreams_lev(strcmp(errorProbStreams_lev.Phase,'Rev1'),:);
prl_error_analysis.rev2.lev_streams = errorProbStreams_lev(strcmp(errorProbStreams_lev.Phase,'Rev2'),:);
prl_error_analysis.rev3.lev_streams = errorProbStreams_lev(strcmp(errorProbStreams_lev.Phase,'Rev3'),:);

errorProbStreams_cue = array2table(master_CueStreams_z);
errorProbStreams_cue = horzcat(IDs,treatList,phaseList,errorProbStreams_cue);
errorProbStreams_cue = sortrows(errorProbStreams_cue,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});
errorProbAmp_cue = array2table(master_CueAmp_z,'VariableNames',{'Amplitude'});
errorProbAmp_cue = horzcat(IDs,treatList,phaseList,errorProbAmp_cue);
errorProbAmp_cue = sortrows(errorProbAmp_cue,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});
errorProbAUC_cue = array2table(master_CueAuc_z,'VariableNames',{'AUC'});
errorProbAUC_cue = horzcat(IDs,treatList,phaseList,errorProbAUC_cue);
errorProbAUC_cue = sortrows(errorProbAUC_cue,{'Phase','Treatment','ID'},{'ascend','descend','ascend'});

prl_error_analysis.acq1.cue_amplitude = errorProbAmp_cue(strcmp(errorProbAmp_cue.Phase,'Acq1'),:);
prl_error_analysis.acq2.cue_amplitude = errorProbAmp_cue(strcmp(errorProbAmp_cue.Phase,'Acq2'),:);
prl_error_analysis.rev1.cue_amplitude = errorProbAmp_cue(strcmp(errorProbAmp_cue.Phase,'Rev1'),:);
prl_error_analysis.rev2.cue_amplitude = errorProbAmp_cue(strcmp(errorProbAmp_cue.Phase,'Rev2'),:);
prl_error_analysis.rev3.cue_amplitude = errorProbAmp_cue(strcmp(errorProbAmp_cue.Phase,'Rev3'),:);
prl_error_analysis.acq1.cue_auc = errorProbAUC_cue(strcmp(errorProbAUC_cue.Phase,'Acq1'),:);
prl_error_analysis.acq2.cue_auc = errorProbAUC_cue(strcmp(errorProbAUC_cue.Phase,'Acq2'),:);
prl_error_analysis.rev1.cue_auc = errorProbAUC_cue(strcmp(errorProbAUC_cue.Phase,'Rev1'),:);
prl_error_analysis.rev2.cue_auc = errorProbAUC_cue(strcmp(errorProbAUC_cue.Phase,'Rev2'),:);
prl_error_analysis.rev3.cue_auc = errorProbAUC_cue(strcmp(errorProbAUC_cue.Phase,'Rev3'),:);
prl_error_analysis.acq1.cue_streams = errorProbStreams_cue(strcmp(errorProbStreams_cue.Phase,'Acq1'),:);
prl_error_analysis.acq2.cue_streams = errorProbStreams_cue(strcmp(errorProbStreams_cue.Phase,'Acq2'),:);
prl_error_analysis.rev1.cue_streams = errorProbStreams_cue(strcmp(errorProbStreams_cue.Phase,'Rev1'),:);
prl_error_analysis.rev2.cue_streams = errorProbStreams_cue(strcmp(errorProbStreams_cue.Phase,'Rev2'),:);
prl_error_analysis.rev3.cue_streams = errorProbStreams_cue(strcmp(errorProbStreams_cue.Phase,'Rev3'),:);

%% Filters Tables to Plot %%
% Acq2 Lever %
streamFilter_VEHLevAcq2 = table2array(prl_error_analysis.acq2.lev_streams(strcmp(prl_error_analysis.acq2.lev_streams.Treatment,'VEH'),4:end));
ampFilter_VEHLevAcq2 = table2array(prl_error_analysis.acq2.lev_amplitude(strcmp(prl_error_analysis.acq2.lev_amplitude.Treatment,'VEH'),4:end));
aucFilter_VEHLevAcq2 = table2array(prl_error_analysis.acq2.lev_auc(strcmp(prl_error_analysis.acq2.lev_auc.Treatment,'VEH'),4:end));

streamFilter_JZL8LevAcq2 = table2array(prl_error_analysis.acq2.lev_streams(strcmp(prl_error_analysis.acq2.lev_streams.Treatment,'JZL8'),4:end));
ampFilter_JZL8LevAcq2 = table2array(prl_error_analysis.acq2.lev_amplitude(strcmp(prl_error_analysis.acq2.lev_amplitude.Treatment,'JZL8'),4:end));
aucFilter_JZL8LevAcq2 = table2array(prl_error_analysis.acq2.lev_auc(strcmp(prl_error_analysis.acq2.lev_auc.Treatment,'JZL8'),4:end));

streamFilter_JZL18LevAcq2 = table2array(prl_error_analysis.acq2.lev_streams(strcmp(prl_error_analysis.acq2.lev_streams.Treatment,'JZL18'),4:end));
ampFilter_JZL18LevAcq2 = table2array(prl_error_analysis.acq2.lev_amplitude(strcmp(prl_error_analysis.acq2.lev_amplitude.Treatment,'JZL18'),4:end));
aucFilter_JZL18LevAcq2 = table2array(prl_error_analysis.acq2.lev_auc(strcmp(prl_error_analysis.acq2.lev_auc.Treatment,'JZL18'),4:end));

% Acq2 Cue %
streamFilter_VEHCueAcq2 = table2array(prl_error_analysis.acq2.cue_streams(strcmp(prl_error_analysis.acq2.cue_streams.Treatment,'VEH'),4:end));
ampFilter_VEHCueAcq2 = table2array(prl_error_analysis.acq2.cue_amplitude(strcmp(prl_error_analysis.acq2.cue_amplitude.Treatment,'VEH'),4:end));
aucFilter_VEHCueAcq2 = table2array(prl_error_analysis.acq2.cue_auc(strcmp(prl_error_analysis.acq2.cue_auc.Treatment,'VEH'),4:end));

streamFilter_JZL8CueAcq2 = table2array(prl_error_analysis.acq2.cue_streams(strcmp(prl_error_analysis.acq2.cue_streams.Treatment,'JZL8'),4:end));
ampFilter_JZL8CueAcq2 = table2array(prl_error_analysis.acq2.cue_amplitude(strcmp(prl_error_analysis.acq2.cue_amplitude.Treatment,'JZL8'),4:end));
aucFilter_JZL8CueAcq2 = table2array(prl_error_analysis.acq2.cue_auc(strcmp(prl_error_analysis.acq2.cue_auc.Treatment,'JZL8'),4:end));

streamFilter_JZL18CueAcq2 = table2array(prl_error_analysis.acq2.cue_streams(strcmp(prl_error_analysis.acq2.cue_streams.Treatment,'JZL18'),4:end));
ampFilter_JZL18CueAcq2 = table2array(prl_error_analysis.acq2.cue_amplitude(strcmp(prl_error_analysis.acq2.cue_amplitude.Treatment,'JZL18'),4:end));
aucFilter_JZL18CueAcq2 = table2array(prl_error_analysis.acq2.cue_auc(strcmp(prl_error_analysis.acq2.cue_auc.Treatment,'JZL18'),4:end));

% Rev1 Lever %
streamFilter_VEHLevRev1 = table2array(prl_error_analysis.rev1.lev_streams(strcmp(prl_error_analysis.rev1.lev_streams.Treatment,'VEH'),4:end));
ampFilter_VEHLevRev1 = table2array(prl_error_analysis.rev1.lev_amplitude(strcmp(prl_error_analysis.rev1.lev_amplitude.Treatment,'VEH'),4:end));
aucFilter_VEHLevRev1 = table2array(prl_error_analysis.rev1.lev_auc(strcmp(prl_error_analysis.rev1.lev_auc.Treatment,'VEH'),4:end));

streamFilter_JZL8LevRev1 = table2array(prl_error_analysis.rev1.lev_streams(strcmp(prl_error_analysis.rev1.lev_streams.Treatment,'JZL8'),4:end));
ampFilter_JZL8LevRev1 = table2array(prl_error_analysis.rev1.lev_amplitude(strcmp(prl_error_analysis.rev1.lev_amplitude.Treatment,'JZL8'),4:end));
aucFilter_JZL8LevRev1 = table2array(prl_error_analysis.rev1.lev_auc(strcmp(prl_error_analysis.rev1.lev_auc.Treatment,'JZL8'),4:end));

streamFilter_JZL18LevRev1 = table2array(prl_error_analysis.rev1.lev_streams(strcmp(prl_error_analysis.rev1.lev_streams.Treatment,'JZL18'),4:end));
ampFilter_JZL18LevRev1 = table2array(prl_error_analysis.rev1.lev_amplitude(strcmp(prl_error_analysis.rev1.lev_amplitude.Treatment,'JZL18'),4:end));
aucFilter_JZL18LevRev1 = table2array(prl_error_analysis.rev1.lev_auc(strcmp(prl_error_analysis.rev1.lev_auc.Treatment,'JZL18'),4:end));

% Rev1 Cue %
streamFilter_VEHCueRev1 = table2array(prl_error_analysis.rev1.cue_streams(strcmp(prl_error_analysis.rev1.cue_streams.Treatment,'VEH'),4:end));
ampFilter_VEHCueRev1 = table2array(prl_error_analysis.rev1.cue_amplitude(strcmp(prl_error_analysis.rev1.cue_amplitude.Treatment,'VEH'),4:end));
aucFilter_VEHCueRev1 = table2array(prl_error_analysis.rev1.cue_auc(strcmp(prl_error_analysis.rev1.cue_auc.Treatment,'VEH'),4:end));

streamFilter_JZL8CueRev1 = table2array(prl_error_analysis.rev1.cue_streams(strcmp(prl_error_analysis.rev1.cue_streams.Treatment,'JZL8'),4:end));
ampFilter_JZL8CueRev1 = table2array(prl_error_analysis.rev1.cue_amplitude(strcmp(prl_error_analysis.rev1.cue_amplitude.Treatment,'JZL8'),4:end));
aucFilter_JZL8CueRev1 = table2array(prl_error_analysis.rev1.cue_auc(strcmp(prl_error_analysis.rev1.cue_auc.Treatment,'JZL8'),4:end));

streamFilter_JZL18CueRev1 = table2array(prl_error_analysis.rev1.cue_streams(strcmp(prl_error_analysis.rev1.cue_streams.Treatment,'JZL18'),4:end));
ampFilter_JZL18CueRev1 = table2array(prl_error_analysis.rev1.cue_amplitude(strcmp(prl_error_analysis.rev1.cue_amplitude.Treatment,'JZL18'),4:end));
aucFilter_JZL18CueRev1 = table2array(prl_error_analysis.rev1.cue_auc(strcmp(prl_error_analysis.rev1.cue_auc.Treatment,'JZL18'),4:end));



%% Acq2 vs Rev1 %%
% Plots Lever %
f1 = figure;
set(f1,'Position',[1660,200,1600,1000]);
subplot(3,2,1)
plot(ts1,mean(streamFilter_VEHLevAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_VEHLevAcq2,'omitnan')/sqrt(height(streamFilter_VEHLevAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_VEHLevAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_VEHLevAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: Vehicle'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,3)
plot(ts1,mean(streamFilter_JZL8LevAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL8LevAcq2,'omitnan')/sqrt(height(streamFilter_JZL8LevAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL8LevAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL8LevAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: JZL8'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,5)
plot(ts1,mean(streamFilter_JZL18LevAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL18LevAcq2,'omitnan')/sqrt(height(streamFilter_JZL18LevAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL18LevAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL18LevAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: JZL18'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,2)
plot(ts1,mean(streamFilter_VEHLevRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_VEHLevRev1,'omitnan')/sqrt(height(streamFilter_VEHLevRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_VEHLevRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_VEHLevRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: Vehicle'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,4)
plot(ts1,mean(streamFilter_JZL8LevRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL8LevRev1,'omitnan')/sqrt(height(streamFilter_JZL8LevRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL8LevRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL8LevRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: JZL8'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,6)
plot(ts1,mean(streamFilter_JZL18LevRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL18LevRev1,'omitnan')/sqrt(height(streamFilter_JZL18LevRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL18LevRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL18LevRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: JZL18'))
ylabel('GrabDA dF/F (z-Score)')

[ax1,h1] = suplabel(TITLE_Lev,'t');
set(h1,'FontSize',20)
% Plots Cue %
f2 = figure;
set(f2,'Position',[1660,200,1600,1000]);
subplot(3,2,1)
plot(ts1,mean(streamFilter_VEHCueAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_VEHCueAcq2,'omitnan')/sqrt(height(streamFilter_VEHCueAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_VEHCueAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_VEHCueAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: Vehicle'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,3)
plot(ts1,mean(streamFilter_JZL8CueAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL8CueAcq2,'omitnan')/sqrt(height(streamFilter_JZL8CueAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL8CueAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL8CueAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: JZL8'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,5)
plot(ts1,mean(streamFilter_JZL18CueAcq2,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL18CueAcq2,'omitnan')/sqrt(height(streamFilter_JZL18CueAcq2));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL18CueAcq2,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL18CueAcq2,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: Last Acq',{' -- '},'Treatment: JZL18'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,2)
plot(ts1,mean(streamFilter_VEHCueRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_VEHCueRev1,'omitnan')/sqrt(height(streamFilter_VEHCueRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_VEHCueRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_VEHCueRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: Vehicle'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,4)
plot(ts1,mean(streamFilter_JZL8CueRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL8CueRev1,'omitnan')/sqrt(height(streamFilter_JZL8CueRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL8CueRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL8CueRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: JZL8'))
ylabel('GrabDA dF/F (z-Score)')

subplot(3,2,6)
plot(ts1,mean(streamFilter_JZL18CueRev1,1,'omitnan'),'b')
hold on
f1std = std(streamFilter_JZL18CueRev1,'omitnan')/sqrt(height(streamFilter_JZL18CueRev1));
XX = [ts1, fliplr(ts1)];
YY = [mean(streamFilter_JZL18CueRev1,1,'omitnan') + f1std, fliplr(mean(streamFilter_JZL18CueRev1,1,'omitnan') - f1std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
yline(0,'LineWidth',2,'Color','black')
title(strcat('Phase: First Rev',{' -- '},'Treatment: JZL18'))
ylabel('GrabDA dF/F (z-Score)')
[ax2,h2] = suplabel(TITLE_Cue,'t');
set(h2,'FontSize',20)

% means for bar data %
AUC_veh_cue_acq2 = mean(aucFilter_VEHCueAcq2,1,'omitnan');
AUC_jzl8_cue_acq2 = mean(aucFilter_JZL8CueAcq2,1,'omitnan');
AUC_jzl18_cue_acq2 = mean(aucFilter_JZL18CueAcq2,1,'omitnan');

AUC_veh_cue_rev1 = mean(aucFilter_VEHCueRev1,1,'omitnan');
AUC_jzl8_cue_rev1 = mean(aucFilter_JZL8CueRev1,1,'omitnan');
AUC_jzl18_cue_rev1 = mean(aucFilter_JZL18CueRev1,1,'omitnan');

AUC_veh_lev_acq2 = mean(aucFilter_VEHLevAcq2,1,'omitnan');
AUC_jzl8_lev_acq2 = mean(aucFilter_JZL8LevAcq2,1,'omitnan');
AUC_jzl18_lev_acq2 = mean(aucFilter_JZL18LevAcq2,1,'omitnan');

AUC_veh_lev_rev1 = mean(aucFilter_VEHLevRev1,1,'omitnan');
AUC_jzl8_lev_rev1 = mean(aucFilter_JZL8LevRev1,1,'omitnan');
AUC_jzl18_lev_rev1 = mean(aucFilter_JZL18LevRev1,1,'omitnan');

% SEM for bar data %
AUC_veh_cue_acq2_SEM = nanstd(aucFilter_VEHCueAcq2) / sqrt(sum(~isnan(aucFilter_VEHCueAcq2)));
AUC_jzl8_cue_acq2_SEM = nanstd(aucFilter_JZL8CueAcq2) / sqrt(sum(~isnan(aucFilter_JZL8CueAcq2)));
AUC_jzl18_cue_acq2_SEM = nanstd(aucFilter_JZL18CueAcq2) / sqrt(sum(~isnan(aucFilter_JZL18CueAcq2)));

AUC_veh_cue_rev1_SEM = nanstd(aucFilter_VEHCueRev1) / sqrt(sum(~isnan(aucFilter_VEHCueRev1)));
AUC_jzl8_cue_rev1_SEM = nanstd(aucFilter_JZL8CueRev1) / sqrt(sum(~isnan(aucFilter_JZL8CueRev1)));
AUC_jzl18_cue_rev1_SEM = nanstd(aucFilter_JZL18CueRev1) / sqrt(sum(~isnan(aucFilter_JZL18CueRev1)));

AUC_veh_lev_acq2_SEM = nanstd(aucFilter_VEHLevAcq2) / sqrt(sum(~isnan(aucFilter_VEHLevAcq2)));
AUC_jzl8_lev_acq2_SEM = nanstd(aucFilter_JZL8LevAcq2) / sqrt(sum(~isnan(aucFilter_JZL8LevAcq2)));
AUC_jzl18_lev_acq2_SEM = nanstd(aucFilter_JZL18LevAcq2) / sqrt(sum(~isnan(aucFilter_JZL18LevAcq2)));

AUC_veh_lev_rev1_SEM = nanstd(aucFilter_VEHLevRev1) / sqrt(sum(~isnan(aucFilter_VEHLevRev1)));
AUC_jzl8_lev_rev1_SEM = nanstd(aucFilter_JZL8LevRev1) / sqrt(sum(~isnan(aucFilter_JZL8LevRev1)));
AUC_jzl18_lev_rev1_SEM = nanstd(aucFilter_JZL18LevRev1) / sqrt(sum(~isnan(aucFilter_JZL18LevRev1)));

% create arrays for bar data/SEM %
barDataCueAcq2 = [AUC_veh_cue_acq2, AUC_jzl8_cue_acq2, AUC_jzl18_cue_acq2];
barDataCueRev1 = [AUC_veh_cue_rev1, AUC_jzl8_cue_rev1, AUC_jzl18_cue_rev1];
barDataLevAcq2 = [AUC_veh_lev_acq2, AUC_jzl8_lev_acq2, AUC_jzl18_lev_acq2];
barDataLevRev1 = [AUC_veh_lev_rev1, AUC_jzl8_lev_rev1, AUC_jzl18_lev_rev1];

barSEMCueAcq2 = [AUC_veh_cue_acq2_SEM, AUC_jzl8_cue_acq2_SEM, AUC_jzl18_cue_acq2_SEM];
barSEMCueRev1 = [AUC_veh_cue_rev1_SEM, AUC_jzl8_cue_rev1_SEM, AUC_jzl18_cue_rev1_SEM];
barSEMLevAcq2 = [AUC_veh_lev_acq2_SEM, AUC_jzl8_lev_acq2_SEM, AUC_jzl18_lev_acq2_SEM];
barSEMLevRev1 = [AUC_veh_lev_rev1_SEM, AUC_jzl8_lev_rev1_SEM, AUC_jzl18_lev_rev1_SEM];

% Plot bar graphs for cue % 
f3 = figure;
subplot(1,2,1)
barHandle = bar(barDataCueAcq2);
hold on;
numBars = numel(barDataCueAcq2);
xCoords = barHandle.XData;
errorbar(xCoords, barDataCueAcq2, barSEMCueAcq2, 'k', 'LineStyle', 'none');

ylabel('Area Under Curve (+/- SEM)')
title('Phase: Last Acquisition');
xticks(1:numBars);
xticklabels({'Vehicle', 'JZL8', 'JZL18'});
hold off;

subplot(1,2,2)
barHandle = bar(barDataCueRev1);
hold on;
numBars = numel(barDataCueRev1);
xCoords = barHandle.XData;
errorbar(xCoords, barDataCueRev1, barSEMCueRev1, 'k', 'LineStyle', 'none');

ylabel('Area Under Curve (+/- SEM)')
title('Phase: First Reversal');
xticks(1:numBars);
xticklabels({'Vehicle', 'JZL8', 'JZL18'});
hold off;
[ax3,h3] = suplabel(TITLE_Cue,'t');
set(h3,'FontSize',20);
% Plot bar graphs for lever % 
f4 = figure;
subplot(1,2,1)
barHandle = bar(barDataLevAcq2);
hold on;
numBars = numel(barDataLevAcq2);
xCoords = barHandle.XData;
errorbar(xCoords, barDataLevAcq2, barSEMLevAcq2, 'k', 'LineStyle', 'none');

ylabel('Area Under Curve (+/- SEM)')
title('Phase: Last Acquisition');
xticks(1:numBars);
xticklabels({'Vehicle', 'JZL8', 'JZL18'});
hold off;

subplot(1,2,2)
barHandle = bar(barDataLevRev1);
hold on;
numBars = numel(barDataLevRev1);
xCoords = barHandle.XData;
errorbar(xCoords, barDataLevRev1, barSEMLevRev1, 'k', 'LineStyle', 'none');

ylabel('Area Under Curve (+/- SEM)')
title('Phase: First Reversal');
xticks(1:numBars);
xticklabels({'Vehicle', 'JZL8', 'JZL18'});
hold off;
[ax4,h4] = suplabel(TITLE_Lev,'t');
set(h4,'FontSize',20);





toc
NERD_STATS(toc,numFiles);
fprintf('Successfully analyzed error probabilities for %d files\n',numFiles)
