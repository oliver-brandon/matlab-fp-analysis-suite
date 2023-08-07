clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0 timeWindow]; % time window to grab amplitude from
auc_window = [-1 timeWindow];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
errorType = 1;
lever = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(...
    '/Users/brandon/My Drive/prl/PRL_GRABDA/test','Choose the .mat files you want to analyze.'...
    ); %gets directory%
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
phaseList = {};
treatList = {};

for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    brokenID = strsplit(name,'_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    phaseList{i} = cellstr(strtrim(brokenID{2}));
    treatList{i} = cellstr(strtrim(brokenID{3}));

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
    
    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;

    [sessionSTREAM, sessionAMP, sessionAUC] = epocExtract(...
        SIGNAL_raw,...
        session_time,...
        ts1,...
        session_identifiers(:,1),...
        baseWindow, ...
        timeWindow, ...
        baseline, ...
        amp_window, ...
        epocArrayLen);
    
    [levRewSTREAM, levRewAMP, levRewAUC] = epocExtract(...
        SIGNAL_raw,...
        session_time,...
        ts1,...
        errorProbLeverTS(:,1),...
        baseWindow, ...
        timeWindow, ...
        baseline, ...
        amp_window, ...
        epocArrayLen);

    [levChoiceSTREAM, levChoiceAMP, levChoiceAUC] = epocExtract(...
        SIGNAL_raw,...
        session_time,...
        ts1,...
        errorProbLeverTS(:,2),...
        baseWindow, ...
        timeWindow, ...
        baseline, ...
        amp_window, ...
        epocArrayLen);

    [cueRewSTREAM, cueRewAMP, cueRewAUC] = epocExtract(...
        SIGNAL_raw,...
        session_time,...
        ts1,...
        errorProbCueTS(:,1),...
        baseWindow, ...
        timeWindow, ...
        baseline, ...
        amp_window, ...
        epocArrayLen);

    [cueChoiceSTREAM, cueChoiceAMP, cueChoiceAUC] = epocExtract(...
        SIGNAL_raw,...
        session_time,...
        ts1,...
        errorProbCueTS(:,2),...
        baseWindow, ...
        timeWindow, ...
        baseline, ...
        amp_window, ...
        epocArrayLen);

end
master_session_STREAMz(i,:) = sessionSTREAM;
master_ampSession(i,:) = sessionAMP;
master_aucSession(i,:) = sessionAUC;

master_levRew_STREAMz(i,:) = levRewSTREAM;
master_ampLevRew(i,:) = levRewAMP;
master_aucLevRew(i,:) = levRewAUC;

master_levChoice_STREAMz(i,:) = levChoiceSTREAM;
master_ampLevChoice(i,:) = levChoiceAMP;
master_aucLevChoice(i,:) = levChoiceAUC;

master_cueRew_STREAMz(i,:) = cueRewSTREAM;
master_ampCueRew(i,:) = cueRewAMP;
master_aucLev(i,:) = cueRewAUC;

master_cueChoice_STREAMz(i,:) = cueChoiceSTREAM;
master_ampLevChoice(i,:) = cueChoiceAMP;
master_aucLevChoice(i,:) = cueChoiceAUC;
