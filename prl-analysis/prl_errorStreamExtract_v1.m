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

    for ii = 1:height(errorProbLeverTS(:,1))
        if isempty(errorProbLeverTS)
            streams_dFF(ii,1:epocArrayLen) = NaN;
            ampdFF(ii) = zeros(1);
            streams_z(ii,1:epocArrayLen) = NaN;
            ampZ(ii) = zeros(1);

            continue
        end
            
        windowStart = errorProbLeverTS(i,1)-baseWindow;
        windowEnd = windowStart+timeWindow+baseWindow;
        [~,windSt] = min(abs(session_time - windowStart));
        [~,windEn] = min(abs(session_time - windowEnd));
        epocSigRaw = SIGNAL_raw(1,windSt:windEn);

        if length(epocSigRaw) < epocArrayLen
            mn = mean(epocSigRaw(1,end-10:end));
            epocSigRaw(1,end:epocArrayLen) = mn;
        elseif length(epocSigRaw) > epocArrayLen
            op = length(epocSigRaw);
            arrayDif = op - epocArrayLen;
            epocSigRaw = epocSigRaw(1,1:end-arrayDif);
        end
        streams_raw(ii,1:epocArrayLen) = epocSigRaw;

    end

    [~,baseSt] = min(abs(ts1 - (baseline(1))));
    [~,baseEn] = min(abs(ts1 - (baseline(2))));

    for j = 1:height(streams_raw)
        % dF/F
        meanBase = mean(streams_raw(j,baseSt:baseEn));
        stdBase = std(streams_raw(j,baseSt:baseEn));
        streams_dFF(j,1:epocArrayLen) = streams_raw(j,1:epocArrayLen) - meanBase;
        streams_dFF(j,1:epocArrayLen) = 100*(streams_dFF(j,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
        stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
        streams_z(j,1:epocArrayLen) = (streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % amplitude
        ampZ(j) = max(streams_z(j,ampSt:ampEn));
    end
end
master_session_STREAMz(i,:) = sessionSTREAM;
master_ampSession(i,:) = sessionAMP;
master_aucSession(i,:) = sessionAUC;

