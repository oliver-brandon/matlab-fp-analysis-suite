clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = 2.0;
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-5 -3]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0.5 1]; % time window to grab amplitude from
auc_window = [0 5];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
filter = 3; %finds runStart times that occur at least 'filt' seconds after onWheel
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
baseAdjust = -5; % adjust baseline of signals to 'baseAdjust' seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("VERSION: %d\n",VERSION)

myDir = uigetdir('/Users/brandon/personal-drive/hot-wheels/wheel-running-mats/test'); %gets directory%
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name}, {'.mat'}));
numFiles = length(myFiles);
IDs = cell(size(numFiles,1));
dayList = cell(size(numFiles,1));
fprintf("Starting batch extraction of %d files...\n",numFiles)
for i = 1:numFiles
    filename = fullfile(myDir, myFiles(i).name);
    [~,name,~] = fileparts(filename);
    fprintf('Analyzing %s (%d of %d)\n', name, i, numFiles)
    load(filename)
    

    brokenID = strsplit(name,'_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    dayList{i} = cellstr(strtrim(brokenID{3}));
    
    DLS_isos = 'x405A';
    DLS = 'x465A';
    NAc_isos = 'x405C';
    NAc = 'x465C';

    runStart = data.epocs.runStart.onset;
    runStop = data.epocs.runStop.onset;
    onWheel = data.epocs.onWheel.onset;
    offWheel = data.epocs.onWheel.onset;
    rSfiltered = [];
    oWfiltered = [];
    for j = 1:length(runStart)
        rS = runStart(j);
        idx = find(onWheel<(rS-filter),1,'last');
        oWfiltered = [oWfiltered;onWheel(idx)];
      
    end
    oWfiltered = unique(oWfiltered,'stable');
    for k = 1:length(oWfiltered)
        oW = oWfiltered(k);
        idx = find(runStart>(oW+filter),1,'first');
        rSfiltered = [rSfiltered;runStart(idx)];
    end

    %time array used for all streams%
    session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    ind = find(session_time>t,1);% find first index of when time crosses threshold
    session_time = session_time(ind:end); % reformat vector to only include allowed time

    DLS_SIGNAL_raw = data.streams.(DLS).data(ind:end);
    DLS_ISOS_raw = data.streams.(DLS_isos).data(ind:end);
    NAc_SIGNAL_raw = data.streams.(NAc).data(ind:end);
    NAc_ISOS_raw = data.streams.(NAc_isos).data(ind:end);

    %downsample streams and time array by N times%
    DLS_ISOS_raw = downsample(DLS_ISOS_raw, N);
    DLS_SIGNAL_raw = downsample(DLS_SIGNAL_raw, N);
    NAc_ISOS_raw = downsample(NAc_ISOS_raw, N);
    NAc_SIGNAL_raw = downsample(NAc_SIGNAL_raw, N);

    minStreamLength = min(length(DLS_ISOS_raw),length(DLS_SIGNAL_raw));
    DLS_ISOS_raw = DLS_ISOS_raw(1:minStreamLength);
    DLS_SIGNAL_raw = DLS_SIGNAL_raw(1:minStreamLength);

    minStreamLength = min(length(NAc_ISOS_raw),length(NAc_SIGNAL_raw));
    NAc_ISOS_raw = NAc_ISOS_raw(1:minStreamLength);
    NAc_SIGNAL_raw = NAc_SIGNAL_raw(1:minStreamLength);

    SIGNAL_raw = data.streams.(SIGNAL).data(ind:end);
    ISOS_raw = data.streams.(ISOS).data(ind:end);

    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(DLS).fs*N;

    % establish baseline windows
    [~,baseSt] = min(abs(ts1 - (baseline(1))));
    [~,baseEn] = min(abs(ts1 - (baseline(2))));
    
    epocList = {};
    %% Streams baslined to before each epoc %%
    for k = 1:numel(epocList)
        epoc = double(epocList{k});
        streams_raw = double(outputSTREAMSraw{k,1});
        streams_dFF = double(outputSTREAMSdFF{k,1});
        streams_z = double(outputSTREAMSz{k,1});
        ampdFF = double(outputAMPdFF{k,1});
        ampZ = double(outputAMPz{k,1});
        aucdFF = double(outputAUCdFF{k,1});
        aucZ = double(outputAUCz{k,1});
    for ii = 1:height(epoc)
        if epoc(ii) == 0
            streams_dFF(ii,1:epocArrayLen) = NaN;
            ampdFF(ii) = NaN;
            streams_z(ii,1:epocArrayLen) = NaN;
            ampZ(ii) = NaN;
            aucdFF(ii) = NaN;
            aucZ(ii) = NaN;
            continue
        end

        windowStart = epoc(ii)-baseWindow;
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
        % adjusts streams to baseline of zero at -0.5s %
        if streams_z(j,idx) < 0
            val = streams_z(j,idx);
            diff = 0 - val;
            streams_z(j,1:epocArrayLen) = streams_z(j,1:epocArrayLen) + abs(diff);
        elseif streams_z(j,idx) > 0
            val = streams_z(j,idx);
            diff = 0 - val;
            streams_z(j,1:epocArrayLen) = streams_z(j,1:epocArrayLen) - abs(diff);
        end

        % amplitude
        ampdFF(j) = max(streams_dFF(j,ampSt:ampEn));
        ampZ(j) = max(streams_z(j,ampSt:ampEn));

        % % Calculate AUC above x=0 (Z-score) %
        % positive_indices = streams_dFF(j,:) > 0;
        % y_pos = streams_dFF(j,positive_indices);
        % x_pos = ts1(1,positive_indices);
        % aucZ(j) = trapz(x_pos,y_pos);
        aucZ(j) = calculateAUC(streams_z(j,aucSt:aucEn),ts1(1,aucSt:aucEn));
        % % Calculate AUC above x=0 (dFF) %
        % positive_indices = streams_z(j,:) > 0;
        % y_pos = streams_z(j,positive_indices);
        % x_pos = ts1(1,positive_indices);
        % aucdFF(j) = trapz(x_pos,y_pos);
        % aucdFF(j) = trapz(streams_dFF(j,aucSt:aucEn),ts1(1,aucSt:aucEn));
        
    
        outputSTREAMSraw{k,1} = streams_raw(:,1:epocArrayLen);
        outputSTREAMSdFF{k,1} = streams_dFF(:,1:epocArrayLen);
        outputSTREAMSz{k,1} = streams_z(:,1:epocArrayLen);
        outputAMPdFF{k,1} = ampdFF;
        outputAMPz{k,1} = ampZ;
        outputAUCdFF{k,1} = aucdFF;
        outputAUCz{k,1} = aucZ;
    end
    end   
end

fprintf("Successfully extracted streams from %d files(s)\n",numFiles)
disp("The DLS signal(s) is stored in DLS_epoc_store")
disp("The NAc signal(s) is stored in NAc_epoc_store")