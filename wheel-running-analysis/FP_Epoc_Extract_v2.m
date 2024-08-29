clear;
VERSION = 2.0;
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % baseline signal to include before TTL 
baseline = [-2 -1]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
filter = 3; %finds runStart times that occur at least 'filt' seconds after onWheel
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
baseAdjust = -2; % adjust baseline of signals to 'baseAdjust' seconds
DLS_isos = 'x405A';
DLS = 'x465A';
NAc_isos = 'x405C';
NAc = 'x465C';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("VERSION: %d\n",VERSION)
[filename,pathname] = uigetfile('/Users/brandon/personal-drive/hot-wheels/wheel-running-mats/unlocked');
% Check if the user selected a file or canceled the dialog
if isequal(filename, 0) || isequal(pathname, 0)
    disp('User canceled the file selection.');
    return
else
    % Full path to the selected file
    fullpath = fullfile(pathname, filename);
    % Load the MAT file
    load(fullpath)
end
tic
numFiles = 1;
[~,name,~] = fileparts(filename);
fprintf('Analyzing %s\n',name);
brokenID = strsplit(name,'_');
ID = cellstr(strtrim(brokenID{1}));
day = cellstr(strtrim(brokenID{3}));

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
if isempty(oWfiltered) || isempty(rSfiltered)
    disp('No valid filtered epocs.');
    return
end
%time array used for all streams%
session_time = (1:length(data.streams.(DLS).data))/data.streams.(DLS).fs;
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


session_time = downsample(session_time, N);
ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(DLS).fs*N;

% establish baseline windows
[~,baseSt] = min(abs(ts1 - (baseline(1))));
[~,baseEn] = min(abs(ts1 - (baseline(2))));
base_idx = find(ts1>baseAdjust,1);
epocList = {oWfiltered;rSfiltered};
%% Streams baslined to before each epoc %%
for k = 1:numel(epocList)
    epoc = double(epocList{k});
    dls_streams_raw = [];
    dls_streams_dFF = [];
    dls_streams_z = [];
    nac_streams_raw = [];
    nac_streams_dFF = [];
    nac_streams_z = [];

    for ii = 1:height(epoc)
        if epoc(ii) < t
            continue
        else
            windowStart = epoc(ii)-baseWindow;
            windowEnd = windowStart+timeWindow+baseWindow;
            [~,windSt] = min(abs(session_time - windowStart));
            [~,windEn] = min(abs(session_time - windowEnd));
            epocSigRaw_DLS = DLS_SIGNAL_raw(1,windSt:windEn);
            epocSigRaw_NAc = NAc_SIGNAL_raw(1,windSt:windEn);
            dls_streams_raw(ii,:) = epocSigRaw_DLS(:,1:epocArrayLen);
            nac_streams_raw(ii,:) = epocSigRaw_NAc(:,1:epocArrayLen);
        end
    end

    for j = 1:height(dls_streams_raw)
        % DLS
        % dF/F
        meanBase = mean(dls_streams_raw(j,baseSt:baseEn));
        stdBase = std(dls_streams_raw(j,baseSt:baseEn));
        dls_streams_dFF(j,1:epocArrayLen) = dls_streams_raw(j,1:epocArrayLen) - meanBase;
        dls_streams_dFF(j,1:epocArrayLen) = 100*(dls_streams_dFF(j,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(dls_streams_dFF(j,baseSt:baseEn));
        stdBase_dFF = std(dls_streams_dFF(j,baseSt:baseEn));
        dls_streams_z(j,1:epocArrayLen) = (dls_streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % adjusts streams to baseline of zero %
        if dls_streams_z(j,base_idx) < 0
            val = dls_streams_z(j,base_idx);
            diff = 0 - val;
            dls_streams_z(j,1:epocArrayLen) = dls_streams_z(j,1:epocArrayLen) + abs(diff);
        elseif dls_streams_z(j,base_idx) > 0
            val = dls_streams_z(j,base_idx);
            diff = 0 - val;
            dls_streams_z(j,1:epocArrayLen) = dls_streams_z(j,1:epocArrayLen) - abs(diff);
        end

        % NAc
        % dF/F
        meanBase = mean(nac_streams_raw(j,baseSt:baseEn));
        stdBase = std(nac_streams_raw(j,baseSt:baseEn));
        nac_streams_dFF(j,1:epocArrayLen) = nac_streams_raw(j,1:epocArrayLen) - meanBase;
        nac_streams_dFF(j,1:epocArrayLen) = 100*(nac_streams_dFF(j,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(nac_streams_dFF(j,baseSt:baseEn));
        stdBase_dFF = std(nac_streams_dFF(j,baseSt:baseEn));
        nac_streams_z(j,1:epocArrayLen) = (nac_streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % adjusts streams to baseline of zero %
        if nac_streams_z(j,base_idx) < 0
            val = nac_streams_z(j,base_idx);
            diff = 0 - val;
            nac_streams_z(j,1:epocArrayLen) = nac_streams_z(j,1:epocArrayLen) + abs(diff);
        elseif nac_streams_z(j,base_idx) > 0
            val = nac_streams_z(j,base_idx);
            diff = 0 - val;
            nac_streams_z(j,1:epocArrayLen) = nac_streams_z(j,1:epocArrayLen) - abs(diff);
        end
        
        dls_outputSTREAMSz{k,1} = dls_streams_z(:,1:epocArrayLen);
        nac_outputSTREAMSz{k,1} = nac_streams_z(:,1:epocArrayLen);
    end
end 

dls_onWheel = dls_outputSTREAMSz{1,1};
dls_runStart = dls_outputSTREAMSz{2,1};
nac_onWheel = nac_outputSTREAMSz{1,1};
nac_runStart = nac_outputSTREAMSz{2,1};

epoc_store.dls.onWheel = dls_onWheel;
epoc_store.dls.runStart = dls_runStart;
epoc_store.nac.onWheel = nac_onWheel;
epoc_store.nac.runStart = nac_runStart;
epoc_store.time = ts1;

toc
NERD_STATS(toc,numFiles);
clearvars -except epoc_store