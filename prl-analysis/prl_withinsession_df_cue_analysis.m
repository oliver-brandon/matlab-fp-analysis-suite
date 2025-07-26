clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore
zeroto = -0.5;
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
dualFiberChannel = 2; % 1 = dual fiber channel A, 2 = dual fiber channel C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(...
    '/Users/brandon/personal-drive/prl/GrabDA/dls-nac/within-session-reversal/good-performance',...
    'Choose the .mat files you want to analyze.'); %gets directory%
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
probability = {};
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    fprintf('Analyzing %s (%d of %d)\n', name, i, numFiles)
    brokenID = strsplit(name,'_');
    tempID = cellstr(brokenID(1));
    tempProb = cellstr(brokenID(2));
    IDs = vertcat(IDs,tempID);
    probability = vertcat(probability,tempProb);
    load(filename)
    if dualFiberChannel == 1
        ISOS = 'x405A';
        SIGNAL = 'x465A';
    elseif dualFiberChannel == 2
        ISOS = 'x405C';
        SIGNAL = 'x465C';
    end
    if isfield(data.epocs, 'St1_')
        cue = data.epocs.St1_.onset;
    else
        disp('Cannot locate TTL')
        break
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

    [~,baseSt] = min(abs(ts1 - baseline(1)));
    [~,baseEn] = min(abs(ts1 - baseline(2)));

    cue_raw = zeros(height(cue),epocArrayLen);
    cue_dFF = zeros(height(cue),epocArrayLen);
    cue_z = zeros(height(cue),epocArrayLen);

    for m = 1:height(cue)
        cueStart = cue(m,1) - baseWindow;
        % if cueStart < t
        %     continue
        % end
        cueEnd = cueStart + timeWindow + baseWindow;
        [~,cSt] = min(abs(session_time - cueStart));
        [~,cEn] = min(abs(session_time - cueEnd));
        cueSigRaw = SIGNAL_raw(1,cSt:cEn);
        if length(cueSigRaw) < epocArrayLen
            mn = mean(cueSigRaw(1,end-10:end));
            cueSigRaw(1,end:epocArrayLen) = mn;
        elseif length(cueSigRaw) > epocArrayLen
            op = length(cueSigRaw);
            arrayDif = op - epocArrayLen;
            cueSigRaw = cueSigRaw(1,1:end-arrayDif);
        end
        cue_raw(m,:) = cueSigRaw;
                
    end
    idx = find(ts1>zeroto,1);
    for n = 1:height(cue_raw)
        % dF/F
        meanBase = mean(cue_raw(n,baseSt:baseEn));
        stdBase = std(cue_raw(n,baseSt:baseEn));
        
        cue_dFF(n,1:epocArrayLen) = cue_raw(n, 1:epocArrayLen) - meanBase;
        cue_dFF(n,1:epocArrayLen) = 100*(cue_dFF(n,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(cue_dFF(n,baseSt:baseEn));
        stdBase_dFF = std(cue_dFF(n,baseSt:baseEn));
        cue_z(n,1:epocArrayLen) = (cue_dFF(n,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
       
        if cue_z(n,idx) < 0
            val = cue_z(n,idx);
            diff = 0 - val;
            cue_z(n,1:epocArrayLen) = cue_z(n,1:epocArrayLen) + abs(diff);
        elseif cue_z(n,idx) > 0
            val = cue_z(n,idx);
            diff = 0 - val;
            cue_z(n,1:epocArrayLen) = cue_z(n,1:epocArrayLen) - abs(diff);
        end
    end
    
    cueAcq_firstFive(:,i) = mean(cue_z(1:5,:))';
    cueAcq_lastFive(:,i) = mean(cue_z(26:30,:))';
    cueRev_firstFive(:,i) = mean(cue_z(31:35,:))';
    cueRev_lastFive(:,i) = mean(cue_z(end-4:end,:))';

end

IDs = cell2table(IDs);
IDs = rows2vars(IDs);
probability = cell2table(probability);
probability = rows2vars(probability);
cueAcq_firstFive_Table = array2table(cueAcq_firstFive);

cueAcq_lastFive_Table = array2table(cueAcq_lastFive);

cueRev_firstFive_Table = array2table(cueRev_firstFive);

cueRev_lastFive_Table = array2table(cueRev_lastFive);




toc
NERD_STATS(toc,numFiles);
clearvars -except cueAcq_firstFive_Table cueAcq_lastFive_Table ...
    cueRev_firstFive_Table cueRev_lastFive_Table IDs probability