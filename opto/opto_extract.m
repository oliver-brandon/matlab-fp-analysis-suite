clear; close all;
warning off

%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % baseline signal to include before TTL 
baseline = [-2 -1]; % baseline signal for dFF/zscore
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
adjustBase = 1; % 1 = yes, 2 = no
baseAdjust = -2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = '/Users/brandon/personal-drive/optomouse-prime/opto-reversal/correct-stim/mat-corstim-SNc/1373M_Rev5CorSham.mat';
load(file)


ISOS = 'x405A';
SIGNAL = 'x465A';
cue = data.epocs.St1_.onset;
correct = data.epocs.CL1_.onset;
incorrect = data.epocs.IL1_.onset;
levers = [correct;incorrect];
levers = sort(levers);
epoc = levers;
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
for m = 1:height(epoc)
    epocTTStart = epoc(m,1) - baseWindow;
    if epocTTStart < t
        continue
    end
    cueTTEnd = epocTTStart + timeWindow + baseWindow;
    [~,cTTSt] = min(abs(session_time - epocTTStart));
    [~,cTTEn] = min(abs(session_time - cueTTEnd));
    epocTTSigRaw = SIGNAL_raw(1,cTTSt:cTTEn);
    if length(epocTTSigRaw) < epocArrayLen
        mn = mean(epocTTSigRaw(1,end-10:end));
        epocTTSigRaw(1,end:epocArrayLen) = mn;
    elseif length(epocTTSigRaw) > epocArrayLen
        op = length(epocTTSigRaw);
        arrayDif = op - epocArrayLen;
        epocTTSigRaw = epocTTSigRaw(1,1:end-arrayDif);
    end
    epocTT_raw(m,:) = epocTTSigRaw;
                
end
for n = 1:height(epoc)
    % dF/F
    meanBase = mean(epocTT_raw(n,baseSt:baseEn));
    stdBase = std(epocTT_raw(n,baseSt:baseEn));
    
    epocTT_dFF(n,1:epocArrayLen) = epocTT_raw(n, 1:epocArrayLen) - meanBase;
    epocTT_dFF(n,1:epocArrayLen) = 100*(epocTT_dFF(n,1:epocArrayLen) / meanBase);
    % z-Score
    meanBase_dFF = mean(epocTT_dFF(n,baseSt:baseEn));
    stdBase_dFF = std(epocTT_dFF(n,baseSt:baseEn));
    epocTT_z(n,1:epocArrayLen) = (epocTT_dFF(n,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
    if adjustBase == 1
        if epocTT_z(n,idx) < 0
            val = epocTT_z(n,idx);
            diff = 0 - val;
            epocTT_z(n,1:epocArrayLen) = epocTT_z(n,1:epocArrayLen) + abs(diff);
        elseif epocTT_z(n,idx) > 0
            val = epocTT_z(n,idx);
            diff = 0 - val;
            epocTT_z(n,1:epocArrayLen) = epocTT_z(n,1:epocArrayLen) - abs(diff);
        end
    else
        disp("")
    end
end