clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore
amp_window = [0 2]; % time window to grab amplitude from
auc_window = [-1 timeWindow];
t = 0; % seconds to clip from start of streams
N = 1; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
removeOutlierTrials = 0; % 1 = remove
figsavepath = 'E:\Google Drive\prl\PRL_GRABDA\cueByTrialFigs\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('E:\Google Drive\prl\PRL_GRABDA','Choose the .mat files you want to analyze.'); %gets directory%
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
    brokenID = strsplit(name,'_');
    tempID = cellstr(brokenID(1));
    tempPhase = cellstr(brokenID(2));
    tempTreat = cellstr(brokenID(3));
    IDs = vertcat(IDs,tempID);
    treatList = vertcat(treatList,tempTreat);
    prl_phase = vertcat(prl_phase,tempPhase);
    load(filename)
    TITLE = strcat(tempID,{' '},tempPhase,{' '},tempTreat);
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
    [~,ampSt] = min(abs(ts1 - amp_window(1)));
    [~,ampEn] = min(abs(ts1 - amp_window(2)));
    [~,aucSt] = min(abs(ts1 - auc_window(1)));
    [~,aucEn] = min(abs(ts1 - auc_window(2)));

    cueArray = session_identifiers(1:2:end-1,:);
    leverArray = session_identifiers(2:2:end,:);

    cuelever_cRew = [];
    cuelever_cNoRew = [];
    cuelever_iRew = [];
    cuelever_iNoRew = [];

    cueTT_dFF = zeros(height(cueArray),epocArrayLen);
    cueTT_z = zeros(height(cueArray),epocArrayLen);
    cueTT_raw = zeros(height(cueArray),epocArrayLen);

    cuelever_Correct = [];
    cuelever_Incorrect = [];

    for m = 1:height(cueArray)
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
    idx = find(ts1>-0.5,1);
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
       
        if cueTT_z(n,idx) < 0
            val = cueTT_z(n,idx);
            diff = 0 - val;
            cueTT_z(n,1:epocArrayLen) = cueTT_z(n,1:epocArrayLen) + abs(diff);
        elseif cueTT_z(n,idx) > 0
            val = cueTT_z(n,idx);
            diff = 0 - val;
            cueTT_z(n,1:epocArrayLen) = cueTT_z(n,1:epocArrayLen) - abs(diff);
        end
        cueAMP(n,1) = max(cueTT_z(n,ampSt:ampEn));
        if removeOutlierTrials == 1
            
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
        elseif leverArray(n,2) == 2
            cuelever_cNoRew = [cuelever_cNoRew;cueTT_z(n,:)];
            cuelever_Correct = [cuelever_Correct;cueTT_z(n,:)];
        elseif leverArray(n,2) == 3
            cuelever_iRew = [cuelever_iRew;cueTT_z(n,:)];
            cuelever_Incorrect = [cuelever_Incorrect;cueTT_z(n,:)];
        elseif leverArray(n,2) == 4
            cuelever_iNoRew = [cuelever_iNoRew;cueTT_z(n,:)];
            cuelever_Incorrect = [cuelever_Incorrect;cueTT_z(n,:)];
        end
    end
    for jj = 1:length(leverArray)
        if leverArray(jj,2) == 1 || leverArray(jj,2) == 2
            trialVal(jj,1) = 1;
        else
            trialVal(jj,1) = 0;
        end
    end
    x = 1:1:length(cueAMP);
    f1 = figure;
    title(TITLE)
    xlabel('Trial')
    yyaxis left
    plot(x(2:end),cueAMP(2:end,1),'r')
    ylabel('Cue GrabDA Amp (zScore)')
    ylim([-5,15]);
    set(gca,'YColor','black','Box', 'on','Color','w')
    set(get(gca, 'YLabel'), 'Color','black') % Rotate the right ylabel
    
    yyaxis right
    plot(x(2:end),trialVal(2:end,1),'Color',[0.5,0.5,0.5])
    ylim([-0.5,1.5]);
    yticks([0,1]);
    ylabel('Choice (Correct = 1')
    set(gca,'YColor','black','Box', 'on','Color','w')
    set(get(gca, 'YLabel'), 'Rotation', -90, 'Color','black') % Rotate the right ylabel
    file_name1 = char(strcat(figsavepath,TITLE,' Cue by Trial','.pdf'));
    print(f1,file_name1,'-dpdf','-vector','-bestfit');
    close all
end
toc
NERD_STATS(toc,numFiles);
