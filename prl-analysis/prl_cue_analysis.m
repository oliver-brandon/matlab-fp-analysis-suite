clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore
amp_window = [0 2]; % time window to grab amplitude from
auc_window = [-1 timeWindow];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/DA_PRL','Choose the .mat files you want to analyze.'); %gets directory%
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
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    [~,treatment,~] = fileparts(myDir);
    brokenID = strsplit(name,'_');
    tempID = cellstr(brokenID(1));
    tempTreat = cellstr(treatment);
    IDs = vertcat(IDs,tempID);
    treatList = vertcat(treatList,tempTreat);
    prl_phase = char(brokenID(2));
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
    % some files end with a cue (0) instead of lever press %
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(3:end-1,:);
        cueArray = cueArray(1:end-1,1);
    end

    for m = 1:height(cueArray)
        cueTTStart = cueArray(m,1) - baseWindow;
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
        
        % sorts the cue psth by the next trial type %
        if leverArray(n,2) == 1
            cuelever_cRew = [cuelever_cRew;cueTT_z(n,:)];
        elseif leverArray(n,2) == 2
            cuelever_cNoRew = [cuelever_cNoRew;cueTT_z(n,:)];
        elseif leverArray(n,2) == 3
            cuelever_iRew = [cuelever_iRew;cueTT_z(n,:)];
        elseif leverArray(n,2) == 4
            cuelever_iNoRew = [cuelever_iNoRew;cueTT_z(n,:)];
        end
    end

    if isempty(cuelever_cRew)
        cuelever_cRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_cRewAMP(i,1) = nan;
        cuelever_cRewAUC(i,1) = nan;
    else
        cuelever_cRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_cRew);
        cuelever_cRewAMP(i,1) = max(cuelever_cRewSTREAMSz(i,ampSt:ampEn));
        cuelever_cRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),cuelever_cRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(cuelever_cNoRew)
        cuelever_cNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_cNoRewAMP(i,1) = nan;
        cuelever_cNoRewAUC(i,1) = nan;
    else
        cuelever_cNoRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_cNoRew);
        cuelever_cNoRewAMP(i,1) = max(cuelever_cNoRewSTREAMSz(i,ampSt:ampEn));
        cuelever_cNoRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),cuelever_cNoRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(cuelever_iRew)
        cuelever_iRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_iRewAMP(i,1) = nan;
        cuelever_iRewAUC(i,1) = nan;
    else
        cuelever_iRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_iRew);
        cuelever_iRewAMP(i,1) = max(cuelever_iRewSTREAMSz(i,ampSt:ampEn));
        cuelever_iRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),cuelever_iRewSTREAMSz(i,aucSt:aucEn));
    end
    if isempty(cuelever_iNoRew)
        cuelever_iNoRewSTREAMSz(i,1:epocArrayLen) = nan;
        cuelever_iNoRewAMP(i,1) = nan;
        cuelever_iNoRewAUC(i,1) = nan;
    else
        cuelever_iNoRewSTREAMSz(i,1:epocArrayLen) = mean(cuelever_iNoRew);
        cuelever_iNoRewAMP(i,1) = max(cuelever_iNoRewSTREAMSz(i,ampSt:ampEn));
        cuelever_iNoRewAUC(i,1) = trapz(ts1(1,aucSt:aucEn),cuelever_iNoRewSTREAMSz(i,aucSt:aucEn));
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
%% cue before lever table %%
cuelever_cRewTable = array2table(cuelever_cRewSTREAMSz);
cuelever_cRewTable = horzcat(IDs,cuelever_cRewTable);
cuelever_cNoRewTable = array2table(cuelever_cNoRewSTREAMSz);
cuelever_cNoRewTable = horzcat(IDs,cuelever_cNoRewTable);
cuelever_iRewTable = array2table(cuelever_iRewSTREAMSz);
cuelever_iRewTable = horzcat(IDs,cuelever_iRewTable);
cuelever_iNoRewTable = array2table(cuelever_iNoRewSTREAMSz);
cuelever_iNoRewTable = horzcat(IDs,cuelever_iNoRewTable);

cueleverTTAMP = [cuelever_cRewAMP cuelever_cNoRewAMP cuelever_iRewAMP cuelever_iNoRewAMP];
cueleverTTAMP = array2table(cueleverTTAMP,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
cueleverTTAMP = horzcat(IDs,cueleverTTAMP);
cueleverTTAUC = [cuelever_cRewAUC cuelever_cNoRewAUC cuelever_iRewAUC cuelever_iNoRewAUC];
cueleverTTAUC = array2table(cueleverTTAUC,'VariableNames',{'cue_cRew','cue_cNoRew','cue_iRew','cue_iNoRew'});
cueleverTTAUC = horzcat(IDs,cueleverTTAUC);

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
sgtitle(sprintf('Cue Before Trial - %s\n',treatment))

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
sgtitle(sprintf('Cue After Trial - %s\n',treatment))
toc
NERD_STATS(toc,numFiles);