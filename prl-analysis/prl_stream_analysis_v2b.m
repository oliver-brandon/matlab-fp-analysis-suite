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
toPlot = 0; % 1 = plot figures, 0 = don't plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(...
    pwd,'Choose the .mat files you want to analyze.'...
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

AMPdFF_analysis = cell(numFiles,5);
AMPz_analysis = cell(numFiles,5);
AUCdFF_analysis = cell(numFiles,5);
AUCz_analysis = cell(numFiles,5);
master_cue_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cNoRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_iRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_iNoRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cue_STREAMz = zeros(numFiles,epocArrayLen);
master_cRew_STREAMz = zeros(numFiles,epocArrayLen);
master_cNoRew_STREAMz = zeros(numFiles,epocArrayLen);
master_iRew_STREAMz = zeros(numFiles,epocArrayLen);
master_iNoRew_STREAMz = zeros(numFiles,epocArrayLen);

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

        cueSTREAMraw = zeros(height(cue),epocArrayLen);
        cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);

        cueSTREAMdFF = zeros(height(cue),epocArrayLen);
        cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
        
        cueSTREAMz = zeros(height(cue),epocArrayLen);
        cRewSTREAMz = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMz = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);

        cueAMPdFF = zeros(height(cue),1);
        cRewAMPdFF = zeros(height(cRew),1);
        cNoRewAMPdFF = zeros(height(cNoRew),1);
        iRewAMPdFF = zeros(height(iRew),1);
        iNoRewAMPdFF = zeros(height(iNoRew),1);

        cueAMPz = zeros(height(cue),1);
        cRewAMPz = zeros(height(cRew),1);
        cNoRewAMPz = zeros(height(cNoRew),1);
        iRewAMPz = zeros(height(iRew),1);
        iNoRewAMPz = zeros(height(iNoRew),1);

        cueAUCdFF = zeros(height(cue),1);
        cRewAUCdFF = zeros(height(cRew),1);
        cNoRewAUCdFF = zeros(height(cNoRew),1);
        iRewAUCdFF = zeros(height(iRew),1);
        iNoRewAUCdFF = zeros(height(iNoRew),1);
        
        cueAUCz = zeros(height(cue),1);
        cRewAUCz = zeros(height(cRew),1);
        cNoRewAUCz = zeros(height(cNoRew),1);
        iRewAUCz = zeros(height(iRew),1);
        iNoRewAUCz = zeros(height(iNoRew),1);



        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
        epocList = {cue;cRew;cNoRew;iRew;iNoRew;session_identifiers};
        outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
            iRewSTREAMraw;iNoRewSTREAMraw};
        outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
            iRewSTREAMdFF;iNoRewSTREAMdFF};
        outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
            iRewSTREAMz;iNoRewSTREAMz};
        outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF};
        outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz};
        outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF};
        outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz};
         
        data.analysis.sessionID = session_identifiers;
        
    elseif isfield(data.streams, 'x405C')
        ISOS = 'x405C';
        SIGNAL = 'x465C';

        cue = data.epocs.St2_.onset;
        cRew = data.epocs.cRewC.onset;
        cNoRew = data.epocs.cNoRewC.onset;
        iRew = data.epocs.iRewC.onset;
        iNoRew = data.epocs.iNoRewC.onset;

        cueSTREAMraw = zeros(height(cue),epocArrayLen);
        cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);

        cueSTREAMdFF = zeros(height(cue),epocArrayLen);
        cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
        
        cueSTREAMz = zeros(height(cue),epocArrayLen);
        cRewSTREAMz = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMz = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);

        cueAMPdFF = zeros(height(cue),1);
        cRewAMPdFF = zeros(height(cRew),1);
        cNoRewAMPdFF = zeros(height(cNoRew),1);
        iRewAMPdFF = zeros(height(iRew),1);
        iNoRewAMPdFF = zeros(height(iNoRew),1);
        
        cueAMPz = zeros(height(cue),1);
        cRewAMPz = zeros(height(cRew),1);
        cNoRewAMPz = zeros(height(cNoRew),1);
        iRewAMPz = zeros(height(iRew),1);
        iNoRewAMPz = zeros(height(iNoRew),1);

        cueAUCdFF = zeros(height(cue),1);
        cRewAUCdFF = zeros(height(cRew),1);
        cNoRewAUCdFF = zeros(height(cNoRew),1);
        iRewAUCdFF = zeros(height(iRew),1);
        iNoRewAUCdFF = zeros(height(iNoRew),1);
        
        cueAUCz = zeros(height(cue),1);
        cRewAUCz = zeros(height(cRew),1);
        cNoRewAUCz = zeros(height(cNoRew),1);
        iRewAUCz = zeros(height(iRew),1);
        iNoRewAUCz = zeros(height(iNoRew),1);

        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
        epocList = {cue;cRew;cNoRew;iRew;iNoRew;session_identifiers};
        outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
            iRewSTREAMraw;iNoRewSTREAMraw};
        outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
            iRewSTREAMdFF;iNoRewSTREAMdFF};
        outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
            iRewSTREAMz;iNoRewSTREAMz};
        outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF};
        outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz};
        outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF};
        outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz};
        
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
    % establish baseline windows
    [~,baseSt] = min(abs(ts1 - (baseline(1))));
    [~,baseEn] = min(abs(ts1 - (baseline(2))));
    [~,ampSt] = min(abs(ts1 - (amp_window(1))));
    [~,ampEn] = min(abs(ts1 - (amp_window(2))));
    [~,aucSt] = min(abs(ts1 - (auc_window(1))));
    [~,aucEn] = min(abs(ts1 - (auc_window(2))));
    
    %% Streams baselined to cue preceding it %%
    cueArray = session_identifiers(1:2:end,:);
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(1:end-1,:);
        cueArray = cueArray(1:end-1,1);
    end
    leverArray = session_identifiers(2:2:end,:);
    levers_z = zeros(height(leverArray),epocArrayLen);
    levers_raw = zeros(height(leverArray),epocArrayLen);
    for m = 1:height(leverArray)
        
        cueBase1 = cueArray(m,1) - 3;
        cueBase2 = cueArray(m,1) - 1;
        [~,cueSt] = min(abs(session_time - cueBase1));
        [~,cueEn] = min(abs(session_time - cueBase2));
        cueBaseMean(m,1) = mean(SIGNAL_raw(1,cueSt:cueEn));
        cueBaseStd(m,1) = std(SIGNAL_raw(1,cueSt:cueEn));

        leverStart = leverArray(m,1) - baseWindow;
        leverEnd = leverStart + timeWindow + baseWindow;
        [~,levSt] = min(abs(session_time - leverStart));
        [~,levEn] = min(abs(session_time - leverEnd));
        leverSigRaw = SIGNAL_raw(1,levSt:levEn);
        if length(leverSigRaw) < epocArrayLen
            mn = mean(leverSigRaw(1,end-10:end));
            leverSigRaw(1,end:epocArrayLen) = mn;
        elseif length(leverSigRaw) > epocArrayLen
            op = length(leverSigRaw);
            arrayDif = op - epocArrayLen;
            leverSigRaw = leverSigRaw(1,1:end-arrayDif);
        end
        levers_raw(m,:) = leverSigRaw;
                
    end
    sessionSTREAMSraw{i,1} = levers_raw;
    amp_cueBase = [];
    auc_cueBase = [];
    for n = 1:height(levers_raw)
        % dF/F
        meanCue = cueBaseMean(n,1);
        stdCue = cueBaseStd(n,1);
        
        levers_dFF(n,1:epocArrayLen) = levers_raw(n, 1:epocArrayLen) - meanCue;
        levers_dFF(n,1:epocArrayLen) = 100*(levers_dFF(n,1:epocArrayLen) / meanCue);
        % z-Score
        meanCueBase_dFF = mean(levers_dFF(n,baseSt:baseEn));
        stdCueBase_dFF = std(levers_dFF(n,baseSt:baseEn));
        levers_z(n,1:epocArrayLen) = (levers_dFF(n,1:epocArrayLen) - meanCueBase_dFF) / stdCueBase_dFF;
        amp_cueBase(n,1) = max(levers_z(n,ampSt:ampEn));
        auc_cueBase(n,1) = abs(trapz(ts1(1,aucSt:aucEn),levers_z(n,aucSt:aucEn)));
    end
%     firstLever(i,1:epocArrayLen) = levers_z(1,:);
%     lastLever(i,1:epocArrayLen) = levers_z(end,:);
    sessionSTREAMSz{i,1} = mean(levers_z);
    trialNames = array2table(session_identifiers(2:2:end,2),'VariableNames',{'Trial_Type'});
    levers_z = array2table(levers_z);
    levers_z = horzcat(trialNames,levers_z);
    cRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 1, 2:end));
    cNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 2, 2:end));
    iRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 3, 2:end));
    iNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 4, 2:end));

    amp_cueBase = array2table(amp_cueBase);
    amp_cueBase = horzcat(trialNames,amp_cueBase);
    auc_cueBase = array2table(auc_cueBase);
    auc_cueBase = horzcat(trialNames,auc_cueBase);

    if ~isempty(cRew_cueBase) && height(cRew_cueBase) > 1
        firstcRew(i,1:epocArrayLen) = cRew_cueBase(1,:);
        firstcRew_AMP(i,1) = max(cRew_cueBase(1,ampSt:ampEn));
        firstcRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cRew_cueBase(1,aucSt:aucEn));
        lastcRew(i,1:epocArrayLen) = cRew_cueBase(end,:);
        lastcRew_AMP(i,1) = max(cRew_cueBase(end,ampSt:ampEn));
        lastcRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cRew_cueBase(end,aucSt:aucEn));

        amp_cRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        auc_cRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
    elseif height(cRew_cueBase) == 1
        firstcRew(i,1:epocArrayLen) = cRew_cueBase(1,:);
        firstcRew_AMP(i,1) = max(cRew_cueBase(1,ampSt:ampEn));
        firstcRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cRew_cueBase(1,aucSt:aucEn));
        lastcRew(i,1:epocArrayLen) = nan;
        lastcRew_AMP(i,1) = nan;
        lastcRew_AUC(i,1) = nan;

        amp_cRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        auc_cRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
    elseif isempty(cRew_cueBase)
        firstcRew(i,1:epocArrayLen) = nan;
        firstcRew_AMP(i,1) = nan;
        firstcRew_AUC(i,1) = nan;
        lastcRew(i,1:epocArrayLen) = nan;
        lastcRew_AMP(i,1) = nan;
        lastcRew_AUC(i,1) = nan;

        amp_cRew_cueBase = nan;
        auc_cRew_cueBase = nan;
    end

    if ~isempty(cNoRew_cueBase) && height(cNoRew_cueBase) > 1
        firstcNoRew(i,1:epocArrayLen) = cNoRew_cueBase(1,:);
        firstcNoRew_AMP(i,1) = max(cNoRew_cueBase(1,ampSt:ampEn));
        firstcNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cNoRew_cueBase(1,aucSt:aucEn));
        lastcNoRew(i,1:epocArrayLen) = cNoRew_cueBase(end,:);
        lastcNoRew_AMP(i,1) = max(cNoRew_cueBase(end,ampSt:ampEn));
        lastcNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cNoRew_cueBase(end,aucSt:aucEn));

        amp_cNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 2,2));
        auc_cNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 2,2));
    elseif height(cNoRew_cueBase) == 1
        firstcNoRew(i,1:epocArrayLen) = cNoRew_cueBase(1,:);
        firstcNoRew_AMP(i,1) = max(cNoRew_cueBase(1,ampSt:ampEn));
        firstcNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),cNoRew_cueBase(1,aucSt:aucEn));
        lastcNoRew(i,1:epocArrayLen) = nan;
        lastcNoRew_AMP(i,1) = nan;
        lastcNoRew_AUC(i,1) = nan;

        amp_cNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 2,2));
        auc_cNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 2,2));
    elseif isempty(cNoRew_cueBase)
        firstcNoRew(i,1:epocArrayLen) = nan;
        firstcNoRew_AMP(i,1) = nan;
        firstcNoRew_AUC(i,1) = nan;
        lastcNoRew(i,1:epocArrayLen) = nan;
        lastcNoRew_AMP(i,1) = nan;
        lastcNoRew_AUC(i,1) = nan;

        amp_cNoRew_cueBase = nan;
        auc_cNoRew_cueBase = nan;
    end

    if ~isempty(iRew_cueBase) && height(iRew_cueBase) > 1
        firstiRew(i,1:epocArrayLen) = iRew_cueBase(1,:);
        firstiRew_AMP(i,1) = max(iRew_cueBase(1,ampSt:ampEn));
        firstiRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iRew_cueBase(1,aucSt:aucEn));
        lastiRew(i,1:epocArrayLen) = iRew_cueBase(end,:);
        lastiRew_AMP(i,1) = max(iRew_cueBase(end,ampSt:ampEn));
        lastiRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iRew_cueBase(end,aucSt:aucEn));

        amp_iRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 3,2));
        auc_iRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 3,2));
    elseif height(iRew_cueBase) == 1
        firstiRew(i,1:epocArrayLen) = iRew_cueBase(1,:);
        firstiRew_AMP(i,1) = max(iRew_cueBase(1,ampSt:ampEn));
        firstiRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iRew_cueBase(1,aucSt:aucEn));
        lastiRew(i,1:epocArrayLen) = nan;
        lastiRew_AMP(i,1) = nan;
        lastiRew_AUC(i,1) = nan;

        amp_iRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 3,2));
        auc_iRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 3,2));
    elseif isempty(iRew_cueBase)
        firstiRew(i,1:epocArrayLen) = nan;
        firstiRew_AMP(i,1) = nan;
        firstiRew_AUC(i,1) = nan;
        lastiRew(i,1:epocArrayLen) = nan;
        lastiRew_AMP(i,1) = nan;
        lastiRew_AUC(i,1) = nan;

        amp_iRew_cueBase = nan;
        auc_iRew_cueBase = nan;
    end

    if ~isempty(iNoRew_cueBase) && height(iNoRew_cueBase) > 1
        firstiNoRew(i,1:epocArrayLen) = iNoRew_cueBase(1,:);
        firstiNoRew_AMP(i,1) = max(iNoRew_cueBase(1,ampSt:ampEn));
        firstiNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iNoRew_cueBase(1,aucSt:aucEn));
        lastiNoRew(i,1:epocArrayLen) = iNoRew_cueBase(end,:);
        lastiNoRew_AMP(i,1) = max(iNoRew_cueBase(end,ampSt:ampEn));
        lastiNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iNoRew_cueBase(end,aucSt:aucEn));

        amp_iNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 4,2));
        auc_iNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 4,2));
    elseif height(iNoRew_cueBase) == 1
        firstiNoRew(i,1:epocArrayLen) = iNoRew_cueBase(1,:);
        firstiNoRew_AMP(i,1) = max(iNoRew_cueBase(1,ampSt:ampEn));
        firstiNoRew_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),iNoRew_cueBase(1,aucSt:aucEn));
        lastiNoRew(i,1:epocArrayLen) = nan;
        lastiNoRew_AMP(i,1) = nan;
        lastiNoRew_AUC(i,1) = nan;

        amp_iNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 4,2));
        auc_iNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 4,2));
    elseif isempty(iNoRew_cueBase)
        firstiNoRew(i,1:epocArrayLen) = nan;
        firstiNoRew_AMP(i,1) = nan;
        firstiNoRew_AUC(i,1) = nan;
        lastiNoRew(i,1:epocArrayLen) = nan;
        lastiNoRew_AMP(i,1) = nan;
        lastiNoRew_AUC(i,1) = nan;

        amp_iNoRew_cueBase = nan;
        auc_iNoRew_cueBase = nan;
    end

    leverLatencies(i,1) = mean(leverArray(:,1) - cueArray(:,1));
%     data.analysis.levers_z = levers_z;

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
                ampdFF(ii) = zeros(1);
                streams_z(ii,1:epocArrayLen) = NaN;
                ampZ(ii) = zeros(1);
                aucdFF(ii) = zeros(1);
                aucZ(ii) = zeros(1);
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
            % amplitude
            ampdFF(j) = max(streams_dFF(j,ampSt:ampEn));
            ampZ(j) = max(streams_z(j,ampSt:ampEn));
            % AUC
            aucdFF(j) = trapz(ts1(1,aucSt:aucEn),streams_dFF(j,aucSt:aucEn));
            aucZ(j) = trapz(ts1(1,aucSt:aucEn),streams_z(j,aucSt:aucEn));
        end
            if k == 1
                firstCue(i,1:epocArrayLen) = streams_z(1,:);
                firstCue_AMP(i,1) = max(streams_z(1,ampSt:ampEn));
                firstCue_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),streams_z(1,aucSt:aucEn));
                lastCue(i,1:epocArrayLen) = streams_z(end,:);
                lastCue_AMP(i,1) = max(streams_z(end,ampSt:ampEn));
                lastCue_AUC(i,1) = trapz(ts1(1,aucSt:aucEn),streams_z(end,aucSt:aucEn));
            
            end
            outputSTREAMSraw{k,1} = streams_raw(:,1:epocArrayLen);
            outputSTREAMSdFF{k,1} = streams_dFF(:,1:epocArrayLen);
            outputSTREAMSz{k,1} = streams_z(:,1:epocArrayLen);
            outputAMPdFF{k,1} = ampdFF;
            outputAMPz{k,1} = ampZ;
            outputAUCdFF{k,1} = aucdFF;
            outputAUCz{k,1} = aucZ;
        
    end
    



    master_cue_STREAMraw(i,:) = mean(outputSTREAMSraw{1,1},1);
    master_cRew_STREAMraw(i,:) = mean(outputSTREAMSraw{2,1},1);
    master_cNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{3,1},1);
    master_iRew_STREAMraw(i,:) = mean(outputSTREAMSraw{4,1},1);
    master_iNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{5,1},1);

    master_cue_STREAMdFF(i,:) = mean(outputSTREAMSdFF{1,1},1);
    master_cRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{2,1},1);
    master_cNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{3,1},1);
    master_iRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{4,1},1);
    master_iNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{5,1},1);

    master_cue_STREAMz(i,:) = mean(outputSTREAMSz{1,1},1);
    master_cRew_STREAMz(i,:) = mean(outputSTREAMSz{2,1},1);
    master_cNoRew_STREAMz(i,:) = mean(outputSTREAMSz{3,1},1);
    master_iRew_STREAMz(i,:) = mean(outputSTREAMSz{4,1},1);
    master_iNoRew_STREAMz(i,:) = mean(outputSTREAMSz{5,1},1);

    master_cRew_cueBase_STREAMz(i,:) = mean(cRew_cueBase,1);
    master_cNoRew_cueBase_STREAMz(i,:) = mean(cNoRew_cueBase,1);
    master_iRew_cueBase_STREAMz(i,:) = mean(iRew_cueBase,1);
    master_iNoRew_cueBase_STREAMz(i,:) = mean(iNoRew_cueBase,1);

    amp_cueBase_analysis(i,1:5) = {mean(outputAMPz{1},1) mean(amp_cRew_cueBase,1)...
        mean(amp_cNoRew_cueBase,1) mean(amp_iRew_cueBase,1) mean(amp_iNoRew_cueBase,1)};
    auc_cueBase_analysis(i,1:5) = {mean(outputAUCz{1},1) mean(auc_cRew_cueBase,1)...
        mean(auc_cNoRew_cueBase,1) mean(auc_iRew_cueBase,1) mean(auc_iNoRew_cueBase,1)};
   
    
    AMPdFF_analysis(i,1:5) = {mean(outputAMPdFF{1},1) mean(outputAMPdFF{2},1)...
        mean(outputAMPdFF{3},1) mean(outputAMPdFF{4},1) mean(outputAMPdFF{5},1)};
    AMPz_analysis(i,1:5) = {mean(outputAMPz{1},1) mean(outputAMPz{2},1)...
        mean(outputAMPz{3},1) mean(outputAMPz{4},1) mean(outputAMPz{5},1)};

    AUCdFF_analysis(i,1:5) = {mean(outputAUCdFF{1},1) mean(outputAUCdFF{2},1)...
        mean(outputAUCdFF{3},1) mean(outputAUCdFF{4},1) mean(outputAUCdFF{5},1)};
    AUCz_analysis(i,1:5) = {mean(outputAUCz{1},1) mean(outputAUCz{2},1)...
        mean(outputAUCz{3},1) mean(outputAUCz{4},1) mean(outputAUCz{5},1)};

    AMPz_fst_lst_analysis(i,1:10) = {mean(firstCue_AMP,1) mean(lastCue_AMP,1) mean(firstcRew_AMP,1)...
        mean(lastcRew_AMP,1) mean(firstcNoRew_AMP,1) mean(lastcNoRew_AMP,1) mean(firstiRew_AMP,1) mean(lastiRew_AMP,1)...
        mean(firstiNoRew_AMP,1) mean(lastiNoRew_AMP,1)};
    AUCz_fst_lst_analysis(i,1:10) = {mean(firstCue_AUC,1) mean(lastCue_AUC,1) mean(firstcRew_AUC,1)...
        mean(lastcRew_AUC,1) mean(firstcNoRew_AUC,1) mean(lastcNoRew_AUC,1) mean(firstiRew_AUC,1) mean(lastiRew_AUC,1)...
        mean(firstiNoRew_AUC,1) mean(lastiNoRew_AUC,1)};
    

    % save(filename,'data');
        
end
idList = cell2table(IDs','VariableNames',{'ID'});
phaseList = cell2table(phaseList','VariableNames',{'Phase'});
treatList = cell2table(treatList','VariableNames',{'Treatment'});


%% Amplitude Table %%
AMPz_analysis_table = cell2table(AMPz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AMPz_analysis_table = horzcat(idList,treatList,AMPz_analysis_table);
amp_cueBase_table = cell2table(amp_cueBase_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
amp_cueBase_table = horzcat(idList,treatList,phaseList,amp_cueBase_table);
%% Area Under Curve Table %%
AUCz_analysis_table = cell2table(AUCz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AUCz_analysis_table = horzcat(idList,treatList,AUCz_analysis_table);
auc_cueBase_table = cell2table(auc_cueBase_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
auc_cueBase_table = horzcat(idList,treatList,phaseList,auc_cueBase_table);
%% Epoc Stream Tables Baselined to Cue %%
a_cue_STREAMz = array2table(master_cue_STREAMz);
a_cue_STREAMz = horzcat(idList,phaseList,treatList,a_cue_STREAMz);
a_cRew_cueBase_STREAMz = array2table(master_cRew_cueBase_STREAMz);
a_cRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_cRew_cueBase_STREAMz);
a_cNoRew_cueBase_STREAMz = array2table(master_cNoRew_cueBase_STREAMz);
a_cNoRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_cNoRew_cueBase_STREAMz);
a_iRew_cueBase_STREAMz = array2table(master_iRew_cueBase_STREAMz);
a_iRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_iRew_cueBase_STREAMz);
a_iNoRew_cueBase_STREAMz = array2table(master_iNoRew_cueBase_STREAMz);
a_iNoRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_iNoRew_cueBase_STREAMz);

a_cue_STREAMz = sortrows(a_cue_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_cRew_cueBase_STREAMz = sortrows(a_cRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_cNoRew_cueBase_STREAMz = sortrows(a_cNoRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_iRew_cueBase_STREAMz = sortrows(a_iRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_iNoRew_cueBase_STREAMz = sortrows(a_iNoRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});


if toPlot == 1
    %% Plots for each epoc %%
    f1 = figure;
    subplot(5,1,1)
    plot(ts1,mean(master_cue_STREAMz,1,'omitnan'),'b')
    hold on
    cue_std = std(master_cue_STREAMz,'omitnan')/sqrt(height(master_cue_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cue_STREAMz,1,'omitnan') + cue_std, fliplr(mean(master_cue_STREAMz,1,'omitnan') - cue_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Cue')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,2)
    plot(ts1,mean(master_cRew_STREAMz,1,'omitnan'),'b')
    hold on
    cRew_std = std(master_cRew_STREAMz,'omitnan')/sqrt(height(master_cRew_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cRew_STREAMz,1,'omitnan') + cRew_std, fliplr(mean(master_cRew_STREAMz,1,'omitnan') - cRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Correct Reward')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,3)
    plot(ts1,mean(master_cNoRew_STREAMz,1,'omitnan'),'b')
    hold on
    cNoRew_std = std(master_cNoRew_STREAMz,'omitnan')/sqrt(height(master_cNoRew_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cNoRew_STREAMz,1,'omitnan') + cNoRew_std, fliplr(mean(master_cNoRew_STREAMz,1,'omitnan') - cNoRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Correct No Reward')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,4)
    plot(ts1,mean(master_iRew_STREAMz,1,'omitnan'),'b')
    hold on
    iRew_std = std(master_iRew_STREAMz,'omitnan')/sqrt(height(master_iRew_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_iRew_STREAMz,1,'omitnan') + iRew_std, fliplr(mean(master_iRew_STREAMz,1,'omitnan') - iRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Incorrect Reward')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,5)
    plot(ts1,mean(master_iNoRew_STREAMz,1,'omitnan'),'b')
    hold on
    iNoRew_std = std(master_iNoRew_STREAMz,'omitnan')/sqrt(height(master_iNoRew_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_iNoRew_STREAMz,1,'omitnan') + iNoRew_std, fliplr(mean(master_iNoRew_STREAMz,1,'omitnan') - iNoRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Incorrect No Reward')
    ylabel('GrabDA dF/F (z-Score)')
    xlabel('Time (s)')
    sgtitle(treatment)
    
    f2 = figure;
    subplot(5,1,1)
    plot(ts1,mean(master_cue_STREAMz,1,'omitnan'),'b')
    hold on
    cue_std = std(master_cue_STREAMz,'omitnan')/sqrt(height(master_cue_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cue_STREAMz,1,'omitnan') + cue_std, fliplr(mean(master_cue_STREAMz,1,'omitnan') - cue_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Cue')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,2)
    plot(ts1,mean(master_cRew_cueBase_STREAMz,1,'omitnan'),'b')
    hold on
    cRew_std = std(master_cRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_cRew_cueBase_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cRew_cueBase_STREAMz,1,'omitnan') + cRew_std, fliplr(mean(master_cRew_cueBase_STREAMz,1,'omitnan') - cRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Correct Reward (F0 = Cue)')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,3)
    plot(ts1,mean(master_cNoRew_cueBase_STREAMz,1,'omitnan'),'b')
    hold on
    cNoRew_std = std(master_cNoRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_cNoRew_cueBase_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_cNoRew_cueBase_STREAMz,1,'omitnan') + cNoRew_std, fliplr(mean(master_cNoRew_cueBase_STREAMz,1,'omitnan') - cNoRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Correct No Reward (F0 = Cue)')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,4)
    plot(ts1,mean(master_iRew_cueBase_STREAMz,1,'omitnan'),'b')
    hold on
    iRew_std = std(master_iRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_iRew_cueBase_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_iRew_cueBase_STREAMz,1,'omitnan') + iRew_std, fliplr(mean(master_iRew_cueBase_STREAMz,1,'omitnan') - iRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Incorrect Reward (F0 = Cue)')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(5,1,5)
    plot(ts1,mean(master_iNoRew_cueBase_STREAMz,1,'omitnan'),'b')
    hold on
    iNoRew_std = std(master_iNoRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_iNoRew_cueBase_STREAMz));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(master_iNoRew_cueBase_STREAMz,1,'omitnan') + iNoRew_std, fliplr(mean(master_iNoRew_cueBase_STREAMz,1,'omitnan') - iNoRew_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Incorrect No Reward (F0 = Cue)')
    ylabel('GrabDA dF/F (z-Score)')
    xlabel('Time (s)')
    sgtitle(treatment)
    
    f3 = figure;
    subplot(1,2,1)
    plot(ts1,mean(firstLever,1,'omitnan'),'b')
    hold on
    Fst_std = std(firstLever,'omitnan')/sqrt(height(firstLever));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(firstLever,1,'omitnan') + Fst_std, fliplr(mean(firstLever,1,'omitnan') - Fst_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('First Lever Press')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(1,2,2)
    plot(ts1,mean(lastLever,1,'omitnan'),'b')
    hold on
    Lst_std = std(lastLever,'omitnan')/sqrt(height(lastLever));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(lastLever,1,'omitnan') + Lst_std, fliplr(mean(lastLever,1,'omitnan') - Lst_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Last Lever Press')
    ylabel('GrabDA dF/F (z-Score)')
    sgtitle(treatment)
    
    f4 = figure;
    subplot(1,2,1)
    plot(ts1,mean(firstCue,1,'omitnan'),'b')
    hold on
    Fst_std = std(firstCue,'omitnan')/sqrt(height(firstCue));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(firstCue,1,'omitnan') + Fst_std, fliplr(mean(firstCue,1,'omitnan') - Fst_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('First Cue')
    ylabel('GrabDA dF/F (z-Score)')
    
    subplot(1,2,2)
    plot(ts1,mean(lastCue,1,'omitnan'),'b')
    hold on
    Lst_std = std(lastCue,'omitnan')/sqrt(height(lastCue));
    XX = [ts1, fliplr(ts1)];
    YY = [mean(lastCue,1,'omitnan') + Lst_std, fliplr(mean(lastCue,1,'omitnan') - Lst_std)];
    % Plot filled standard error bands.
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    % Plot line at x=0
    xline(0,'LineWidth',2,'Color','black')
    title('Last Cue')
    ylabel('GrabDA dF/F (z-Score)')
    sgtitle(treatment)
end
amp_cueBase_table = sortrows(amp_cueBase_table,{'Phase','Treatment'},{'ascend','descend'});
auc_cueBase_table = sortrows(auc_cueBase_table,{'Phase','Treatment'},{'ascend','descend'});

prl_stream_analysis.acq1.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Acq1'),:);
prl_stream_analysis.acq2.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Acq2'),:);
prl_stream_analysis.rev1.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev1'),:);
prl_stream_analysis.rev2.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev2'),:);
prl_stream_analysis.rev3.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev3'),:);
prl_stream_analysis.acq1.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Acq1'),:);
prl_stream_analysis.acq2.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Acq2'),:);
prl_stream_analysis.rev1.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev1'),:);
prl_stream_analysis.rev2.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev2'),:);
prl_stream_analysis.rev3.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev3'),:);

prl_stream_analysis.acq1.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Acq1'),:);

prl_stream_analysis.acq2.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Acq2'),:);

prl_stream_analysis.rev1.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev1'),:);

prl_stream_analysis.rev2.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev2'),:);

prl_stream_analysis.rev3.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev3'),:);

prl_stream_analysis.info.time = ts1;

save('../data-files/prl_stream_analysis.mat','prl_stream_analysis')
toc

disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);