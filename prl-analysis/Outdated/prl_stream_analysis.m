clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ_cue = [3 1];
amp_window = [-1 2];
baselineZ_lever = [3 1];
N = 10; %Downsample N times
minArrayLen = 713; 
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample) 
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
AUCz_analysis = cell(numFiles,5);
AMPz_analysis = cell(numFiles,5);
master_cue_STREAMz = zeros(numFiles,minArrayLen);
master_cRew_STREAMz = zeros(numFiles,minArrayLen);
master_cNoRew_STREAMz = zeros(numFiles,minArrayLen);
master_iRew_STREAMz = zeros(numFiles,minArrayLen);
master_iNoRew_STREAMz = zeros(numFiles,minArrayLen);
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    [~,treatment,~] = fileparts(myDir);
    brokenID = strsplit(name,'_');
    prl_phase = char(brokenID(2));
    load(filename)
    sessionSTREAMdff = [];
    cueSTREAMdff = [];
    cRewSTREAMdff = [];
    cNoRewSTREAMdff = [];
    iRewSTREAMdff = [];
    iNoRewSTREAMdff = [];
    sessionSTREAMz = [];
    cueSTREAMz = [];
    cRewSTREAMz = [];
    cNoRewSTREAMz = [];
    iRewSTREAMz = [];
    iNoRewSTREAMz = [];
    if isfield(data.streams, 'x405A')
        ISOS = 'x405A';
        SIGNAL = 'x465A';
        cue = data.epocs.St1_.onset;
        cRew = data.epocs.cRewA.onset;
        cNoRew = data.epocs.cNoRewA.onset;
        iRew = data.epocs.iRewA.onset;
        iNoRew = data.epocs.iNoRewA.onset;
        [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
    elseif isfield(data.streams, 'x405C')
        ISOS = 'x405C';
        SIGNAL = 'x465C';
        cue = data.epocs.St2_.onset;
        cRew = data.epocs.cRewC.onset;
        cNoRew = data.epocs.cNoRewC.onset;
        iRew = data.epocs.iRewC.onset;
        iNoRew = data.epocs.iNoRewC.onset;
        [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
    end
    %time array used for all streams%
    session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    t = 0; % time threshold below which we will discard
    ind = find(session_time>t,1);% find first index of when time crosses threshold
    session_time = session_time(ind:end); % reformat vector to only include allowed time
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
   
    %downsample streams and time array by N times%
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
    minLength = min(length(data.streams.(ISOS).data),length(data.streams.(SIGNAL).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(1:minLength);

    ISOS_raw = data.streams.(ISOS).data;
    SIGNAL_raw = data.streams.(SIGNAL).data;

    session_time = downsample(session_time, N);
    ts1 = -baseline + (1:minArrayLen) / data.streams.(SIGNAL).fs*N;
    
    %detrend & dFF%
    bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,1);
    Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
    Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    detrend_465 = detrend(dFF);
    z465 = zscore(detrend_465);



    for e = 1:height(session_ts)
        e1 = session_ts(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        e3 = session_ts(e,1)-baselineZ_cue(:,2);
        e4 = session_ts(e,1)-baselineZ_cue(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        sessionSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        sessionAMPdff(e,:) = max(GRABDA_dffSignal);
        sessionSTREAMz(e,:) = zfinal(1:minArrayLen);
        sessionAMPz(e,:) = max(zfinal);
        
    end
    for e = 1:height(cue)
        e1 = data.epocs.St1_.onset(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        e3 = data.epocs.St1_.onset(e,1)-baselineZ_cue(:,2);
        e4 = data.epocs.St1_.onset(e,1)-baselineZ_cue(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        ind5 = find(ts1>amp_window(1),1);
        ind6 = find(ts1>amp_window(2),1);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        cueSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        cueAMPdff(e,:) = max(GRABDA_dffSignal(1, ind5:ind6));
        cueSTREAMz(e,:) = zfinal(1:minArrayLen);
        cueAMPz(e,:) = max(zfinal(1, ind5:ind6));
        
    end

    for e = 1:height(cRew)
        e1 = cRew(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        if cRew == 0
            cRewSTREAMdff(e,1:minArrayLen) = zeros;
            cRewAMPdff = zeros(1);
            cRewSTREAMz(e,1:minArrayLen) = zeros;
            cRewAMPz = zeros(1);
            break
        end
        e3 = cRew(e,1)-baselineZ_lever(:,2);
        e4 = cRew(e,1)-baselineZ_lever(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        ind5 = find(ts1>amp_window(1),1);
        ind6 = find(ts1>amp_window(2),1);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        cRewSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        cRewAMPdff(e,:) = max(GRABDA_dffSignal(1, ind5:ind6));
        cRewSTREAMz(e,:) = zfinal(1:minArrayLen);
        cRewAMPz(e,:) = max(zfinal(1, ind5:ind6));
        
    end

    for e = 1:height(cNoRew)
        e1 = cNoRew(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        if cNoRew == 0
            cNoRewSTREAMdff(e,1:minArrayLen) = zeros;
            cNoRewAMPdff = zeros(1);
            cNoRewSTREAMz(e,1:minArrayLen) = zeros;
            cNoRewAMPz = zeros(1);
            break
        end
        e3 = cNoRew(e,1)-baselineZ_lever(:,2);
        e4 = cNoRew(e,1)-baselineZ_lever(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        ind5 = find(ts1>amp_window(1),1);
        ind6 = find(ts1>amp_window(2),1);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        cNoRewSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        cNoRewAMPdff(e,:) = max(GRABDA_dffSignal(1, ind5:ind6));
        cNoRewSTREAMz(e,:) = zfinal(1:minArrayLen);
        cNoRewAMPz(e,:) = max(zfinal(1, ind5:ind6));
    end

    for e = 1:height(iRew)
        e1 = iRew(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        if iRew == 0
            iRewSTREAMdff(e,1:minArrayLen) = zeros;
            iRewAMPdff = zeros(1);
            iRewSTREAMz(e,1:minArrayLen) = zeros;
            iRewAMPz = zeros(1);
            break
        end
        e3 = iRew(e,1)-baselineZ_lever(:,2);
        e4 = iRew(e,1)-baselineZ_lever(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        ind5 = find(ts1>amp_window(1),1);
        ind6 = find(ts1>amp_window(2),1);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        iRewSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        iRewAMPdff(e,:) = max(GRABDA_dffSignal(1, ind5:ind6));
        iRewSTREAMz(e,:) = zfinal(1:minArrayLen);
        iRewAMPz(e,:) = max(zfinal(1, ind5:ind6));
    end

    for e = 1:height(iNoRew)
        e1 = iNoRew(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        if iNoRew == 0
            iNoRewSTREAMdff(e,1:minArrayLen) = zeros;
            iNoRewAMPdff = zeros(1);
            iNoRewSTREAMz(e,1:minArrayLen) = zeros;
            iNoRewAMPz = zeros(1);
            break
        end
        e3 = iNoRew(e,1)-baselineZ_lever(:,2);
        e4 = iNoRew(e,1)-baselineZ_lever(:,1);
        [~,ind1] = min(abs(session_time-e1));
        [~,ind2] = min(abs(session_time-e2));
        [~,ind3] = min(abs(session_time-e3));
        [~,ind4] = min(abs(session_time-e4));
        GRABDA_dffBase = mean(detrend_465(1,ind4:ind3));
        GRABDA_dffSignal = detrend_465(1,ind1:ind2) - GRABDA_dffBase;
        GRABDA_zBase = z465(1,ind4:ind3);
        GRABDA_zSignal = z465(1,ind1:ind2);
        GRABDA_time = session_time(1,ind1:ind2);
        ind5 = find(ts1>amp_window(1),1);
        ind6 = find(ts1>amp_window(2),1);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_zSignal - zb)/zsd;
        if length(zfinal) < minArrayLen || length(GRABDA_dffSignal) < minArrayLen
            continue
        end
        iNoRewSTREAMdff(e,:) = GRABDA_dffSignal(1:minArrayLen);
        iNoRewAMPdff(e,:) = max(GRABDA_dffSignal(1, ind5:ind6));
        iNoRewSTREAMz(e,:) = zfinal(1:minArrayLen);
        iNoRewAMPz(e,:) = max(zfinal(1, ind5:ind6));
    end

    

    master_cue_STREAMz(i,:) = mean(cueSTREAMz,1);
    master_cRew_STREAMz(i,:) = mean(cRewSTREAMz,1);
    master_cNoRew_STREAMz(i,:) = mean(cNoRewSTREAMz,1);
    master_iRew_STREAMz(i,:) = mean(iRewSTREAMz,1);
    master_iNoRew_STREAMz(i,:) = mean(iNoRewSTREAMz,1);
    cueAMPz(cueAMPz == 0) = nan;
    cRewAMPz(cRewAMPz == 0) = nan;
    cNoRewAMPz(cNoRewAMPz == 0) = nan;
    iRewAMPz(iRewAMPz == 0) = nan;
    iNoRewAMPz(iNoRewAMPz == 0) = nan;
    AMPz_analysis(i,1:5) = {mean(cueAMPz,1) mean(cRewAMPz,1) mean(cNoRewAMPz,1) mean(iRewAMPz,1) mean(iNoRewAMPz,1)};

    
        
end
AMPdff_analysis(cellfun(@(x) x==0,AMPdff_analysis)) = {NaN};
AMPdff_analysis_table = cell2table(AMPdff_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AMPz_analysis(cellfun(@(x) x==0,AMPz_analysis)) = {NaN};
AMPz_analysis_table = cell2table(AMPz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});

toc

disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);