clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ_cue = [5 1];
baselineZ_lever = [3 1];
N = 1; %Downsample N times
minArrayLen = 7121; 
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('Z:\DA_PRL\PRL_Mats_Split\Good_Sig','Choose the .mat files you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
try
    tic
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
    myFiles = myFiles(endsWith({myFiles.name},'.mat'));
    numFiles = length(myFiles);
    LOAD_BAR = waitbar(0,'1','Name','Analyzing...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(LOAD_BAR,'canceling',0)
    AUC_analysis = cell(numFiles,5);
    AMP_analysis = cell(numFiles,7);
    master_cue_STREAM = zeros(numFiles,minArrayLen);
    master_cRew_STREAM = zeros(numFiles,minArrayLen);
    master_cNoRew_STREAM = zeros(numFiles,minArrayLen);
    master_iRew_STREAM = zeros(numFiles,minArrayLen);
    master_iNoRew_STREAM = zeros(numFiles,minArrayLen);
    master_winStay_STREAM = zeros(numFiles,minArrayLen);
    master_loseShift_STREAM = zeros(numFiles,minArrayLen);
    for i = 1:numFiles
        filename = fullfile(myDir,myFiles(i).name);
        [~,name,~] = fileparts(filename);
        [~,treatment,~] = fileparts(myDir);
        brokenID = strsplit(name,'_');
        prl_phase = char(brokenID(2));
        load(filename)
        sessionSTREAM = [];
        cueSTREAM = [];
        cRewSTREAM = [];
        cNoRewSTREAM = [];
        iRewSTREAM = [];
        iNoRewSTREAM = [];
        if isfield(data.streams, 'x405A')
            ISOS = 'x405A';
            GRABDA = 'x465A';
            cueTS = data.epocs.St1_.onset;
            correct_rewarded = data.epocs.cRewA.onset;
            correct_noreward = data.epocs.cNoRewA.onset;
            incorrect_rewarded = data.epocs.iRewA.onset;
            incorrect_noreward = data.epocs.iNoRewA.onset;
        elseif isfield(data.streams, 'x405C')
            ISOS = 'x405C';
            GRABDA = 'x465C';
            cueTS = data.epocs.St2_.onset;
            correct_rewarded = data.epocs.cRewC.onset;
            correct_noreward = data.epocs.cNoRewC.onset;
            incorrect_rewarded = data.epocs.iRewC.onset;
            incorrect_noreward = data.epocs.iNoRewC.onset;
        else
            disp('No streams detected')
            
        end
        %time array used for all streams%
        time = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
        %removes the first (t) seconds where the data is wild due to turning on LEDs%
        t = 0; % time threshold below which we will discard
        ind = find(time>t,1);% find first index of when time crosses threshold
        time = time(ind:end); % reformat vector to only include allowed time
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(ind:end);
        data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
        
        %downsample streams and time array by N times%
        data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
        data.streams.(GRABDA).data = downsample(data.streams.(GRABDA).data, N);
        minLength = min(length(data.streams.(ISOS).data),length(data.streams.(GRABDA).data));
        data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(1:minLength);
        time = downsample(time, N);
        ts1 = -baseline + (1:minArrayLen) / data.streams.(GRABDA).fs*N;
        
        %detrend & dFF%
        bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
        Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
        Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
        dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(dFF));
        detrend_465 = detrend(dFF);
        z465 = zscore(detrend_465);
    
        [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cueTS,correct_rewarded,...
            correct_noreward,incorrect_rewarded,incorrect_noreward);
        winStayTS = errorProbExtract(trial_type,session_ts,1,1);
        winShiftTS = errorProbExtract(trial_type,session_ts,2,1);
        loseStayTS = errorProbExtract(trial_type,session_ts,3,1);
        loseShiftTS = errorProbExtract(trial_type,session_ts,4,1);
    
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            session_ts, ...
            baseline, ...
            timeWindow, ...
            baselineZ_cue, ...
            [0,5], ...
            minArrayLen ...
            );
        sessionSTREAM = epocSTREAM;
        sessionAMP = epocAMP;
        sessionAUC = epocAUC;
    
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            cueTS, ...
            baseline, ...
            timeWindow, ...
            baselineZ_cue, ...
            [0,2], ...
            minArrayLen ...
            );
        cueSTREAM = epocSTREAM;
        cueAMP = epocAMP;
        cueAUC = epocAUC;
        
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            correct_rewarded, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [0,2], ...
            minArrayLen ...
            );
        cRewSTREAM = epocSTREAM;
        cRewAMP = epocAMP;
        cRewAUC = epocAUC;
        
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            correct_noreward, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [2,4], ...
            minArrayLen ...
            );
        cNoRewSTREAM = epocSTREAM;
        cNoRewAMP = epocAMP;
        cNoRewAUC = epocAUC;
        
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            incorrect_rewarded, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [0,2], ...
            minArrayLen ...
            );
        iRewSTREAM = epocSTREAM;
        iRewAMP = epocAMP;
        iRewAUC = epocAUC;
     
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            incorrect_noreward, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [2,4], ...
            minArrayLen ...
            );
        iNoRewSTREAM = epocSTREAM;
        iNoRewAMP = epocAMP;
        iNoRewAUC = epocAUC;
    
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            winStayTS, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [0,2], ...
            minArrayLen ...
            );
        winStaySTREAM = epocSTREAM;
        winStayAMP = epocAMP;
        winStayAUC = epocAUC;
    
        [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
            detrend_465, ...
            time, ...
            loseShiftTS, ...
            baseline, ...
            timeWindow, ...
            baselineZ_lever, ...
            [0,2], ...
            minArrayLen ...
            );
        loseShiftSTREAM = epocSTREAM;
        loseShiftAMP = epocAMP;
        loseShiftAUC = epocAUC;
    
        
        master_cue_STREAM(i,:) = mean(cueSTREAM,1);
        master_cRew_STREAM(i,:) = mean(cRewSTREAM,1);
        master_cNoRew_STREAM(i,:) = mean(cNoRewSTREAM,1);
        master_iRew_STREAM(i,:) = mean(iRewSTREAM,1);
        master_iNoRew_STREAM(i,:) = mean(iNoRewSTREAM,1);
        if ~isempty(winStaySTREAM)
            master_winStay_STREAM(i,:) = mean(winStaySTREAM,1);
        else
            continue
        end
        if ~isempty(loseShiftSTREAM)
            master_loseShift_STREAM(i,:) = mean(loseShiftSTREAM,1);
        else
            continue
        end
        cueAMP(cueAMP == 0) = nan;
        cRewAMP(cRewAMP == 0) = nan;
        cNoRewAMP(cNoRewAMP == 0) = nan;
        iRewAMP(iRewAMP == 0) = nan;
        iNoRewAMP(iNoRewAMP == 0) = nan;
        winStayAMP(winStayAMP == 0) = nan;
        loseShiftAMP(loseShiftAMP == 0) = nan;
        AUC_analysis(i,1:5) = {mean(cueAUC,1) mean(cRewAUC,1) mean(cNoRewAUC,1) mean(iRewAUC,1) mean(iNoRewAUC,1)};
    
        AMP_analysis(i,1:5) = {mean(cueAMP,1) mean(cRewAMP,1) mean(cNoRewAMP,1) mean(iRewAMP,1) mean(iNoRewAMP,1)};
    
        
        waitbar(i/numFiles,LOAD_BAR,sprintf('Progress: %d %%',floor(i/numFiles*100)));
        pause(0.1)        
    end
    AUC_analysis(cellfun(@(x) x==0,AUC_analysis)) = {NaN};
    AMP_analysis(cellfun(@(x) x==0,AMP_analysis)) = {NaN};
    AUC_analysis_table = cell2table(AUC_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
    AMP_analysis_table = cell2table(AMP_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
    % master_cue_STREAM = master_cue_STREAM;
    % master_cRew_STREAM = master_cRew_STREAM;
    % master_cNoRew_STREAM = master_cNoRew_STREAM;
    % master_iRew_STREAM = master_iRew_STREAM;
    % master_iNoRew_STREAM = master_iNoRew_STREAM;
    
    toc
    delete(LOAD_BAR)
    disp("Successfully analyzed .mat files")
    fprintf("Files analyzed: %d\n", numFiles)
    
    NERD_STATS(toc,numFiles);
catch ME
    delete(LOAD_BAR)
    rethrow(ME)
end