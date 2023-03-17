clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ = [5 1];
N = 100; %Downsample N times
minArrayLen = 72; 
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('Z:\DA_PRL','Choose the .mat files you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
LOAD_BAR = waitbar(0,'1','Name','Analyzing...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(LOAD_BAR,'canceling',0)
AUC_analysis = cell(numFiles,5);
AMP_analysis = cell(numFiles,5);
master_cue_STREAM = zeros(numFiles,minArrayLen);
master_cRew_STREAM = zeros(numFiles,minArrayLen);
master_cNoRew_STREAM = zeros(numFiles,minArrayLen);
master_iRew_STREAM = zeros(numFiles,minArrayLen);
master_iNoRew_STREAM = zeros(numFiles,minArrayLen);
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    load(filename)
    cueSTREAM = [];
    cRewSTREAM = [];
    cNoRewSTREAM = [];
    iRewSTREAM = [];
    iNoRewSTREAM = [];
    if isfield(data.streams, 'x405A')
        ISOS = 'x405A';
        GRABDA = 'x465A';
        %time array used for all streams%
        time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
        %removes the first (t) seconds where the data is wild due to turning on LEDs%
        t = 0; % time threshold below which we will discard
        ind = find(time1>t,1);% find first index of when time crosses threshold
        time1 = time1(ind:end); % reformat vector to only include allowed time
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(ind:end);
        data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
        
        %downsample streams and time array by N times%
        data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
        data.streams.(GRABDA).data = downsample(data.streams.(GRABDA).data, N);
        minLength = min(length(data.streams.(ISOS).data),length(data.streams.(GRABDA).data));
        data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(1:minLength);
        time1 = downsample(time1, N);
        ts1 = -baseline + (1:minArrayLen) / data.streams.(GRABDA).fs*N;
        
        %detrend & dFF%
        bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
        Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
        Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
        dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(dFF));
        detrend_465 = detrend(dFF);
        z465 = zscore(data.streams.(GRABDA).data);
        cueTSA = data.epocs.St1_.onset;
        correct_rewardedA = data.epocs.cRewA.onset;
        correct_norewardA = data.epocs.cNoRewA.onset;
        incorrect_rewardedA = data.epocs.iRewA.onset;
        incorrect_norewardA = data.epocs.iNoRewA.onset;

        session_ts = sessionArraySort(correct_rewardedA,...
            correct_norewardA,incorrect_rewardedA,incorrect_norewardA);

        for e = 1:height(session_ts)
            e1 = session_ts(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = session_ts(e,1)-baselineZ(:,2);
            e4 = session_ts(e,1)-baselineZ(:,1);
            [~,ind1] = min(abs(time1-e1));
            [~,ind2] = min(abs(time1-e2));
            [~,ind3] = min(abs(time1-e3));
            [~,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                break
            end
            sessionSTREAM(e,:) = zfinal(1:minArrayLen);
            sessionAMP(e,:) = max(zfinal);
            sessionAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        for e = 1:height(cueTSA)
            e1 = data.epocs.St1_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St1_.onset(e,1)-baselineZ(:,2);
            e4 = data.epocs.St1_.onset(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                break
            end
            cueSTREAM(e,:) = zfinal(1:minArrayLen);
            cueAMP(e,:) = max(zfinal);
            cueAUC(e,:) = trapz(GRABDA_time,zfinal);
        end

        for e = 1:height(correct_rewardedA)
            e1 = correct_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedA(e,1)-baselineZ(:,2);
            e4 = correct_rewardedA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if correct_rewardedA == 0
                cRewSTREAM(1:minArrayLen,1) = zeros;
                cRewAMP = zeros(1);
                cRewAUC = zeros(1);
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cRewSTREAM(e,:) = zfinal(1:minArrayLen);
            cRewAMP(e,:) = max(zfinal);
            cRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end

        for e = 1:height(correct_norewardA)
            e1 = correct_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardA(e,1)-baselineZ(:,2);
            e4 = correct_norewardA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if correct_norewardA == 0
                cNoRewSTREAM(1:minArrayLen,1) = zeros;
                cNoRewAMP = zeros(1);
                cNoRewAUC = zeros(1);
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            cNoRewAMP(e,:) = max(zfinal);
            cNoRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end

        for e = 1:height(incorrect_rewardedA)
            e1 = incorrect_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedA(e,1)-baselineZ(:,2);
            e4 = incorrect_rewardedA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if incorrect_rewardedA == 0
                iRewSTREAM(1:minArrayLen,1) = zeros;
                iRewAMP = zeros(1);
                iRewAUC = zeros(1);
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            iRewSTREAM(e,:) = zfinal(1:minArrayLen);
            iRewAMP(e,:) = max(zfinal);
            iRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end

        for e = 1:height(incorrect_norewardA)
            e1 = incorrect_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardA(e,1)-baselineZ(:,2);
            e4 = incorrect_norewardA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            iNoRewAMP(e,:) = max(zfinal);
            iNoRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end

    else
        ISOS = 'x405C';
        GRABDA = 'x465C';
        %time array used for all streams%
        time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
        %removes the first (t) seconds where the data is wild due to turning on LEDs%
        t = 0; % time threshold below which we will discard
        ind = find(time1>t,1);% find first index of when time crosses threshold
        time1 = time1(ind:end); % reformat vector to only include allowed time
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(ind:end);
        data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
        minLength = min(length(data.streams.(ISOS).data),length(data.streams.(GRABDA).data));
        data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
        data.streams.(GRABDA).data = data.streams.(GRABDA).data(1:minLength);
        
        %downsample streams and time array by N times%
        data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
        data.streams.(GRABDA).data = downsample(data.streams.(GRABDA).data, N);
        time1 = downsample(time1, N);
        ts1 = -baseline + (1:minArrayLen) / data.streams.(GRABDA).fs*N;
        
        %detrend & dFF%
        bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
        Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
        Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
        dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(dFF));
        detrend_465 = detrend(dFF);
        z465 = zscore(data.streams.(GRABDA).data);
        cueTSC = data.epocs.St2_.onset;
        correct_rewardedC = data.epocs.cRewC.onset;
        correct_norewardC = data.epocs.cNoRewC.onset;
        incorrect_rewardedC = data.epocs.iRewC.onset;
        incorrect_norewardC = data.epocs.iNoRewC.onset;

        session_ts = sessionArraySort(correct_rewardedC,...
            correct_norewardC,incorrect_rewardedC,incorrect_norewardC);

        for e = 1:height(session_ts)
            e1 = session_ts(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = session_ts(e,1)-baselineZ(:,2);
            e4 = session_ts(e,1)-baselineZ(:,1);
            [~,ind1] = min(abs(time1-e1));
            [~,ind2] = min(abs(time1-e2));
            [~,ind3] = min(abs(time1-e3));
            [~,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                break
            end
            sessionSTREAM(e,:) = zfinal(1:minArrayLen);
            sessionAMP(e,:) = max(zfinal);
            sessionAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        for e = 1:height(cueTSC)
            e1 = data.epocs.St2_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St2_.onset(e,1)-baselineZ(:,2);
            e4 = data.epocs.St2_.onset(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                break
            end
            cueSTREAM(e,:) = zfinal(1:minArrayLen);
            cueAMP(e,:) = max(zfinal);
            cueAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        
        for e = 1:height(correct_rewardedC)
            e1 = correct_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedC(e,1)-baselineZ(:,2);
            e4 = correct_rewardedC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if correct_rewardedA == 0
                cRewSTREAM(1:minArrayLen,1) = zeros;
                cRewAMP = zeros(1);
                cRewAUC = zeros(1);
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cRewSTREAM(e,:) = zfinal(1:minArrayLen);
            cRewAMP(e,:) = max(zfinal);
            cRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        
        for e = 1:height(correct_norewardC)
            e1 = correct_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardC(e,1)-baselineZ(:,2);
            e4 = correct_norewardC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if correct_norewardA == 0
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            cNoRewAMP(e,:) = max(zfinal);
            cNoRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        for e = 1:height(incorrect_rewardedC)
            e1 = incorrect_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedC(e,1)-baselineZ(:,2);
            e4 = incorrect_rewardedC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if incorrect_rewardedC == 0
                iRewSTREAM(1:minArrayLen,1) = zeros;
                iRewAMP = zeros(1);
                iRewAUC = zeros(1);
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            iRewSTREAM(e,:) = zfinal(1:minArrayLen);
            iRewAMP(e,:) = max(zfinal);
            iRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        
        for e = 1:height(incorrect_norewardC)
            e1 = incorrect_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardC(e,1)-baselineZ(:,2);
            e4 = incorrect_norewardC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if incorrect_norewardC == 0
                iNoRewSTREAM(1:minArrayLen,1) = zeros;
                iNoRewAMP = zeros(1);
                iNoRewAUC = zeros(1);
                break
            end
            iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            iNoRewAMP(e,:) = max(zfinal);
            iNoRewAUC(e,:) = trapz(GRABDA_time,zfinal);
        end
        
    end
    master_cue_STREAM(i,:) = mean(cueSTREAM);
    master_cRew_STREAM(i,:) = mean(cRewSTREAM);
    master_cNoRew_STREAM(i,:) = mean(cNoRewSTREAM);
    master_iRew_STREAM(i,:) = mean(iRewSTREAM);
    master_iNoRew_STREAM(i,:) = mean(iNoRewSTREAM);

    AUC_analysis(i,1:5) = {mean(cueAUC) mean(cRewAUC) mean(cNoRewAUC) mean(iRewAUC) mean(iNoRewAUC)};

    AMP_analysis(i,1:5) = {mean(cueAMP) mean(cRewAMP) mean(cNoRewAMP) mean(iRewAMP) mean(iNoRewAMP)};

    
    waitbar(i/numFiles,LOAD_BAR,sprintf('Progress: %d %%',floor(i/numFiles*100)));
    pause(0.1)        
end
AUC_analysis_table = cell2table(AUC_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
master_cue_STREAM = master_cue_STREAM.';
master_cRew_STREAM = master_cRew_STREAM.';
master_cNoRew_STREAM = master_cNoRew_STREAM.';
master_iRew_STREAM = master_iRew_STREAM.';
master_iNoRew_STREAM = master_iNoRew_STREAM.';

toc
delete(LOAD_BAR)
disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);