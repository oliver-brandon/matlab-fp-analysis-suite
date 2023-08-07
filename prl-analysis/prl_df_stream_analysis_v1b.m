clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0 timeWindow]; % time window to grab amplitude from
auc_window = [-1 timeWindow];
t = 10; % seconds to clip from start of streams
N = 1; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
toPlot = 0; % 1 = plot figures, 0 = don't plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/My Drive/prl/dual_fiber/tanks','Choose a folder containing the tank(s) you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a folder with tank(s) to start")
    return
end
tic

myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(~endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
IDs = {};
phaseList = {};
treatList = {};

for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    data = TDTbin2mat(filename,'TYPE',{'streams','epocs'});
    [~,name,~] = fileparts(filename);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    tempID = char(brokenID{1});
    if strcmp(tempID,emptyID)
        IDs{i} = cellstr(strtrim(brokenID{4}));
        phaseList{i} = cellstr(strtrim(brokenID{5}));
        treatList{i} = cellstr(strtrim(brokenID{6}));
        TTLs = 2;
    elseif ~strcmp(tempID,emptyID)
        IDs{i} = cellstr(strtrim(brokenID{1}));
        phaseList{i} = cellstr(strtrim(brokenID{2}));
        treatList{i} = cellstr(strtrim(brokenID{3}));
        TTLs = 1;
    end
    data = prl_df_epocs(data,TTLs);
    if TTLs == 1
        cue = data.epocs.St1_.onset;
        cRew = data.epocs.cRewA.onset;
        cNoRew = data.epocs.cNoRewA.onset;
        iRew = data.epocs.iRewA.onset;
        iNoRew = data.epocs.iNoRewA.onset;
    elseif TTLs == 2
        cue = data.epocs.St2_.onset;
        cRew = data.epocs.cRewC.onset;
        cNoRew = data.epocs.cNoRewC.onset;
        iRew = data.epocs.iRewC.onset;
        iNoRew = data.epocs.iNoRewC.onset;
    end
    
    [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(1:end-1,:);
    end
    [errorProbLeverTS, errorProbCueTS] = errorProbExtract(session_identifiers, 1, 1);
    cueArray = session_identifiers(1:2:end,:);
    leverArray = session_identifiers(2:2:end,:);
    DLS_levers_z = zeros(height(leverArray), epocArrayLen);
    DLS_levers_raw = zeros(height(leverArray), epocArrayLen);
    NAc_levers_z = zeros(height(leverArray), epocArrayLen);
    NAc_levers_raw = zeros(height(leverArray), epocArrayLen);
    for k = 1:2
        if k == 1
            ISOS = 'x405A';
            SIGNAL = 'x465A';
            ISOS_raw = data.streams.(ISOS).data;
            SIGNAL_raw = data.streams.(SIGNAL).data;
            %time array used for all streams%
            session_time = (1:length(SIGNAL_raw))/data.streams.(SIGNAL).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            ind = find(session_time>t,1);% find first index of when time crosses threshold
            session_time = session_time(ind:end); % reformat vector to only include allowed time
            SIGNAL_raw = data.streams.(SIGNAL).data(ind:end);
            ISOS_raw = data.streams.(ISOS).data(ind:end);
            
            %downsample streams and time array by N times%
            ISOS_raw = downsample(ISOS_raw, N);
            SIGNAL_raw = downsample(SIGNAL_raw, N);
            minLength = min(length(ISOS_raw),length(SIGNAL_raw));
            ISOS_raw = ISOS_raw(1:minLength);
            SIGNAL_raw = SIGNAL_raw(1:minLength);
            session_time = downsample(session_time, N);
            ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;

            % establish baseline windows
            [~,baseSt] = min(abs(ts1 - (baseline(1))));
            [~,baseEn] = min(abs(ts1 - (baseline(2))));
            [~,ampSt] = min(abs(ts1 - (amp_window(1))));
            [~,ampEn] = min(abs(ts1 - (amp_window(2))));
            [~,aucSt] = min(abs(ts1 - (auc_window(1))));
            [~,aucEn] = min(abs(ts1 - (auc_window(2))));

            %detrend & dFF%
            % bls = polyfit(ISOS_raw,SIGNAL_raw,1);
            % Y_fit_all = bls(1) .* ISOS_raw + bls(2);
            % Y_dF_all = SIGNAL_raw - Y_fit_all; %dF (units mV) is not dFF
            % dFF = 100*(Y_dF_all)./Y_fit_all;
            % std_dFF = std(double(dFF));
            % detrend_465 = detrend(dFF);
            % z465 = zscore(SIGNAL_raw);
            detrend_465 = detrend(SIGNAL_raw);

            for e = 1:height(leverArray)
                cueBase1 = cueArray(e,1) - (baseline(2));
                cueBase2 = cueArray(e,1) - (baseline(1));
                [~,cueSt] = min(abs(session_time - cueBase1));
                [~,cueEn] = min(abs(session_time - cueBase2));
                cueBaseMean(e,1) = mean(SIGNAL_raw(1,cueSt:cueEn));
                cueBaseStd(e,1) = std(SIGNAL_raw(1,cueSt:cueEn));
                leverStart = leverArray(e,1) - baseWindow;
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
                DLS_levers_raw(e,:) = leverSigRaw;
            end
            DLS_sessionSTREAMSraw{i,1} = DLS_levers_raw;
            DLS_amp_cueBase = [];
            DLS_auc_cueBase = [];
            for f = 1:height(DLS_levers_raw)
                % dF/F
                meanCue = cueBaseMean(f,1);
                stdCue = cueBaseStd(f,1);
                
                DLS_levers_dFF(f,1:epocArrayLen) = DLS_levers_raw(f, 1:epocArrayLen) - meanCue;
                DLS_levers_dFF(f,1:epocArrayLen) = 100*(DLS_levers_dFF(f, 1:epocArrayLen) / meanCue);
                % z-Score
                meanCueBase_dFF = mean(DLS_levers_dFF(f,baseSt:baseEn));
                stdCueBase_dFF = std(DLS_levers_dFF(f,baseSt:baseEn));
                DLS_levers_z(f,1:epocArrayLen) = (DLS_levers_dFF(f,1:epocArrayLen) - meanCueBase_dFF) / stdCueBase_dFF;
                DLS_amp_cueBase(f,1) = max(DLS_levers_z(f,ampSt:ampEn));
                DLS_auc_cueBase(f,1) = abs(trapz(ts1(1,aucSt:aucEn),DLS_levers_z(f,aucSt:aucEn)));
            end
            DLS_sessionSTREAMSz{i,1} = mean(DLS_levers_z);
            trialNames = array2table(session_identifiers(2:2:end,2),'VariableNames',{'Trial_Type'});
            DLS_levers_z = array2table(DLS_levers_z);
            DLS_levers_z = horzcat(trialNames,DLS_levers_z);
            DLS_cRew_cueBase = table2array(DLS_levers_z(DLS_levers_z.Trial_Type == 1, 2:end));
            DLS_cNoRew_cueBase = table2array(DLS_levers_z(DLS_levers_z.Trial_Type == 2, 2:end));
            DLS_iRew_cueBase = table2array(DLS_levers_z(DLS_levers_z.Trial_Type == 3, 2:end));
            DLS_iNoRew_cueBase = table2array(DLS_levers_z(DLS_levers_z.Trial_Type == 4, 2:end));
        
            DLS_amp_cueBase = array2table(DLS_amp_cueBase);
            DLS_amp_cueBase = horzcat(trialNames,DLS_amp_cueBase);
            DLS_auc_cueBase = array2table(DLS_auc_cueBase);
            DLS_auc_cueBase = horzcat(trialNames,DLS_auc_cueBase);

            for g = 1:height(cueArray)
                windowStart = cueArray(g,:) - baseWindow;
                windowEnd = windowStart + timeWindow + baseWindow;
                [~,windSt] = min(abs(session_time - windowStart(1)));
                [~,windEn] = min(abs(session_time - windowEnd(1)));
                cueSigRaw = SIGNAL_raw(1, windSt:windEn);
                if length(cueSigRaw) < epocArrayLen
                    mn = mean(cueSigRaw(1,end-10:end));
                    cueSigRaw(1,end:epocArrayLen) = mn;
                elseif length(cueSigRaw) > epocArrayLen
                    op = length(cueSigRaw);
                    arrayDif = op - epocArrayLen;
                    cueSigRaw = cueSigRaw(1,1:end-arrayDif);
                end
                DLS_cue_raw(g,1:epocArrayLen) = cueSigRaw;
            end
            [~,baseSt] = min(abs(ts1 - (baseline(1))));
            [~,baseEn] = min(abs(ts1 - (baseline(2))));

            for h = 1:height(DLS_cue_raw)
                % dF/F
                meanBase = mean(DLS_cue_raw(h,baseSt:baseEn));
                stdBase = std(DLS_cue_raw(h,baseSt:baseEn));
                DLS_cue_dFF(h,1:epocArrayLen) = DLS_cue_raw(h,1:epocArrayLen) - meanBase;
                DLS_cue_dFF(h,1:epocArrayLen) = 100*(DLS_cue_raw(h,1:epocArrayLen) / meanBase);
                % z-score
                meanBase_dFF = mean(DLS_cue_dFF(h,baseSt:baseEn));
                stdBase_dFF = std(DLS_cue_dFF(h,baseSt:baseEn));
                DLS_cue_z(h,1:epocArrayLen) = (DLS_cue_dFF(h,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
                % amplitude
                DLS_cue_amp_z(h) = max(DLS_cue_z(h,:));
                % AUC
                DLS_cue_auc_z(h) = trapz(ts1,DLS_cue_z(h,:));
            end
            DLS_cue_zALL{i,1} = DLS_cue_z(:,1:epocArrayLen);
            DLS_cue_ampALL{i,1} = DLS_cue_amp_z;
            DLS_cue_aucALL{i,1} = DLS_cue_auc_z;
        elseif k == 2
            ISOS = 'x405C';
            SIGNAL = 'x465C';
            ISOS_raw = data.streams.(ISOS).data;
            SIGNAL_raw = data.streams.(SIGNAL).data;
            %time array used for all streams%
            session_time = (1:length(SIGNAL_raw))/data.streams.(SIGNAL).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            ind = find(session_time>t,1);% find first index of when time crosses threshold
            session_time = session_time(ind:end); % reformat vector to only include allowed time
            SIGNAL_raw = data.streams.(SIGNAL).data(ind:end);
            ISOS_raw = data.streams.(ISOS).data(ind:end);
            minLength = min(length(ISOS_raw),length(SIGNAL_raw));
            ISOS_raw = ISOS_raw(1:minLength);
            SIGNAL_raw = SIGNAL_raw(1:minLength);
            
            %downsample streams and time array by N times%
            ISOS_raw = downsample(ISOS_raw, N);
            SIGNAL_raw = downsample(SIGNAL_raw, N);
            session_time = downsample(session_time, N);
            ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;

            % establish baseline windows
            [~,baseSt] = min(abs(ts1 - (baseline(1))));
            [~,baseEn] = min(abs(ts1 - (baseline(2))));
            [~,ampSt] = min(abs(ts1 - (amp_window(1))));
            [~,ampEn] = min(abs(ts1 - (amp_window(2))));
            [~,aucSt] = min(abs(ts1 - (auc_window(1))));
            [~,aucEn] = min(abs(ts1 - (auc_window(2))));
            
            %detrend & dFF%
            % bls = polyfit(ISOS_raw,SIGNAL_raw,1);
            % Y_fit_all = bls(1) .* ISOS_raw + bls(2);
            % Y_dF_all = SIGNAL_raw - Y_fit_all; %dF (units mV) is not dFF
            % dFF = 100*(Y_dF_all)./Y_fit_all;
            % std_dFF = std(double(dFF));
            % detrend_465 = detrend(dFF);
            % z465 = zscore(SIGNAL_raw);
            detrend_465 = detrend(SIGNAL_raw);

            for e = 1:height(leverArray)
                cueBase1 = cueArray(e,1) - (baseline(2));
                cueBase2 = cueArray(e,1) - (baseline(1));
                [~,cueSt] = min(abs(session_time - cueBase1));
                [~,cueEn] = min(abs(session_time - cueBase2));
                cueBaseMean(e,1) = mean(SIGNAL_raw(1,cueSt:cueEn));
                cueBaseStd(e,1) = std(SIGNAL_raw(1,cueSt:cueEn));
                leverStart = leverArray(e,1) - baseWindow;
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
                NAc_levers_raw(e,:) = leverSigRaw;
            end
            NAc_sessionSTREAMSraw{i,1} = NAc_levers_raw;
            NAc_amp_cueBase = [];
            NAc_auc_cueBase = [];
            for f = 1:height(NAc_levers_raw)
                % dF/F
                meanCue = cueBaseMean(f,1);
                stdCue = cueBaseStd(f,1);
                
                NAc_levers_dFF(f,1:epocArrayLen) = NAc_levers_raw(f, 1:epocArrayLen) - meanCue;
                NAc_levers_dFF(f,1:epocArrayLen) = 100*(NAc_levers_dFF(f,1:epocArrayLen) / meanCue);
                % z-Score
                meanCueBase_dFF = mean(NAc_levers_dFF(f,baseSt:baseEn));
                stdCueBase_dFF = std(NAc_levers_dFF(f,baseSt:baseEn));
                NAc_levers_z(f,1:epocArrayLen) = (NAc_levers_dFF(f,1:epocArrayLen) - meanCueBase_dFF) / stdCueBase_dFF;
                NAc_amp_cueBase(f,1) = max(NAc_levers_z(f,ampSt:ampEn));
                NAc_auc_cueBase(f,1) = abs(trapz(ts1(1,aucSt:aucEn),NAc_levers_z(f,aucSt:aucEn)));
            end
            NAc_sessionSTREAMSz{i,1} = mean(NAc_levers_z);
            trialNames = array2table(session_identifiers(2:2:end,2),'VariableNames',{'Trial_Type'});
            NAc_levers_z = array2table(NAc_levers_z);
            NAc_levers_z = horzcat(trialNames,NAc_levers_z);
            NAc_cRew_cueBase = table2array(NAc_levers_z(NAc_levers_z.Trial_Type == 1, 2:end));
            NAc_cNoRew_cueBase = table2array(NAc_levers_z(NAc_levers_z.Trial_Type == 2, 2:end));
            NAc_iRew_cueBase = table2array(NAc_levers_z(NAc_levers_z.Trial_Type == 3, 2:end));
            NAc_iNoRew_cueBase = table2array(NAc_levers_z(NAc_levers_z.Trial_Type == 4, 2:end));
            
            NAc_amp_cueBase = array2table(NAc_amp_cueBase);
            NAc_amp_cueBase = horzcat(trialNames,NAc_amp_cueBase);
            NAc_auc_cueBase = array2table(NAc_auc_cueBase);
            NAc_auc_cueBase = horzcat(trialNames,NAc_auc_cueBase);
            for g = 1:height(cueArray)
                windowStart = cueArray(g,:)-baseWindow;
                windowEnd = windowStart + timeWindow + baseWindow;
                [~,windSt] = min(abs(session_time - windowStart(1)));
                [~,windEn] = min(abs(session_time - windowEnd(1)));
                cueSigRaw = SIGNAL_raw(1, windSt:windEn);
                if length(cueSigRaw) < epocArrayLen
                    mn = mean(cueSigRaw(1,end-10:end));
                    cueSigRaw(1,end:epocArrayLen) = mn;
                elseif length(cueSigRaw) > epocArrayLen
                    op = length(cueSigRaw);
                    arrayDif = op - epocArrayLen;
                    cueSigRaw = cueSigRaw(1,1:end-arrayDif);
                end
                NAc_cue_raw(g,1:epocArrayLen) = cueSigRaw;
            end
            [~,baseSt] = min(abs(ts1 - (baseline(1))));
            [~,baseEn] = min(abs(ts1 - (baseline(2))));
            
            for h = 1:height(NAc_cue_raw)
                % dF/F
                meanBase = mean(NAc_cue_raw(h,baseSt:baseEn));
                stdBase = std(NAc_cue_raw(h,baseSt:baseEn));
                NAc_cue_dFF(h,1:epocArrayLen) = NAc_cue_raw(h,1:epocArrayLen) - meanBase;
                NAc_cue_dFF(h,1:epocArrayLen) = 100*(NAc_cue_raw(h,1:epocArrayLen) / meanBase);
                % z-score
                meanBase_dFF = mean(NAc_cue_dFF(h,baseSt:baseEn));
                stdBase_dFF = std(NAc_cue_dFF(h,baseSt:baseEn));
                NAc_cue_z(h,1:epocArrayLen) = (NAc_cue_dFF(h,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
                % amplitude
                NAc_cue_amp_z(h) = max(NAc_cue_z(h,:));
                % AUC
                NAc_cue_auc_z(h) = trapz(ts1,NAc_cue_z(h,:));
            end
            NAc_cue_zALL{i,1} = NAc_cue_z(:,1:epocArrayLen);
            NAc_cue_ampALL{i,1} = NAc_cue_amp_z;
            NAc_cue_aucALL{i,1} = NAc_cue_auc_z;
            
        end
    end
   
end
idList = cell2table(IDs', 'VariableNames',{'ID'});
phaseList = cell2table(phaseList', 'VariableNames',{'Phase'});
treatList = cell2table(treatList', 'VariableNames',{'Treatment'});

