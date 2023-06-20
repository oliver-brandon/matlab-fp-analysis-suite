clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [3 1]; % baseline signal for dFF/zscore
amp_window = [0 2]; % time window to grab amplitude from
auc_window = [-1 timeWindow];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/Library/CloudStorage/GoogleDrive-boliv018@ucr.edu/My Drive/prl/dual_fiber/tanks','Choose a folder containing the tank(s) you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a folder with tanks to start")
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
    for k = 1:2
        if k == 1
            ISOS = 'x405A';
            GRABDA = 'x465A';
            %time array used for all streams%
            time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            t = 5; % time threshold below which we will discard
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
            % ts1 = -baseline + (1:epocArrayLen) / data.streams.(GRABDA).fs*N;
            
            %detrend & dFF%
            bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
            Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
            Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
            dFF = 100*(Y_dF_all)./Y_fit_all;
            std_dFF = std(double(dFF));
            detrend_465 = detrend(dFF);
            z465 = zscore(data.streams.(GRABDA).data);
            if TTLs == 1
                cueTS = data.epocs.St1_.onset;
                correct_rewarded = data.epocs.cRewA.onset;
                correct_noreward = data.epocs.cNoRewA.onset;
                incorrect_rewarded = data.epocs.iRewA.onset;
                incorrect_noreward = data.epocs.iNoRewA.onset;
            elseif TTLs == 2
                cueTS = data.epocs.St2_.onset;
                correct_rewarded = data.epocs.cRewC.onset;
                correct_noreward = data.epocs.cNoRewC.onset;
                incorrect_rewarded = data.epocs.iRewC.onset;
                incorrect_noreward = data.epocs.iNoRewC.onset;
            end
            for e = 1:height(cueTS)
                e1 = cueTS(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                e3 = cueTS(e,1)-baseline(:,2);
                e4 = cueTS(e,1)-baseline(:,1);
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
                if length(zfinal) < epocArrayLen
                    continue
                end
                DLS_cueSTREAM(e,:) = zfinal(1:epocArrayLen);
                
            end
        
            for e = 1:height(correct_rewarded)
                e1 = correct_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_rewarded(e,1)-baseline(:,2);
                e4 = correct_rewarded(e,1)-baseline(:,1);
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
                if correct_rewarded == 0
                    %cRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                DLS_cRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
        
            for e = 1:height(correct_noreward)
                e1 = correct_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_noreward(e,1)-baseline(:,2);
                e4 = correct_noreward(e,1)-baseline(:,1);
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
                if correct_noreward == 0
                    %cNoRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                DLS_cNoRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
        
            for e = 1:height(incorrect_rewarded)
                e1 = incorrect_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_rewarded(e,1)-baseline(:,2);
                e4 = incorrect_rewarded(e,1)-baseline(:,1);
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
                if incorrect_rewarded == 0
                    %iRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                DLS_iRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
        
            for e = 1:height(incorrect_noreward)
                e1 = incorrect_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_noreward(e,1)-baseline(:,2);
                e4 = incorrect_noreward(e,1)-baseline(:,1);
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
                if incorrect_noreward == 0
                    %iNoRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                DLS_iNoRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
            epoc_streamA = {};
            epoc_nameA = {};
            if exist('cueSTREAM','var')
                epoc_streamA{end+1} = double(DLS_cueSTREAM);
                epoc_nameA{end+1} = char('Cue');
            end
            if correct_rewarded > 0
                epoc_streamA{end+1} = double(DLS_cRewSTREAM);
                epoc_nameA{end+1} = char('C+');
            else
                continue
            end
            if correct_noreward > 0
                epoc_streamA{end+1} = double(DLS_cNoRewSTREAM);
                epoc_nameA{end+1} = char('C-');
            else
                continue
            end
            if incorrect_rewarded > 0
                epoc_streamA{end+1} = double(DLS_iRewSTREAM);
                epoc_nameA{end+1} = char('I+');
            else
                continue
            end
            if incorrect_noreward > 0
                epoc_streamA{end+1} = double(DLS_iNoRewSTREAM);
                epoc_nameA{end+1} = char('I-');
            else
                continue
            end
            
    
        elseif k == 2
            ISOS = 'x405C';
            GRABDA = 'x465C';
            %time array used for all streams%
            time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            t = 5;% time threshold below which we will discard
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
            % ts1 = -baseline + (1:epocArrayLen) / data.streams.(GRABDA).fs*N;
            
            %detrend & dFF%
            bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
            Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
            Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
            dFF = 100*(Y_dF_all)./Y_fit_all;
            std_dFF = std(double(dFF));
            detrend_465 = detrend(dFF);
            z465 = zscore(data.streams.(GRABDA).data);
            if TTLs == 1
                cueTS = cueTS;
                correct_rewarded = data.epocs.cRewA.onset;
                correct_noreward = data.epocs.cNoRewA.onset;
                incorrect_rewarded = data.epocs.iRewA.onset;
                incorrect_noreward = data.epocs.iNoRewA.onset;
            elseif TTLs == 2
                cueTS = cueTS;
                correct_rewarded = data.epocs.cRewC.onset;
                correct_noreward = data.epocs.cNoRewC.onset;
                incorrect_rewarded = data.epocs.iRewC.onset;
                incorrect_noreward = data.epocs.iNoRewC.onset;
            end
            for e = 1:height(cueTS)
                e1 = cueTS(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                e3 = cueTS(e,1)-baseline(:,2);
                e4 = cueTS(e,1)-baseline(:,1);
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
                if length(zfinal) < epocArrayLen
                    continue
                end
                NAc_cueSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
            
            for e = 1:height(correct_rewarded)
                e1 = correct_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_rewarded(e,1)-baseline(:,2);
                e4 = correct_rewarded(e,1)-baseline(:,1);
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
                if correct_rewarded == 0
                    %cRewSTREAM(1:epocArrayLen,1) = zeros;
        
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                NAc_cRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
            
            for e = 1:height(correct_noreward)
                e1 = correct_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_noreward(e,1)-baseline(:,2);
                e4 = correct_noreward(e,1)-baseline(:,1);
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
                if correct_noreward == 0
                    %cNoRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                NAc_cNoRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
            for e = 1:height(incorrect_rewarded)
                e1 = incorrect_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_rewarded(e,1)-baseline(:,2);
                e4 = incorrect_rewarded(e,1)-baseline(:,1);
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
                if incorrect_rewarded == 0
                    %iRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < epocArrayLen
                    continue
                end
                NAc_iRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
            
            for e = 1:height(incorrect_noreward)
                e1 = incorrect_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_noreward(e,1)-baseline(:,2);
                e4 = incorrect_noreward(e,1)-baseline(:,1);
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
                if incorrect_noreward == 0
                    %iNoRewSTREAM(1:epocArrayLen,1) = zeros;
                    continue
                end
                NAc_iNoRewSTREAM(e,:) = zfinal(1:epocArrayLen);
            end
        

            epoc_streamC = {};
            epoc_nameC = {};
            if exist('cueSTREAM','var')
                epoc_streamC{end+1} = double(DLS_cueSTREAM);
                epoc_nameC{end+1} = char('Cue');
            end
            if correct_rewarded > 0
                epoc_streamC{end+1} = double(DLS_cRewSTREAM);
                epoc_nameC{end+1} = char('C+');
            else
                continue
            end
            if correct_noreward > 0
                epoc_streamC{end+1} = double(DLS_cNoRewSTREAM);
                epoc_nameC{end+1} = char('C-');
            else
                continue
            end
            if incorrect_rewarded > 0
                epoc_streamC{end+1} = double(DLS_iRewSTREAM);
                epoc_nameC{end+1} = char('I+');
            else
                continue
            end
            if incorrect_noreward > 0
                epoc_streamC{end+1} = double(DLS_iNoRewSTREAM);
                epoc_nameC{end+1} = char('I-');
            else
                continue
            end
        end
    end
     
end