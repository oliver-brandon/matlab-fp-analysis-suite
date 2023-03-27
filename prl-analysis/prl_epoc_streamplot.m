clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ_cue = [5 1];
baselineZ_lever = [3 1];
N = 100; %Downsample N times
minArrayLen = 72; %timeWindow = 5 - 72, timeWindow = 10 - 123
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
figure_savepath = 'Z:\DA_PRL\Figs\';
savetype = ".pdf";
VERSION = '1.0';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(pwd,"Select a folder containing one or more files"); 
fprintf("prl_epoc_plotsaver Version: %s\n",VERSION)
tic
if myDir == 0
    disp("Select a folder containing one or more files")
    return
end
myFiles = dir(myDir);
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv'}));
numFiles = length(myFiles);
fig = cell(length(myFiles),1);
TITLE = cell(length(myFiles),1);
if length(myFiles) > 1
    disp("Starting batch plot...")
elseif length(myFiles) == 1
    disp("Starting single plot")
end

for batch = 1:numFiles
    fprintf('Loading file %d of %d...\n',batch,length(myFiles))
    filename = fullfile(myDir, myFiles(batch).name);
    load(filename)
    [~,name,~] = fileparts(filename);
    brokenID = strsplit(name,'_');
    animalID = char(brokenID{1});
    task = char(brokenID{2});
    treatment = char(brokenID{3});

    if isfield(data.streams, 'x405A')
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
        for e = 1:height(cueTSA)
            e1 = data.epocs.St1_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St1_.onset(e,1)-baselineZ_cue(:,2);
            e4 = data.epocs.St1_.onset(e,1)-baselineZ_cue(:,1);
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
                continue
            end
            cueSTREAM(e,:) = zfinal(1:minArrayLen);
            
        end
    
        for e = 1:height(correct_rewardedA)
            e1 = correct_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedA(e,1)-baselineZ_lever(:,2);
            e4 = correct_rewardedA(e,1)-baselineZ_lever(:,1);
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
                %cRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
    
        for e = 1:height(correct_norewardA)
            e1 = correct_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardA(e,1)-baselineZ_lever(:,2);
            e4 = correct_norewardA(e,1)-baselineZ_lever(:,1);
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
                %cNoRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
    
        for e = 1:height(incorrect_rewardedA)
            e1 = incorrect_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedA(e,1)-baselineZ_lever(:,2);
            e4 = incorrect_rewardedA(e,1)-baselineZ_lever(:,1);
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
                %iRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            iRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
    
        for e = 1:height(incorrect_norewardA)
            e1 = incorrect_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardA(e,1)-baselineZ_lever(:,2);
            e4 = incorrect_norewardA(e,1)-baselineZ_lever(:,1);
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
            if incorrect_norewardA == 0
                %iNoRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
    
    else
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
        for e = 1:height(cueTSC)
            e1 = data.epocs.St2_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St2_.onset(e,1)-baselineZ_cue(:,2);
            e4 = data.epocs.St2_.onset(e,1)-baselineZ_cue(:,1);
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
                continue
            end
            cueSTREAM(e,:) = zfinal(1:minArrayLen);
        end
        
        for e = 1:height(correct_rewardedC)
            e1 = correct_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedC(e,1)-baselineZ_lever(:,2);
            e4 = correct_rewardedC(e,1)-baselineZ_lever(:,1);
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
            if correct_rewardedC == 0
                %cRewSTREAM(1:minArrayLen,1) = zeros;
    
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
        
        for e = 1:height(correct_norewardC)
            e1 = correct_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardC(e,1)-baselineZ_lever(:,2);
            e4 = correct_norewardC(e,1)-baselineZ_lever(:,1);
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
            if correct_norewardC == 0
                %cNoRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
        for e = 1:height(incorrect_rewardedC)
            e1 = incorrect_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedC(e,1)-baselineZ_lever(:,2);
            e4 = incorrect_rewardedC(e,1)-baselineZ_lever(:,1);
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
                %iRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            if length(zfinal) < minArrayLen
                break
            end
            iRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
        
        for e = 1:height(incorrect_norewardC)
            e1 = incorrect_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardC(e,1)-baselineZ_lever(:,2);
            e4 = incorrect_norewardC(e,1)-baselineZ_lever(:,1);
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
                %iNoRewSTREAM(1:minArrayLen,1) = zeros;
                break
            end
            iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
        end
    end

    epoc_stream = {};
    epoc_name = {};
    if exist('cueSTREAM','var')
        epoc_stream{end+1} = double(cueSTREAM);
        epoc_name{end+1} = char('Cue');
    end
    if exist('cRewSTREAM','var')
        epoc_stream{end+1} = double(cRewSTREAM);
        epoc_name{end+1} = char('C+');
    end
    if exist('cNoRewSTREAM','var')
        epoc_stream{end+1} = double(cNoRewSTREAM);
        epoc_name{end+1} = char('C-');
    end
    if exist('iRewSTREAM','var')
        epoc_stream{end+1} = double(iRewSTREAM);
        epoc_name{end+1} = char('I+');
    end
    if exist('iNoRewSTREAM','var')
        epoc_stream{end+1} = double(iNoRewSTREAM);
        epoc_name{end+1} = char('I-');
    end
    pmchar = char(177);

    % Loop through the cell of arrays and plot mean and standard error
    for i = 1:length(epoc_stream)
        fig{i} = figure;
        % Mean and standard error
        mean_data = mean(epoc_stream{i}, 1);
        se_data = std(epoc_stream{i}, 0, 1) / sqrt(size(epoc_stream{i}, 1));
        dc_data = mean(mean_data);
        % Plot mean with blue line
        plot(ts1, mean_data, 'b-', 'LineWidth', 3);
        xlim([-2 max(ts1)]);
        hold on
        
        % Plot standard error with red line (50% translucence)
        XX = [ts1, fliplr(ts1)];
        YY = [mean_data + se_data, fliplr(mean_data - se_data)];
        
        % Plot filled standard error bands.
        h = fill(XX, YY, 'r');
        set(h, 'facealpha',.25,'edgecolor','none')
        % Plot line at x=0
        xline(0,'LineWidth',2,'Color','black')
        
        % Axis labels
        xlabel('Time (s)', 'FontSize',12);
        ylabel(['zScore GrabDA ' pmchar 'SEM'],'FontSize',12);
        TITLE{i} = strcat(animalID," ",task," ",epoc_name{i});


        % Legend
        title(TITLE{i},'FontSize',18);
        legend({'Mean', 'Standard Error'}, 'Location', 'best');
        savepath = strcat(figure_savepath,animalID,"\");
        file_name = strcat(savepath,TITLE{i},savetype);
        if not(isfolder(savepath))
            mkdir(savepath);
        end
        saveas(fig{i},file_name)
        
        
    end
end
toc
NERD_STATS(toc,numFiles);