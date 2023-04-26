% prl_df_epoc_streamplot created by Brandon L. Oliver, M.A.
% (boliv018@ucr.edu)
%
% ABOUT:
% This script will plot dual fiber PRL epoc averaged signal data and save
% to a user designated directory 'figure_savepath'. To use this, files must
% be named either 'Empty_NA_NA_ID_Task_Treatment' (recorded in boxes 2, 4,
% or 6) OR 'ID_Task_Treatment_Empty_NA_NA' (recorded in boxes 1, 3, or 5)
% depending on which box the animal was recorded in. Currently, this script
% only accepts TDT tanks.
%
% IMPORTANT:
% This script uses functions outside of this script that are included in
% the 'functions' directory contained within the entire repositiory. Thus, 
% it is recommended to download the entire repository this script is apart 
% of by visiting https://github.com/oliver-brandon/matlab-fp-analysis-suite. 
% Once downloaded, add the entire matlab-fp-analysis-suite folder to your 
% MATLAB path. 
%
% INSTRUCTIONS:
% To begin, click 'Run' and select a folder containing one or more TDT
% tanks. Make sure you set your figure save path. It is recommended to
% create a specific folder on your computer for the figures to output to.

clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ_cue = [5 1];
baselineZ_lever = [5 1];
N = 100; %Downsample N times
minArrayLen = 72; %timeWindow = 5 - 72, timeWindow = 10 - 123
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
figure_savepath = '/Volumes/CUDADRIVE/DA_PRL/PRL_DF_EpocPlots/test/';% must include backslash at end of path
savetype = '.pdf'; % can set to desired file type
stream_A = 'DLS';
stream_C = 'NAc';
VERSION = 'v1.1b';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(pwd,"Select a folder containing one or more files"); 
fprintf("prl_df_epoc_plotsaver Version: %s\n",VERSION)
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
    disp("Starting single plot...")
end

for batch = 1:numFiles
    fprintf('Loading file %d of %d...\n',batch,length(myFiles))
    filename = fullfile(myDir, myFiles(batch).name);
    data = TDTbin2mat(filename,'TYPE',{'streams','epocs'});
    [~,name,~] = fileparts(filename);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    animalIDA = char(brokenID{1});
    animalIDC = char(brokenID{4});
    emptylogicA = strcmp(animalIDA,emptyID);
    emptylogicC = strcmp(animalIDC,emptyID);
    if emptylogicA == 0
        TTLs = 1;
        animalID = animalIDA;
        task = char(brokenID{2});
        treatment = char(brokenID{3});
    elseif emptylogicA == 1
        TTLs = 2;
        animalID = animalIDC;
        task = char(brokenID{5});
        treatment = char(brokenID{6});
    end
    data = prl_df_epocs(data,TTLs);

    for k = 1:2
        if k == 1
            % cueTS = [];
            % correct_rewarded = [];
            % correct_noreward = [];
            % incorrect_rewarded = [];
            % incorrect_noreward = [];
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
                e3 = cueTS(e,1)-baselineZ_cue(:,2);
                e4 = cueTS(e,1)-baselineZ_cue(:,1);
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
        
            for e = 1:height(correct_rewarded)
                e1 = correct_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_rewarded(e,1)-baselineZ_lever(:,2);
                e4 = correct_rewarded(e,1)-baselineZ_lever(:,1);
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
                    %cRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                cRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
        
            for e = 1:height(correct_noreward)
                e1 = correct_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_noreward(e,1)-baselineZ_lever(:,2);
                e4 = correct_noreward(e,1)-baselineZ_lever(:,1);
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
                    %cNoRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
        
            for e = 1:height(incorrect_rewarded)
                e1 = incorrect_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_rewarded(e,1)-baselineZ_lever(:,2);
                e4 = incorrect_rewarded(e,1)-baselineZ_lever(:,1);
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
                    %iRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                iRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
        
            for e = 1:height(incorrect_noreward)
                e1 = incorrect_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_noreward(e,1)-baselineZ_lever(:,2);
                e4 = incorrect_noreward(e,1)-baselineZ_lever(:,1);
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
                    %iNoRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
            epoc_streamA = {};
            epoc_nameA = {};
            if exist('cueSTREAM','var')
                epoc_streamA{end+1} = double(cueSTREAM);
                epoc_nameA{end+1} = char('Cue');
            end
            if correct_rewarded > 0
                epoc_streamA{end+1} = double(cRewSTREAM);
                epoc_nameA{end+1} = char('C+');
            else
                continue
            end
            if correct_noreward > 0
                epoc_streamA{end+1} = double(cNoRewSTREAM);
                epoc_nameA{end+1} = char('C-');
            else
                continue
            end
            if incorrect_rewarded > 0
                epoc_streamA{end+1} = double(iRewSTREAM);
                epoc_nameA{end+1} = char('I+');
            else
                continue
            end
            if incorrect_noreward > 0
                epoc_streamA{end+1} = double(iNoRewSTREAM);
                epoc_nameA{end+1} = char('I-');
            else
                continue
            end
            
    
        elseif k == 2
            % cueTS = [];
            % correct_rewarded = [];
            % correct_noreward = [];
            % incorrect_rewarded = [];
            % incorrect_noreward = [];
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
                e3 = cueTS(e,1)-baselineZ_cue(:,2);
                e4 = cueTS(e,1)-baselineZ_cue(:,1);
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
            
            for e = 1:height(correct_rewarded)
                e1 = correct_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_rewarded(e,1)-baselineZ_lever(:,2);
                e4 = correct_rewarded(e,1)-baselineZ_lever(:,1);
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
                    %cRewSTREAM(1:minArrayLen,1) = zeros;
        
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                cRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
            
            for e = 1:height(correct_noreward)
                e1 = correct_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = correct_noreward(e,1)-baselineZ_lever(:,2);
                e4 = correct_noreward(e,1)-baselineZ_lever(:,1);
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
                    %cNoRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                cNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
            for e = 1:height(incorrect_rewarded)
                e1 = incorrect_rewarded(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_rewarded(e,1)-baselineZ_lever(:,2);
                e4 = incorrect_rewarded(e,1)-baselineZ_lever(:,1);
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
                    %iRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                if length(zfinal) < minArrayLen
                    continue
                end
                iRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
            
            for e = 1:height(incorrect_noreward)
                e1 = incorrect_noreward(e,1)-baseline;
                e2 = e1+timeWindow+baseline;
                if e1 == 0
                    continue
                end
                e3 = incorrect_noreward(e,1)-baselineZ_lever(:,2);
                e4 = incorrect_noreward(e,1)-baselineZ_lever(:,1);
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
                    %iNoRewSTREAM(1:minArrayLen,1) = zeros;
                    continue
                end
                iNoRewSTREAM(e,:) = zfinal(1:minArrayLen);
            end
        

            epoc_streamC = {};
            epoc_nameC = {};
            if exist('cueSTREAM','var')
                epoc_streamC{end+1} = double(cueSTREAM);
                epoc_nameC{end+1} = char('Cue');
            end
            if correct_rewarded > 0
                epoc_streamC{end+1} = double(cRewSTREAM);
                epoc_nameC{end+1} = char('C+');
            else
                continue
            end
            if correct_noreward > 0
                epoc_streamC{end+1} = double(cNoRewSTREAM);
                epoc_nameC{end+1} = char('C-');
            else
                continue
            end
            if incorrect_rewarded > 0
                epoc_streamC{end+1} = double(iRewSTREAM);
                epoc_nameC{end+1} = char('I+');
            else
                continue
            end
            if incorrect_noreward > 0
                epoc_streamC{end+1} = double(iNoRewSTREAM);
                epoc_nameC{end+1} = char('I-');
            else
                continue
            end
        end
    end
    pmchar = char(177);
        
    % Loop through the cell of arrays and plot mean and standard error
    for i = 1:length(epoc_streamA)
        fig{i} = figure;
        % Mean and standard error
        mean_dataA = mean(epoc_streamA{i}, 1);
        mean_dataC = mean(epoc_streamC{i}, 1);
        se_dataA = std(epoc_streamA{i}, 0, 1) / sqrt(size(epoc_streamA{i}, 1));
        se_dataC = std(epoc_streamC{i}, 0, 1) / sqrt(size(epoc_streamC{i}, 1));
        dc_dataA = mean(mean_dataA);
        dc_dataC = mean(mean_dataC);
        % Plot mean with blue line
        plot(ts1, mean_dataA, 'b-', 'LineWidth', 3, 'DisplayName','Mean DLS');
        hold on
        plot(ts1, mean_dataC, 'r-', 'LineWidth', 3, 'DisplayName','Mean NAc');
       
        xlim([-2 max(ts1)]);
        
        
        % Plot standard error with red line (50% translucence)
        XX = [ts1, fliplr(ts1)];
        YY_A = [mean_dataA + se_dataA, fliplr(mean_dataA - se_dataA)];
        YY_C = [mean_dataC + se_dataC, fliplr(mean_dataC - se_dataC)];
        
        % Plot filled standard error bands.
        h = fill(XX, YY_A, 'r', 'DisplayName', 'SEM DLS');
        g = fill(XX, YY_C, 'b', 'DisplayName','SEM NAc');
        set(h, 'facealpha',.20,'edgecolor','none')
        set(g, 'facealpha',.20,'edgecolor','none')
        
        % Axis labels
        xlabel('Time (s)', 'FontSize',12);
        ylabel(['GrabDA \DeltaF/F ' pmchar 'SEM'],'FontSize',12);
        if TTLs == 1
            epoc_name = epoc_nameA;
        elseif TTLs == 2
            epoc_name = epoc_nameC;
        end
        % Plot line at x=0
        xline(0,'LineWidth',2,'Color','black','DisplayName','Epoc Onset')
        % Create graph title
        TITLE{i} = strcat(animalID," ",task," ",epoc_name{i});
    
    
        % Legend
        title(TITLE{i},'FontSize',18);
        legend('Location', 'best');
        savepath = strcat(figure_savepath,animalID,"/");
        file_name = strcat(savepath,TITLE{i},savetype);
        if not(isfolder(savepath))
            mkdir(savepath);
        end
        saveas(fig{i},file_name)
        
    end
     
end
toc
NERD_STATS(toc,numFiles);