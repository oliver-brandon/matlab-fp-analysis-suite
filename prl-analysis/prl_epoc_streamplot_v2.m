% prl_epoc_streamplot created by Brandon L. Oliver, M.A.
% (boliv018@ucr.edu)
%
% ABOUT:
% This script will plot single fiber PRL epoc averaged signal data and save
% to a user designated directory 'figure_savepath'. To use this, files must
% converted and separated into .mat files using the scripts in 'tank2mat'
% (see README in tank2mat).
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
% To begin, click 'Run' and select a folder containing one or more .mat
% files. Make sure you set your figure save path. It is recommended to
% create a specific folder on your computer for the figures to output to.
clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 5; % baseline signal to include before TTL 
baselineZ_cue = [3 1];
baselineZ_lever = [3 1];
N = 1; %Downsample N times
sigHz = 1017/N;
minArrayLen = round(sigHz * (timeWindow + baseline));
figure_savepath = '/Users/brandon/Library/CloudStorage/GoogleDrive-boliv018@ucr.edu/My Drive/prl/PRL_GRABDA/testFigs/';
savetype = '.pdf';
VERSION = '2.0';
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

    [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cueTS,correct_rewarded,...
        correct_noreward,incorrect_rewarded,incorrect_noreward);
    % winStayTS = errorProbExtract(trial_type,session_ts,1,1);
    % winShiftTS = errorProbExtract(trial_type,session_ts,2,1);
    % loseStayTS = errorProbExtract(trial_type,session_ts,3,1);
    % loseShiftTS = errorProbExtract(trial_type,session_ts,4,1);

    % [epocSTREAM,~,~] = epocExtract( ...
    %     detrend_465, ...
    %     time, ...
    %     session_ts, ...
    %     baseline, ...
    %     timeWindow, ...
    %     baselineZ_cue, ...
    %     [0,5], ...
    %     minArrayLen ...
    %     );
    % if ~isempty(epocSTREAM)
    %     sessionSTREAM = epocSTREAM;
    % else
    %     disp('')
    % end
    % 
    [epocSTREAM,~,~] = epocExtract( ...
        z465, ...
        time, ...
        cueTS, ...
        baseline, ...
        timeWindow, ...
        baselineZ_cue, ...
        [0,2], ...
        minArrayLen ...
        );
    if ~isempty(epocSTREAM)
        cueSTREAM = epocSTREAM;
    else
        disp('')
    end

    [epocSTREAM,~,~] = epocExtract( ...
        detrend_465, ...
        time, ...
        correct_rewarded, ...
        baseline, ...
        timeWindow, ...
        baselineZ_lever, ...
        [0,2], ...
        minArrayLen ...
        );
    if ~isempty(epocSTREAM)
        cRewSTREAM = epocSTREAM;
    else
        disp('')
    end

    [epocSTREAM,~,~] = epocExtract( ...
        detrend_465, ...
        time, ...
        correct_noreward, ...
        baseline, ...
        timeWindow, ...
        baselineZ_lever, ...
        [2,4], ...
        minArrayLen ...
        );
    if ~isempty(epocSTREAM)
        cNoRewSTREAM = epocSTREAM;
    else
        disp('')
    end

    [epocSTREAM,~,~] = epocExtract( ...
        detrend_465, ...
        time, ...
        incorrect_rewarded, ...
        baseline, ...
        timeWindow, ...
        baselineZ_lever, ...
        [0,2], ...
        minArrayLen ...
        );
    if ~isempty(epocSTREAM)
        iRewSTREAM = epocSTREAM;
    else
        disp('')
    end

    [epocSTREAM,~,~] = epocExtract( ...
        detrend_465, ...
        time, ...
        incorrect_noreward, ...
        baseline, ...
        timeWindow, ...
        baselineZ_lever, ...
        [2,4], ...
        minArrayLen ...
        );
    if ~isempty(epocSTREAM)
        iNoRewSTREAM = epocSTREAM;
    else
        disp('')
    end

    % [epocSTREAM,~,~] = epocExtract( ...
    %     detrend_465, ...
    %     time, ...
    %     winStayTS, ...
    %     baseline, ...
    %     timeWindow, ...
    %     baselineZ_lever, ...
    %     [0,2], ...
    %     minArrayLen ...
    %     );
    % if ~isempty(epocSTREAM)
    %     winStaySTREAM = epocSTREAM;
    % else
    %     disp('')
    % end
    % 
    % [epocSTREAM,~,~] = epocExtract( ...
    %     detrend_465, ...
    %     time, ...
    %     loseShiftTS, ...
    %     baseline, ...
    %     timeWindow, ...
    %     baselineZ_lever, ...
    %     [0,2], ...
    %     minArrayLen ...
    %     );
    % if ~isempty(epocSTREAM)
    %     loseShiftSTREAM = epocSTREAM;
    % else
    %     disp('')
    % end

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
    if exist('winStaySTREAM','var')
        epoc_stream{end+1} = double(winStaySTREAM);
        epoc_name{end+1} = char('winStay');
    end
    if exist('loseShiftSTREAM','var')
        epoc_stream{end+1} = double(loseShiftSTREAM);
        epoc_name{end+1} = char('loseShift');
    end
    pmchar = char(177);
    
    % Loop through the cell of arrays and plot mean and standard error
    for i = 1:length(epoc_stream)
        if height(epoc_stream{i}) > 1
            fig{i} = figure;
            % Mean and standard error
            mean_data = mean(epoc_stream{i}, 1);
            se_data = std(epoc_stream{i}, 0, 1) / sqrt(size(epoc_stream{i}, 1));
            dc_data = mean(mean_data);
            % Plot mean with blue line
            plot(ts1, mean_data, 'b-', 'LineWidth', 3);
            xlim([-baseline max(ts1)]);
            hold on
            
            % Plot standard error with red line (50% translucence)
            XX = [ts1, fliplr(ts1)];
            YY = [mean_data + se_data, fliplr(mean_data - se_data)];
            
            % Plot filled standard error bands.
            h = fill(XX, YY, 'r');
            set(h, 'facealpha',.25,'edgecolor','none')
            % Plot line at x=0
            xline(0,'LineWidth',2,'Color','black')
            yline(0,'LineWidth',2,'Color','black')
            
            % Axis labels
            xlabel('Time (s)', 'FontSize',12);
            ylabel(['zScore GrabDA ' pmchar 'SEM'],'FontSize',12);
            TITLE{i} = strcat(animalID," ",task," ",epoc_name{i}," ",treatment);
        
        
            % Legend
            title(TITLE{i},'FontSize',18);
            legend({'Mean', 'Standard Error'}, 'Location', 'best');
            savepath = strcat(figure_savepath,animalID,"/");
            file_name = strcat(savepath,TITLE{i},savetype);
            if not(isfolder(savepath))
                mkdir(savepath);
            end
            print(fig{i}, file_name, '-dpdf', '-bestfit');
            % saveas(fig{i},file_name)
        else
            continue
        end
    end
    
end
toc
NERD_STATS(toc,numFiles);