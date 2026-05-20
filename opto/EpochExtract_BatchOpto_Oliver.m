clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Paramaters to Edit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STREAM_STORE1 = 'x405A';
STREAM_STORE2 = 'x465A';
smoothFactor = 20;
setBaseline = 1; % 1 = yes, 0 = no (adjusts signals to zero indicated by baseAdjust)
baseAdjust = -2; % seconds on x axis to adjust baseline to

TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-3 -1]; % baseline period before epoc

REF_EPOC = 'St1/'; % Stimulation event to center on
REF_EPOC2 = 'levers';

dataType = 2; % 1 = tank, 2 = mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Leave Code Below As Is %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir = uigetdir('/Users/brandon/personal-drive/optomouse-prime');
if dataDir == 0
    disp('Select a folder containing one or more .mat files to start');
    return
end
tic
dataFiles = dir(dataDir);
dataFiles = dataFiles(~startsWith({dataFiles.name},{'.','..','._'}));
dataFiles = dataFiles(endsWith({dataFiles.name},{'.mat'}));
numFiles = numel(dataFiles);
signals = struct("cues",[],"levers",[],"meta",[]);
allCues = cell(size(numFiles));
allLevers = cell(size(numFiles));
for k = 1:numFiles
    filename = fullfile(dataDir,dataFiles(k).name);
    [~,name,~] = fileparts(filename); % gets name of tank
    fprintf('Analyzing %s (%d of %d)\n', name, k, numFiles)
    
    if dataType == 1
        data = TDTbin2mat(filename, 'TYPE', {'epocs', 'streams'});% TDT function for extracting data to struct 'data'
    elseif dataType == 2
        load(filename)
    end
        
    if isfield(data.epocs, "events")
        data.epocs = rmfield(data.epocs, 'events');
    end
    % Use TDTfilter to extract data around our epoc event
    % Using the 'TIME' parameter extracts data only from the time range around
    % our epoc event. Use the 'VALUES' parameter to specify allowed values of
    % the REF_EPOC to extract.  For stream events, the chunks of data are 
    % stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
    data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE); % extracts data around epoc of interest
    
    
    %%
    % Applying a time filter to a uniformly sampled signal means that the
    % length of each segment could vary by one sample.  Let's find the minimum
    % length so we can trim the excess off before calculating the mean.
    minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
    minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
    data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), ...
        data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
    data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), ...
        data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
    
    allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');
    
    % downsample 10x and average 405 signal
    N = 10;
    F405 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength1 = size(F405,2);   
    
    % downsample 10x and average 465 signal
    allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
    F465 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength2 = size(F465,2);
    
    %% Plot Epoch Averaged Response
    
    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
    ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
    
    bls = polyfit(F465(1:end), F405(1:end), 1);
    Y_fit_all = bls(1) .* F405 + bls(2);
    Y_dF_all = F465 - Y_fit_all;
    
    zall = zeros(size(Y_dF_all));
    for i = 1:size(Y_dF_all,1)
        ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
        zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
        zsd = std(Y_dF_all(i,ind)); % baseline period stdev
        zall(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
    end
    
    % Smoothes the z score signal traces
    zallSmooth = zeros(size(zall));
    for ii = 1:height(zall)
        zallSmooth(ii,:) = smoothdata(zall(ii,:),'movmean',smoothFactor);
    end
    
    
    % Baseline correction for smoothed data
    if setBaseline == 1
        idx = find(ts1>baseAdjust,1);
        for base = 1:height(zall)
            if zallSmooth(base,idx) < 0
                val = zallSmooth(base,idx);
                diff = 0 - val;
                zallSmooth(base,:) = zallSmooth(base,:) + abs(diff);
            elseif zallSmooth(base,idx) > 0
                val = zallSmooth(base,idx);
                diff = 0 - val;
                zallSmooth(base,:) = zallSmooth(base,:) - abs(diff);
            end
        end
    else
        disp('')
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%LEVERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if dataType == 1
        data = TDTbin2mat(filename, 'TYPE', {'epocs', 'streams'});% TDT function for extracting data to struct 'data'
    elseif dataType == 2
        load(filename)
    end
    
    if isfield(data.epocs, "events")
        data.epocs = rmfield(data.epocs, 'events');
    end
    
    % Use TDTfilter to extract data around our epoc event
    % Using the 'TIME' parameter extracts data only from the time range around
    % our epoc event. Use the 'VALUES' parameter to specify allowed values of
    % the REF_EPOC to extract.  For stream events, the chunks of data are 
    % stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
    data = TDTfilter(data, REF_EPOC2, 'TIME', TRANGE); % extracts data around epoc of interest 
    [~,name,~] = fileparts(filename); % gets name of tank
    brokenID = strsplit(name,'_'); % splits tank name into parts separated by '_'
    % Edit below to have ID and task included in figure title
    
    ID = brokenID(1); % integer following brokenID can be changed depending on what position the ID is in the file name
    task = brokenID(2); % interger following brokenID can be changed depending on what position the task is in the file name
    
    %%
    % Applying a time filter to a uniformly sampled signal means that the
    % length of each segment could vary by one sample.  Let's find the minimum
    % length so we can trim the excess off before calculating the mean.
    minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
    minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
    data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), ...
        data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
    data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), ...
        data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
    
    allSignals2 = cell2mat(data.streams.(STREAM_STORE1).filtered');
    
    % downsample 10x and average 405 signal
    N = 10;
    F405 = zeros(size(allSignals2(:,1:N:end-N+1)));
    for ii = 1:size(allSignals2,1)
        F405(ii,:) = arrayfun(@(i) mean(allSignals2(ii,i:i+N-1)),1:N:length(allSignals2)-N+1);
    end
    minLength1 = size(F405,2);
    
    % downsample 10x and average 465 signal
    allSignals2 = cell2mat(data.streams.(STREAM_STORE2).filtered');
    F465 = zeros(size(allSignals2(:,1:N:end-N+1)));
    for ii = 1:size(allSignals2,1)
        F465(ii,:) = arrayfun(@(i) mean(allSignals2(ii,i:i+N-1)),1:N:length(allSignals2)-N+1);
    end
    minLength2 = size(F465,2);
    
    %% Plot Epoch Averaged Response
    
    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
    ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
    
    bls = polyfit(F465(1:end), F405(1:end), 1);
    Y_fit_all = bls(1) .* F405 + bls(2);
    Y_dF_all = F465 - Y_fit_all;
    
    zall2 = zeros(size(Y_dF_all));
    for i = 1:size(Y_dF_all,1)
        ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
        zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
        zsd = std(Y_dF_all(i,ind)); % baseline period stdev
        zall2(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
    end
    
    % Smoothes the z score signal traces
    zallSmooth2 = zeros(size(zall2));
    for ii = 1:height(zall2)
        zallSmooth2(ii,:) = smoothdata(zall2(ii,:),'movmean',smoothFactor);
    end
    
    
    % Baseline correction for smoothed data
    if setBaseline == 1
        idx = find(ts1>baseAdjust,1);
        for base = 1:height(zall2)
            if zallSmooth2(base,idx) < 0
                val = zallSmooth2(base,idx);
                diff = 0 - val;
                zallSmooth2(base,:) = zallSmooth2(base,:) + abs(diff);
            elseif zallSmooth2(base,idx) > 0
                val = zallSmooth2(base,idx);
                diff = 0 - val;
                zallSmooth2(base,:) = zallSmooth2(base,:) - abs(diff);
            end
        end
    else
        disp('')
    end
   

    allCues{k,1} = zallSmooth;
    allCues{k,2} = name;
    allLevers{k,1} = zallSmooth2;
    allLevers{k,2} = name;
end

signals.cues = allCues;
signals.levers = allLevers;
signals.meta.files = dataFiles;
signals.meta.TRANGE = TRANGE;
signals.meta.BASELINE_PER = BASELINE_PER;
signals.meta.smoothFactor = smoothFactor;
toc
disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);

clearvars -except signals