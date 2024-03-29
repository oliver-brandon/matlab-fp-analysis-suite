
%% Housekeeping
% Clear workspace and close existing figures.
% Add data
clear all; clc; close all;
FR = 1; % Fixed ratio
tanks2analyze = 2; % 1=batch, 2=single
streamAorC = 1; % 1=465A, 2=465C
% epoc = {'aHL/','bHL/'};

remove_infusion_np = 0; % 1 = remove infusion nosepokes
remove_free_inf = 1; % 1 = remove nosepokes from free infusion at start
% epoc = {'aReward','bReward'}; % enter the epoc of interest (a and b)
% epoc = {'aActiveRew','bActiveRew'};
% epoc = {'aActiveTimeout','bActiveTimeout'};
epoc = {'aRL/','bRL/'};
TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-3 -1]; % baseline period before stim

% epoc = {'aRw/','bRw/'};
% TRANGE = [-2 22]; %window size [start time relative to epoc onset, entire duration]
% BASELINE_PER = [-2 0]; % baseline period before stim

if tanks2analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
    myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv'}));
    numFiles = length(myFiles);
    epoc_stream_store = cell(numFiles);
    disp("Starting batch tank extract...")
elseif tanks2analyze == 2
    myDir = uigetdir('/Users/brandon/personal-drive/self_admin/coc_sa/PrL-aIC/FR1-Cocaine/tanks');
    numFiles = 1;
    epoc_stream_store = cell(numFiles);
    disp("Starting single tank extract...")
end
for batch = 1:numFiles
    fprintf('Loading tank %d of %d...\n',batch,numFiles)
    if tanks2analyze == 1
        BLOCKPATH = fullfile(myDir, myFiles(batch).name);
    elseif tanks2analyze == 2
        BLOCKPATH = myDir;
    end 
    [~,name,~] = fileparts(BLOCKPATH);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    animalIDA = char(brokenID{1});
    animalIDC = char(brokenID{3});
    taskA = char(brokenID{2});
    taskC = char(brokenID{4});
    if streamAorC == 1
        emptylogicA = strcmp(animalIDA,emptyID);
        if emptylogicA == 1
            disp("Stream A is empty")
            continue
        elseif emptylogicA == 0
            data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs','streams'});
            [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(data.epocs.aRL_.onset, 20);
            [data] = createEpoc(data, rewardTimestamps, 'aActiveRew');
            [data] = createEpoc(data, rewardTimeout, 'aRewTimeout');
            [data] = createEpoc(data, timeoutTimestamps, 'aActiveTimeout');
            [data] = createEpoc(data, data.epocs.aRw_.offset, 'aReward');
        end
    end
    
    if streamAorC == 2
        emptylogicC = strcmp(animalIDC,emptyID);
        if emptylogicC == 1
            disp("Stream C is empty")
            continue
        elseif emptylogicC == 0
            data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs','streams'});
            [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(data.epocs.bRL_.onset, 20);
            [data] = createEpoc(data, rewardTimestamps, 'bActiveRew');
            [data] = createEpoc(data, rewardTimeout, 'bRewTimeout');
            [data] = createEpoc(data, timeoutTimestamps, 'bActiveTimeout');
            [data] = createEpoc(data, data.epocs.bRw_.offset, 'bReward');
        end
    end
    
    % Creates reward epoc from offset instead of onset
    % data.epocs.aReward.onset = data.epocs.aRw_.offset;
    % data.epocs.aReward.offset = data.epocs.aRw_.onset;
    % data.epocs.aReward.name = 'aReward';
    % data.epocs.aReward.data = ones(height(data.epocs.aRw_.offset)) * 10;
    % data.epocs.bReward.onset = data.epocs.bRw_.offset;
    % data.epocs.bReward.offset = data.epocs.bRw_.onset;
    % data.epocs.bReward.name = 'bReward';
    % data.epocs.bReward.data = ones(height(data.epocs.bRw_.offset)) * 20;
    REF_EPOC = char(epoc(streamAorC));
    % Removes the nosepokes that are from the simulated responses (free
    % infusion)
    if remove_free_inf == 1
        if streamAorC == 1
            data.epocs.aRL_.onset = data.epocs.aRL_.onset(FR+1:end,:);
            data.epocs.aActiveRew.onset = data.epocs.aActiveRew.onset(FR+1:end,:);
            data.epocs.aActiveTimeout.onset = data.epocs.aActiveTimeout.onset(FR+1:end,:);
        elseif streamAorC == 2
            data.epocs.bRL_.onset = data.epocs.bRL_.onset(FR+1:end,:);
            data.epocs.bActiveRew.onset = data.epocs.bActiveRew.onset(FR+1:end,:);
            data.epocs.bActiveTimeout.onset = data.epocs.bActiveTimeout.onset(FR+1:end,:);
        end
    end
    % Removes the nosepoke epoc that results in an infusion 
    if remove_infusion_np == 1
        if streamAorC == 1
            indices_to_remove = 1:FR:size(data.epocs.aRL_.onset,1);
            data.epocs.aRL_.onset = data.epocs.aRL_.onset(setdiff(1:end,indices_to_remove),:);
        elseif streamAorC == 2
            indices_to_remove = 1:FR:size(data.epocs.bRL_.onset,1);
            data.epocs.bRL_.onset = data.epocs.bRL_.onset(setdiff(1:end,indices_to_remove),:);
        end
    end

    
    
    ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
    ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
    if streamAorC == 1
        disp("Extracting stream A")
        STREAM_STORE1 = 'x405A'; % name of the 405 store
        STREAM_STORE2 = 'x465A'; % name of the 465 store
    elseif streamAorC == 2
        disp("Extracting stream C")
        STREAM_STORE1 = 'x405C'; % name of the 405 store
        STREAM_STORE2 = 'x465C'; % name of the 465 store
    end
    % Use TDTfilter to extract data around our epoc event
    % Using the 'TIME' parameter extracts data only from the time range around
    % our epoc event. Use the 'VALUES' parameter to specify allowed values of
    % the REF_EPOC to extract.  For stream events, the chunks of data are 
    % stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
    data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);
    % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
    % below -ARTIFACT level, remove it from the data set.
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
    good = ~art1 & ~art2;
    data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);
    
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
    good2 = ~art1 & ~art2;
    data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);
    
    numArtifacts = sum(~good) + sum(~good2);
    
    %%
    % Applying a time filter to a uniformly sampled signal means that the
    % length of each segment could vary by one sample.  Let's find the minimum
    % length so we can trim the excess off before calculating the mean.
    minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
    minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
    data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
    data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
    
    allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');
    
    % downsample 10x and average 405 signal
    N = 10;
    F405 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength1 = size(F405,2);
    
    % Create mean signal, standard error of signal, and DC offset of 405 signal
    meanSignal1 = mean(F405);
    stdSignal1 = std(double(F405))/sqrt(size(F405,1));
    dcSignal1 = mean(meanSignal1);
    
    % downsample 10x and average 465 signal
    allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
    F465 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength2 = size(F465,2);
    
    % Create mean signal, standard error of signal, and DC offset of 465 signal
    meanSignal2 = mean(F465);
    stdSignal2 = std(double(F465))/sqrt(size(F465,1));
    dcSignal2 = mean(meanSignal2);
    
    %% Plot Epoch Averaged Response
    
    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
    ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
    
    % Subtract DC offset to get signals on top of one another
    meanSignal1 = meanSignal1 - dcSignal1;
    meanSignal2 = meanSignal2 - dcSignal2;

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

    if tanks2analyze == 1
        epoc_stream_store{batch} = zall;
    elseif tanks2analyze == 2
        epoc_stream_store = zall;
        f1 = figure;
        imagesc(ts1, 1, zall);
        colormap('jet'); colorbar; 
        title(sprintf('Z-Score/Trial', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 20);
        xlabel('Time, s', 'FontSize', 16);
        ylabel('Trial', 'FontSize', 16); 
    end
  
end



fprintf("Finished extracting epocs from %d tanks...\n",numFiles)
