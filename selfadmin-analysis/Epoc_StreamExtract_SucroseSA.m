
%% Housekeeping
% Clear workspace and close existing figures.
% Add data
clear all; clc; close all;
tanks2analyze = 1; % 1=batch, 2=single
streamAorC = 2; % 1=465A, 2=465C
epoc = {'aRL/','bRL/'}; % enter the epoc of interest (a and b)




if tanks2analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
    myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv'}));
    numFiles = length(myFiles);
    epoc_stream_store = cell(numFiles);
    disp("Starting batch tank extract...")
elseif tanks2analyze == 2
    myDir = uigetdir;
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
    animalIDA = char(brokenID{2});
    if streamAorC == 1
        emptylogicA = strcmp(animalIDA,emptyID);
        if emptylogicA == 1
            disp("Stream A is empty")
            continue
        elseif emptylogicA == 0
            disp('')
        end
    end
    animalIDC = char(brokenID{5});
    if streamAorC == 2
        emptylogicC = strcmp(animalIDC,emptyID);
        if emptylogicC == 1
            disp("Stream C is empty")
            continue
        elseif emptylogicC == 0
            disp('')
        end
    end
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs','streams'});
    REF_EPOC = char(epoc(streamAorC));
    TRANGE = [-2 12]; %window size [start time relative to epoc onset, entire duration]
    BASELINE_PER = [-5 -1]; % baseline period before stim
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
    if tanks2analyze == 1
        epoc_stream_store{batch} = F465;
    elseif tanks2analyze == 2
        epoc_stream_store = F465;
    end
       
end
fprintf("Finished extracting epocs from %d tanks...\n",numFiles)
