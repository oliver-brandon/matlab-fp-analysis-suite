%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_Epoc_Extract.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)



%Extracts stream snippets around an epoc (TTL) from TDT fiber photometry
%tanks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = 1.0;
batch_analyze = 2;%1 = batch of tanks in folder, 2 = single tank analysis
EPOC = 'runStop'; % Stimulation event to center on
DLS_ISOS = 'x405A'; % name of the 405 store
DLS_GRABDA = 'x465A'; % name of the 465 store
NAc_ISOS = 'x405C'; % name of the 405 store
NAc_GRABDA = 'x465C'; % name of the 465 store
TRANGE = [-2 12]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-5 -1]; % baseline period before stim
ARTIFACT405 = Inf;
ARTIFACT465 = Inf;
fprintf("VERSION: %d\n",VERSION)
if batch_analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~ismember({myFiles.name},{'.','..','.DS_Store'}));
    numFiles = length(myFiles);
    fprintf("Starting batch extraction of %d files...\n",numFiles)
    for i = 1:numFiles
        BLOCKPATH = fullfile(myDir, myFiles(i).name);
        fprintf("Extracting tank %d of %d...\n",i,numFiles)
        data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams','epocs'});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% Running Epoc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% makes new wheel data epocs %%
        %% epocs created in TDT OpenScope save in the notes of Cam1 %% 
        
        %combine index with timestamp data from Cam1 notes%
        ind = double(data.epocs.Cam1.notes.index);
        ts = data.epocs.Cam1.notes.ts;
        var1 = [ind ts];
        
        %separate by index to make separate epocs%
        % onWheel = var1(ismember(var1(:,1),[1]),:);
        % offWheel = var1(ismember(var1(:,1),[2]),:);
        runStart = var1(ismember(var1(:,1),[3]),:);
        runStop = var1(ismember(var1(:,1),[4]),:);
        
        %extract time stamps%
        % onWheelTs = onWheel(:,2);
        % offWheelTs = offWheel(:,2);
        runStartTs = runStart(:,2);
        runStopTs = runStop(:,2);
        %extract indicies%
        % onWheelInd = onWheel(:,1);
        % offWheelInd = offWheel(:,1);
        runStartInd = runStart(:,1);
        runStopInd = runStop(:,1);
        %make a new epoc structure based on Cam1 notes extracted data%
        %onWheel%
        % data.epocs.onWheel.onset = onWheelTs;
        % data.epocs.onWheel.offset = onWheelTs + 0.01;
        % data.epocs.onWheel.name = 'onWheel';
        % data.epocs.onWheel.data = onWheelInd;
        %offWheel%
        % data.epocs.offWheel.onset = offWheelTs;
        % data.epocs.offWheel.offset = offWheelTs + 0.01;
        % data.epocs.offWheel.name = 'offWheel';
        % data.epocs.offWheel.data = offWheelInd;
        %runStart%
        data.epocs.runStart.onset = runStartTs;
        data.epocs.runStart.offset = runStartTs + 0.01;
        data.epocs.runStart.name = 'runStart';
        data.epocs.runStart.data = runStartInd;
        %runStop%
        data.epocs.runStop.onset = runStopTs;
        data.epocs.runStop.offset = runStopTs + 0.01;
        data.epocs.runStop.name = 'runStop';
        data.epocs.runStop.data = runStopInd;
        REF_EPOC = EPOC;
        data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%   DLS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
        % below -ARTIFACT level, remove it from the data set.
        art1_DLS = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false));
        art2_DLS = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false));
        good_DLS = ~art1_DLS & ~art2_DLS;
        data.streams.(DLS_ISOS).filtered = data.streams.(DLS_ISOS).filtered(good_DLS);
        art1_DLS = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false));
        art2_DLS = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false));
        good2_DLS = ~art1_DLS & ~art2_DLS;
        data.streams.(DLS_GRABDA).filtered = data.streams.(DLS_GRABDA).filtered(good2_DLS);
        numArtifacts_DLS = sum(~good_DLS) + sum(~good2_DLS);
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.
        minLength1 = min(cellfun('prodofsize', data.streams.(DLS_ISOS).filtered));
        minLength2 = min(cellfun('prodofsize', data.streams.(DLS_GRABDA).filtered));
        data.streams.(DLS_ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false);
        data.streams.(DLS_GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false);
        allSignals_DLS = cell2mat(data.streams.(DLS_ISOS).filtered');
        % downsample 10x and average 405 signal
        N = 10;
        F405_DLS = zeros(size(allSignals_DLS(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_DLS,1)
            F405_DLS(ii,:) = arrayfun(@(i) mean(allSignals_DLS(ii,i:i+N-1)),1:N:length(allSignals_DLS)-N+1);
        end
        minLength1_DLS = size(F405_DLS,2);
        % Create mean signal, standard error of signal, and DC offset of 405 signal
        meanSignal1_DLS = mean(F405_DLS);
        stdSignal1_DLS = std(double(F405_DLS))/sqrt(size(F405_DLS,1));
        dcSignal1_DLS = mean(meanSignal1_DLS);
        
        % downsample 10x and average 465 signal
        allSignals_DLS = cell2mat(data.streams.(DLS_GRABDA).filtered');
        F465_DLS = zeros(size(allSignals_DLS(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_DLS,1)
            F465_DLS(ii,:) = arrayfun(@(i) mean(allSignals_DLS(ii,i:i+N-1)),1:N:length(allSignals_DLS)-N+1);
        end
        minLength2_DLS = size(F465_DLS,2);
        
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2_DLS = mean(F465_DLS);
        stdSignal2_DLS = std(double(F465_DLS))/sqrt(size(F465_DLS,1));
        dcSignal2_DLS = mean(meanSignal2_DLS);
        % Create the time vector for each stream store
        ts1_DLS = TRANGE(1) + (1:minLength1_DLS) / data.streams.(DLS_ISOS).fs*N;
        ts2_DLS = TRANGE(1) + (1:minLength2_DLS) / data.streams.(DLS_GRABDA).fs*N;

        meanSignal1_DLS = meanSignal1_DLS - dcSignal1_DLS;
        meanSignal2_DLS = meanSignal2_DLS - dcSignal2_DLS;
        
        DLS_epoc_store(i,:) = meanSignal2_DLS;
        DLS_epoc_amp(i,:) = max(meanSignal2_DLS(1,204:713));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%   NAc   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
        % below -ARTIFACT level, remove it from the data set.
        art1_NAc = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false));
        art2_NAc = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false));
        good_NAc = ~art1_NAc & ~art2_NAc;
        data.streams.(NAc_ISOS).filtered = data.streams.(NAc_ISOS).filtered(good_DLS);
        art1_NAc = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false));
        art2_NAc = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false));
        good2_NAc = ~art1_NAc & ~art2_NAc;
        data.streams.(NAc_GRABDA).filtered = data.streams.(NAc_GRABDA).filtered(good2_NAc);
        numArtifacts_NAc = sum(~good_NAc) + sum(~good2_NAc);
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.
        minLength1 = min(cellfun('prodofsize', data.streams.(NAc_ISOS).filtered));
        minLength2 = min(cellfun('prodofsize', data.streams.(NAc_GRABDA).filtered));
        data.streams.(NAc_ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false);
        data.streams.(NAc_GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false);
        allSignals_NAc = cell2mat(data.streams.(NAc_ISOS).filtered');
        % downsample 10x and average 405 signal
        N = 10;
        F405_NAc = zeros(size(allSignals_NAc(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_NAc,1)
            F405_NAc(ii,:) = arrayfun(@(i) mean(allSignals_NAc(ii,i:i+N-1)),1:N:length(allSignals_NAc)-N+1);
        end
        minLength1_NAc = size(F405_NAc,2);
        % Create mean signal, standard error of signal, and DC offset of 405 signal
        meanSignal1_NAc = mean(F405_NAc);
        stdSignal1_NAc = std(double(F405_NAc))/sqrt(size(F405_NAc,1));
        dcSignal1_NAc = mean(meanSignal1_NAc);
        
        % downsample 10x and average 465 signal
        allSignals_NAc = cell2mat(data.streams.(NAc_GRABDA).filtered');
        F465_NAc = zeros(size(allSignals_NAc(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_NAc,1)
            F465_NAc(ii,:) = arrayfun(@(i) mean(allSignals_NAc(ii,i:i+N-1)),1:N:length(allSignals_NAc)-N+1);
        end
        minLength2_NAc = size(F465_NAc,2);
        
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2_NAc = mean(F465_NAc);
        stdSignal2_NAc = std(double(F465_NAc))/sqrt(size(F465_NAc,1));
        dcSignal2_NAc = mean(meanSignal2_NAc);
        
        % Create the time vector for each stream store
        ts1_NAc = TRANGE(1) + (1:minLength1_NAc) / data.streams.(NAc_ISOS).fs*N;
        ts2_NAc = TRANGE(1) + (1:minLength2_NAc) / data.streams.(NAc_GRABDA).fs*N;
        
        meanSignal1_NAc = meanSignal1_NAc - dcSignal1_NAc;
        meanSignal2_NAc = meanSignal2_NAc - dcSignal2_NAc;

        NAc_epoc_store(i,:) = meanSignal2_NAc;
        NAc_epoc_amp(i,:) = max(meanSignal2_NAc(1,204:713));
    end
elseif batch_analyze == 2
        disp("Starting single tank extraction...")
        myDir = uigetdir;
        numFiles = 1;
        BLOCKPATH = myDir;
        disp("Extracting...")
        data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams','epocs'});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% Running Epoc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% makes new wheel data epocs %%
        %% epocs created in TDT OpenScope save in the notes of Cam1 %% 
        
        %combine index with timestamp data from Cam1 notes%
        ind = double(data.epocs.Cam1.notes.index);
        ts = data.epocs.Cam1.notes.ts;
        var1 = [ind ts];
        
        %separate by index to make separate epocs%
        % onWheel = var1(ismember(var1(:,1),[1]),:);
        % offWheel = var1(ismember(var1(:,1),[2]),:);
        runStart = var1(ismember(var1(:,1),[3]),:);
        runStop = var1(ismember(var1(:,1),[4]),:);
        
        %extract time stamps%
        % onWheelTs = onWheel(:,2);
        % offWheelTs = offWheel(:,2);
        runStartTs = runStart(:,2);
        runStopTs = runStop(:,2);
        %extract indicies%
        % onWheelInd = onWheel(:,1);
        % offWheelInd = offWheel(:,1);
        runStartInd = runStart(:,1);
        runStopInd = runStop(:,1);
        %make a new epoc structure based on Cam1 notes extracted data%
        %onWheel%
        % data.epocs.onWheel.onset = onWheelTs;
        % data.epocs.onWheel.offset = onWheelTs + 0.01;
        % data.epocs.onWheel.name = 'onWheel';
        % data.epocs.onWheel.data = onWheelInd;
        %offWheel%
        % data.epocs.offWheel.onset = offWheelTs;
        % data.epocs.offWheel.offset = offWheelTs + 0.01;
        % data.epocs.offWheel.name = 'offWheel';
        % data.epocs.offWheel.data = offWheelInd;
        %runStart%
        data.epocs.runStart.onset = runStartTs;
        data.epocs.runStart.offset = runStartTs + 0.01;
        data.epocs.runStart.name = 'runStart';
        data.epocs.runStart.data = runStartInd;
        %runStop%
        data.epocs.runStop.onset = runStopTs;
        data.epocs.runStop.offset = runStopTs + 0.01;
        data.epocs.runStop.name = 'runStop';
        data.epocs.runStop.data = runStopInd;
        REF_EPOC = EPOC;
        data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%   DLS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
        % below -ARTIFACT level, remove it from the data set.
        art1_DLS = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false));
        art2_DLS = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false));
        good_DLS = ~art1_DLS & ~art2_DLS;
        data.streams.(DLS_ISOS).filtered = data.streams.(DLS_ISOS).filtered(good_DLS);
        art1_DLS = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false));
        art2_DLS = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false));
        good2_DLS = ~art1_DLS & ~art2_DLS;
        data.streams.(DLS_GRABDA).filtered = data.streams.(DLS_GRABDA).filtered(good2_DLS);
        numArtifacts_DLS = sum(~good_DLS) + sum(~good2_DLS);
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.
        minLength1 = min(cellfun('prodofsize', data.streams.(DLS_ISOS).filtered));
        minLength2 = min(cellfun('prodofsize', data.streams.(DLS_GRABDA).filtered));
        data.streams.(DLS_ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(DLS_ISOS).filtered, 'UniformOutput',false);
        data.streams.(DLS_GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(DLS_GRABDA).filtered, 'UniformOutput',false);
        allSignals_DLS = cell2mat(data.streams.(DLS_ISOS).filtered');
        % downsample 10x and average 405 signal
        N = 10;
        F405_DLS = zeros(size(allSignals_DLS(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_DLS,1)
            F405_DLS(ii,:) = arrayfun(@(i) mean(allSignals_DLS(ii,i:i+N-1)),1:N:length(allSignals_DLS)-N+1);
        end
        minLength1_DLS = size(F405_DLS,2);
        % Create mean signal, standard error of signal, and DC offset of 405 signal
        meanSignal1_DLS = mean(F405_DLS);
        stdSignal1_DLS = std(double(F405_DLS))/sqrt(size(F405_DLS,1));
        dcSignal1_DLS = mean(meanSignal1_DLS);
        
        % downsample 10x and average 465 signal
        allSignals_DLS = cell2mat(data.streams.(DLS_GRABDA).filtered');
        F465_DLS = zeros(size(allSignals_DLS(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_DLS,1)
            F465_DLS(ii,:) = arrayfun(@(i) mean(allSignals_DLS(ii,i:i+N-1)),1:N:length(allSignals_DLS)-N+1);
        end
        minLength2_DLS = size(F465_DLS,2);
        
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2_DLS = mean(F465_DLS);
        stdSignal2_DLS = std(double(F465_DLS))/sqrt(size(F465_DLS,1));
        dcSignal2_DLS = mean(meanSignal2_DLS);
        % Create the time vector for each stream store
        ts1_DLS = TRANGE(1) + (1:minLength1_DLS) / data.streams.(DLS_ISOS).fs*N;
        ts2_DLS = TRANGE(1) + (1:minLength2_DLS) / data.streams.(DLS_GRABDA).fs*N;

        meanSignal1_DLS = meanSignal1_DLS - dcSignal1_DLS;
        meanSignal2_DLS = meanSignal2_DLS - dcSignal2_DLS;
    
        DLS_epoc_store = meanSignal2_DLS;
        DLS_epoc_amp = max(meanSignal2_DLS(1,204:713));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%   NAc   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
        % below -ARTIFACT level, remove it from the data set.
        art1_NAc = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false));
        art2_NAc = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false));
        good_NAc = ~art1_NAc & ~art2_NAc;
        data.streams.(NAc_ISOS).filtered = data.streams.(NAc_ISOS).filtered(good_DLS);
        art1_NAc = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false));
        art2_NAc = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false));
        good2_NAc = ~art1_NAc & ~art2_NAc;
        data.streams.(NAc_GRABDA).filtered = data.streams.(NAc_GRABDA).filtered(good2_NAc);
        numArtifacts_NAc = sum(~good_NAc) + sum(~good2_NAc);
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.
        minLength1 = min(cellfun('prodofsize', data.streams.(NAc_ISOS).filtered));
        minLength2 = min(cellfun('prodofsize', data.streams.(NAc_GRABDA).filtered));
        data.streams.(NAc_ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(NAc_ISOS).filtered, 'UniformOutput',false);
        data.streams.(NAc_GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(NAc_GRABDA).filtered, 'UniformOutput',false);
        allSignals_NAc = cell2mat(data.streams.(NAc_ISOS).filtered');
        % downsample 10x and average 405 signal
        N = 10;
        F405_NAc = zeros(size(allSignals_NAc(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_NAc,1)
            F405_NAc(ii,:) = arrayfun(@(i) mean(allSignals_NAc(ii,i:i+N-1)),1:N:length(allSignals_NAc)-N+1);
        end
        minLength1_NAc = size(F405_NAc,2);
        % Create mean signal, standard error of signal, and DC offset of 405 signal
        meanSignal1_NAc = mean(F405_NAc);
        stdSignal1_NAc = std(double(F405_NAc))/sqrt(size(F405_NAc,1));
        dcSignal1_NAc = mean(meanSignal1_NAc);
        
        % downsample 10x and average 465 signal
        allSignals_NAc = cell2mat(data.streams.(NAc_GRABDA).filtered');
        F465_NAc = zeros(size(allSignals_NAc(:,1:N:end-N+1)));
        for ii = 1:size(allSignals_NAc,1)
            F465_NAc(ii,:) = arrayfun(@(i) mean(allSignals_NAc(ii,i:i+N-1)),1:N:length(allSignals_NAc)-N+1);
        end
        minLength2_NAc = size(F465_NAc,2);
        
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2_NAc = mean(F465_NAc);
        stdSignal2_NAc = std(double(F465_NAc))/sqrt(size(F465_NAc,1));
        dcSignal2_NAc = mean(meanSignal2_NAc);
        
        % Create the time vector for each stream store
        ts1_NAc = TRANGE(1) + (1:minLength1_NAc) / data.streams.(NAc_ISOS).fs*N;
        ts2_NAc = TRANGE(1) + (1:minLength2_NAc) / data.streams.(NAc_GRABDA).fs*N;
        
        meanSignal1_NAc = meanSignal1_NAc - dcSignal1_NAc;
        meanSignal2_NAc = meanSignal2_NAc - dcSignal2_NAc;

        NAc_epoc_store = meanSignal2_NAc;
        NAc_epoc_amp = max(meanSignal2_NAc(1,204:713));
    
end
fprintf("Successfully extracted streams from %d tank(s)\n",numFiles)
disp("The DLS signal(s) is stored in DLS_epoc_store")
disp("The NAc signal(s) is stored in NAc_epoc_store")