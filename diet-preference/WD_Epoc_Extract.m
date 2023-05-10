%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%WD_Epoc_Extract.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)

VERSION = "1.0";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
disp(VERSION)

batch_analyze = 1;%1 = batch of tanks in folder, 2 = single tank analysis
EPOC = 'WDApp'; % Stimulation event to center on (WD_app,WD_awy,SD_app,SD_awy)
arena_number = 1;% arena 1 = first ID in tank name, arena 2 = second ID
if arena_number == 1
        ISOS = 'x405A'; % name of the 405 store
        GRABDA = 'x465A'; % name of the 465 store
elseif arena_number == 2
        ISOS = 'x405C'; % name of the 405 store
        GRABDA = 'x465C'; % name of the 465 store
end
TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-5 -1]; % baseline period before stim
ARTIFACT405 = Inf;
ARTIFACT465 = Inf;
disp(VERSION)
if batch_analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~ismember({myFiles.name},{'.','..','.DS_Store','*.xlsx'}));
    numFiles = length(myFiles);
    fprintf("Starting batch extraction of %d files...\n",numFiles)
    for i = 1:numFiles
        BLOCKPATH = fullfile(myDir, myFiles(i).name);
        fprintf("Extracting tank %d of %d...\n",i,numFiles)
        data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams','epocs'});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% Create Epocs %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%%%% SD Approach and Away %%%%%
        data.epocs.SDApp.onset = exceldata(:, 1);
        data.epocs.SDApp.offset = data.epocs.SDApp.onset + 2;
        data.epocs.SDApp.name = 'SDApp';
        data.epocs.SDApp.data = ones(length(data.epocs.SDApp.onset), 1);
        data.epocs.SDAwy.onset = exceldata(:, 2);
        data.epocs.SDAwy.offset = data.epocs.SDAwy.onset + 2;
        data.epocs.SDAwy.name = 'SDAwy';
        data.epocs.SDAwy.data = ones(length(data.epocs.SDAwy.onset), 1) * 2;
        %%%%% WD Approach and Away %%%%%
        data.epocs.WDApp.onset = exceldata(:, 3);
        data.epocs.WDApp.offset = data.epocs.WDApp.onset + 2;
        data.epocs.WDApp.name = 'WDApp';
        data.epocs.WDApp.data = ones(length(data.epocs.WDApp.onset), 1) * 3;
        data.epocs.WDAwy.onset = exceldata(:, 4);
        data.epocs.WDAwy.offset = data.epocs.WDAwy.onset + 2;
        data.epocs.WDAwy.name = 'WDAwy';
        data.epocs.WDAwy.data = ones(length(data.epocs.WDAwy.onset), 1) * 4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%   Signal Processing   %%%%%%%%%%%%%%%%%%%%
        % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
        % below -ARTIFACT level, remove it from the data set.
        art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(ISOS).filtered, 'UniformOutput',false));
        art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(ISOS).filtered, 'UniformOutput',false));
        good = ~art1 & ~art2;
        data.streams.(ISOS).filtered = data.streams.(ISOS).filtered(good);
        art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(GRABDA).filtered, 'UniformOutput',false));
        art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(GRABDA).filtered, 'UniformOutput',false));
        good2 = ~art1 & ~art2;
        data.streams.(GRABDA).filtered = data.streams.(GRABDA).filtered(good2);
        numArtifacts_DLS = sum(~good) + sum(~good2);
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.
        minLength1 = min(cellfun('prodofsize', data.streams.(ISOS).filtered));
        minLength2 = min(cellfun('prodofsize', data.streams.(GRABDA).filtered));
        data.streams.(ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(ISOS).filtered, 'UniformOutput',false);
        data.streams.(GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(GRABDA).filtered, 'UniformOutput',false);
        allSignals = cell2mat(data.streams.(ISOS).filtered');
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
        allSignals = cell2mat(data.streams.(GRABDA).filtered');
        F465 = zeros(size(allSignals(:,1:N:end-N+1)));
        for ii = 1:size(allSignals,1)
            F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
        end
        minLength2_DLS = size(F465,2);
        
        % Create mean signal, standard error of signal, and DC offset of 465 signal
        meanSignal2 = mean(F465);
        stdSignal2 = std(double(F465))/sqrt(size(F465,1));
        dcSignal2 = mean(meanSignal2);
        % Create the time vector for each stream store
        ts1 = TRANGE(1) + (1:minLength1) / data.streams.(ISOS).fs*N;
        ts2 = TRANGE(1) + (1:minLength2_DLS) / data.streams.(GRABDA).fs*N;

        meanSignal1 = meanSignal1 - dcSignal1;
        meanSignal2 = meanSignal2 - dcSignal2;
        
        epoc_stream_store(i,:) = meanSignal2;
        epoc_amp(i,:) = max(meanSignal2(1,204:713));
        
    end
elseif batch_analyze == 2
    disp("Starting single tank extraction...")
    myDir = uigetdir;
    numFiles = 1;
    BLOCKPATH = myDir;
    disp("Extracting...")
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams','epocs'});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Create Epocs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %combine index with timestamp data from Cam1 notes%
    ind = double(data.epocs.Cam1.notes.index);
    ts = data.epocs.Cam1.notes.ts;
    var1 = [ind ts];
    %separate by index to make separate epocs%
    SD_app = var1(ismember(var1(:,1),[1]),:);
    SD_awy = var1(ismember(var1(:,1),[2]),:);
    WD_app = var1(ismember(var1(:,1),[3]),:);
    WD_awy = var1(ismember(var1(:,1),[4]),:);
    %extract time stamps%
    SD_app_ts = SD_app(:,2);
    SD_awy_ts = SD_awy(:,2);
    WD_app_ts = WD_app(:,2);
    WD_awy_ts = WD_awy(:,2);
    %extract indicies%
    SD_app_ind = SD_app(:,1);
    SD_awy_ind = SD_awy(:,1);
    WD_app_ind = WD_app(:,1);
    WD_awy_ind = WD_awy(:,1);
    %make a new epoc structure based on Cam1 notes extracted data%
    %SD Approach%
    data.epocs.SDapp.onset = SD_app_ts;
    data.epocs.SDapp.offset = SD_app_ts + 0.01;
    data.epocs.SDapp.name = 'SD_app';
    data.epocs.SDapp.data = SD_app_ind;
    %SD Away%
    data.epocs.SDawy.onset = SD_awy_ts;
    data.epocs.SDawy.offset = SD_awy_ts + 0.01;
    data.epocs.SDawy.name = 'SD_awy';
    data.epocs.SDawy.data = SD_awy_ind;
    %WD Approach%
    data.epocs.WDapp.onset = WD_app_ts;
    data.epocs.WDapp.offset = WD_app_ts + 0.01;
    data.epocs.WDapp.name = 'WD_app';
    data.epocs.WDapp.data = WD_app_ind;
    %WD Away%
    data.epocs.WDawy.onset = WD_awy_ts;
    data.epocs.WDawy.offset = WD_awy_ts + 0.01;
    data.epocs.WDawy.name = 'WD_awy';
    data.epocs.WDawy.data = WD_awy_ind;
    REF_EPOC = EPOC;
    data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%   Signal Processing   %%%%%%%%%%%%%%%%%%%%%%%%
    % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
    % below -ARTIFACT level, remove it from the data set.
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(ISOS).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(ISOS).filtered, 'UniformOutput',false));
    good = ~art1 & ~art2;
    data.streams.(ISOS).filtered = data.streams.(ISOS).filtered(good);
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(GRABDA).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(GRABDA).filtered, 'UniformOutput',false));
    good2 = ~art1 & ~art2;
    data.streams.(GRABDA).filtered = data.streams.(GRABDA).filtered(good2);
    numArtifacts_DLS = sum(~good) + sum(~good2);
    % Applying a time filter to a uniformly sampled signal means that the
    % length of each segment could vary by one sample.  Let's find the minimum
    % length so we can trim the excess off before calculating the mean.
    minLength1 = min(cellfun('prodofsize', data.streams.(ISOS).filtered));
    minLength2 = min(cellfun('prodofsize', data.streams.(GRABDA).filtered));
    data.streams.(ISOS).filtered = cellfun(@(x) x(1:minLength1), data.streams.(ISOS).filtered, 'UniformOutput',false);
    data.streams.(GRABDA).filtered = cellfun(@(x) x(1:minLength2), data.streams.(GRABDA).filtered, 'UniformOutput',false);
    allSignals = cell2mat(data.streams.(ISOS).filtered');
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
    allSignals = cell2mat(data.streams.(GRABDA).filtered');
    F465 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength2_DLS = size(F465,2);
    
    % Create mean signal, standard error of signal, and DC offset of 465 signal
    meanSignal2 = mean(F465);
    stdSignal2 = std(double(F465))/sqrt(size(F465,1));
    dcSignal2 = mean(meanSignal2);
    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(ISOS).fs*N;
    ts2 = TRANGE(1) + (1:minLength2_DLS) / data.streams.(GRABDA).fs*N;

    meanSignal1 = meanSignal1 - dcSignal1;
    meanSignal2 = meanSignal2 - dcSignal2;
    
    epoc_stream_store = meanSignal2;
    epoc_amp = max(meanSignal2(1,204:713));
    
    % Plot the 405 and 465 average signals
    figure(1);
    subplot(3,1,1);
    plot(ts1, meanSignal1, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
    plot(ts2, meanSignal2, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); 
    
    % Plot vertical line at epoch onset, time = 0
    line([0 0], [(min(F465(:) - dcSignal2)), ((max(F465(:)) - dcSignal2))], 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
    
    % Create the standard error bands for the 405 signal
    XX = [ts1, fliplr(ts1)];
    YY = [meanSignal1 + stdSignal1, fliplr(meanSignal1 - stdSignal1)];
    
    % Plot filled standard error bands.
    h = fill(XX, YY, 'g');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Repeat for 465
    XX = [ts2, fliplr(ts2)];
    YY = [meanSignal2 + stdSignal2, fliplr(meanSignal2 - stdSignal2)];
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('mV', 'FontSize', 12)
    title(sprintf(TITLE1, numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    set(gcf, 'Position',[100, 100, 800, 500])
    
    % Heat Map based on z score of 405 fit subtracted 465
    % Fitting 405 channel onto 465 channel to detrend signal bleaching
    % Scale and fit data
    % Algorithm sourced from Tom Davidson's Github:
    % https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
    
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
    
    % Standard error of the z-score
    meanZall = mean(zall);
    zerror = std(zall)/sqrt(size(zall,1));
    
    % Plot heat map
    subplot(3,1,2);
    imagesc(ts2, 1, zall);
    colormap('jet'); % c1 = colorbar; 
    title(sprintf('Z-Score Heat Map', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14);
    ylabel('Trials', 'FontSize', 12);
    
    % Fill band values for second subplot. Doing here to scale onset bar
    % correctly
    XX = [ts2, fliplr(ts2)];
    YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
    
    subplot(3,1,3)
    plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('Z-score', 'FontSize', 12)
    title(sprintf('465 nm Z-Score', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    %c2 = colorbar;
    %%
    figure(2)
    plot(ts2, zall)
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('Z-score', 'FontSize', 12)
    title(sprintf('465 nm Z-Score', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    
    %%
    figure(3)
    subplot(2,3,[1,2,4,5])
    plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    h = fill(XX, YY, 'b');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',18)
    ylabel('Z-score +/- SEM', 'FontSize', 18)
    title(sprintf(TITLE1, numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 18)
    box off
    
    subplot(2,3,6);
    imagesc(ts2, 1, zall);
    colormap('jet'); colorbar; 
    title(sprintf('Z-Score/Trial', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 16);
    xlabel('Time, s', 'FontSize', 12);
    ylabel('Trial', 'FontSize', 12);
    
    % Fill band values for second subplot. Doing here to scale onset bar
    % correctly
    XX = [ts2, fliplr(ts2)];
    YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
    
end
fprintf("Successfully extracted streams from %d tank(s)\n",numFiles)
fprintf("The %d epoc signal(s) is stored in epoc_stream_store",REF_EPOC)
