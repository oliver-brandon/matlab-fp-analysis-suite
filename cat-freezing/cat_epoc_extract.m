clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To have proper graph titles, files need to be named ID_Task_ID_Task     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Paramaters to Edit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setBaseline = 1; % 1 = yes, 0 = no (adjusts signals to zero indicated by baseAdjust)
baseAdjust = -2; % seconds on x axis to adjust baseline to
plotSmoothed = 1; % 1 = yes, 0 = no
smoothFactor = 20;
TRANGE = [-2 12]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-10 -1]; % baseline period before epoc
N = 10; % downsample
BLOCKPATH = '/Users/brandon/ucr-drive/collabs/Yaminaka/20250213_FPSF1-05_FPSF1-06'; % path to TDT data tank (folder containing TDT data)
channel = 1; % 1 = mouse on A channel, 2 = mouse on C channel
dataType = 1; % 1 = tank, 2 = .mat

useExcelTS = 1; % 1 = yes, 0 = no
excelPath = '/Users/brandon/ucr-drive/collabs/Yaminaka/FPSF1-05 Behavior Timestamps.xlsx';
excelSheetName = '-05 Control';
excelRange = 'G2:H9';
epocName = 'ControlFreeze';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Leave Code Below As Is %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dataType == 1
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});% TDT function for extracting data to struct 'data'
elseif dataType == 2
    load(BLOCKPATH)
end
% Load timestamps from excel file
if useExcelTS == 1
    excelTS = readmatrix(excelPath, 'Sheet', excelSheetName, 'Range', excelRange);
    onset = excelTS(:,1);
    offset = excelTS(:,2);

    data.epocs.(epocName).onset = onset;
    data.epocs.(epocName).offset = offset;
    data.epocs.(epocName).name = epocName;
    data.epocs.(epocName).data = ones(length(onset));
else
    disp('')
end
ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
if channel == 1
    STREAM_STORE1 = 'x405A';
    STREAM_STORE2 = 'x465A';
elseif channel == 2
    STREAM_STORE1 = 'x405C';
    STREAM_STORE2 = 'x465C';
end

% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event. Use the 'VALUES' parameter to specify allowed values of
% the REF_EPOC to extract.  For stream events, the chunks of data are 
% stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
data = TDTfilter(data, epocName, 'TIME', TRANGE); % extracts data around epoc of interest 

% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), ...
    data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), ...
    data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), ...
    data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), ...
    data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);

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
F405 = blockMeanDownsample(allSignals, N);
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);

% downsample 10x and average 465 signal
allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F465 = blockMeanDownsample(allSignals, N);
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

bls = polyfit(F405(1:end), F465(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;

zall = zeros(size(Y_dF_all));
for i = 1:size(Y_dF_all,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    zall(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
end

dFFall = zeros(size(Y_dF_all));
for i = 1:size(Y_dF_all,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
    dFFall(i,:)=(Y_dF_all(i,:) - zb);
    dFFall(i,:)=(Y_dF_all(i,:)/zb);
end

% Smoothes the z score signal traces
zallSmooth = zeros(size(zall));
for k = 1:height(zall)
    zallSmooth(k,:) = smoothdata(zall(k,:),'movmean',smoothFactor);
end
meanZall_smooth = mean(zallSmooth);

% Smoothes the z score signal traces
dFFallSmooth = zeros(size(dFFall));
for k = 1:height(dFFall)
    dFFallSmooth(k,:) = smoothdata(dFFall(k,:),'movmean',smoothFactor);
end
meandFF_smooth = mean(dFFallSmooth);

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
    disp('Baseline correction applied to smoothed signals')
else
    disp('No baseline correction applied to smoothed signals')
end

% Baseline correction for smoothed data
if setBaseline == 1
    idx = find(ts1>baseAdjust,1);
    for base = 1:height(zall)
        if dFFallSmooth(base,idx) < 0
            val = dFFallSmooth(base,idx);
            diff = 0 - val;
            dFFallSmooth(base,:) = dFFallSmooth(base,:) + abs(diff);
        elseif dFFallSmooth(base,idx) > 0
            val = dFFallSmooth(base,idx);
            diff = 0 - val;
            dFFallSmooth(base,:) = dFFallSmooth(base,:) - abs(diff);
        end
    end
    disp('Baseline correction applied to smoothed signals')
else
    disp('No baseline correction applied to smoothed signals')
end

% Baseline correction
if setBaseline == 1
    idx = find(ts1>baseAdjust,1);
    for base = 1:height(zall)
        if zall(base,idx) < 0
            val = zall(base,idx);
            diff = 0 - val;
            zall(base,:) = zall(base,:) + abs(diff);
        elseif zall(base,idx) > 0
            val = zall(base,idx);
            diff = 0 - val;
            zall(base,:) = zall(base,:) - abs(diff);
        end
    end
    disp('Baseline correction applied')
else
    disp('No baseline correction applied')
end



if plotSmoothed == 0
    % Standard error of the z-score
    meanZall = mean(zall);
    zerror = std(zall)/sqrt(size(zall,1));
    
    % Plot heat map
    subplot(3,1,2);
    imagesc(ts2, 1, zall);
    colormap('jet'); % c1 = colorbar; 
    title(sprintf('Z-Score Heat Map', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14);
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
    title(sprintf('465 nm Z-Score', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    %c2 = colorbar;
    %%
    figure(2)
    plot(ts2, zall)
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('Z-score', 'FontSize', 12)
    title(sprintf('465 nm Z-Score', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    
    %%
    f3 = figure(3);
    subplot(2,3,[1,2,4,5])
    plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    h = fill(XX, YY, 'b');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',18)
    ylabel('Z-score +/- SEM', 'FontSize', 18)
    box off
    
    subplot(2,3,6);
    imagesc(ts2, 1, zall);
    colormap('jet'); colorbar; 
    title(sprintf('Z-Score/Trial', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 16);
    xlabel('Time, s', 'FontSize', 12);
    ylabel('Trial', 'FontSize', 12);
    
    % Fill band values for second subplot. Doing here to scale onset bar
    % correctly
    XX = [ts2, fliplr(ts2)];
    YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
    
    % Saves figure
    if saveFig == 1
        % Save figure 3 to saveFigPath
        file_name1 = char(strcat(figSavePath,TITLE,'.pdf'));
        orient(f3,'landscape');
        print(f3,file_name1,'-dpdf','-vector','-bestfit','');
    else
        disp('Figure not saved')
    end
elseif plotSmoothed == 1
    % Standard error of the z-score
    meanZall = mean(zallSmooth);
    zerror = std(zallSmooth)/sqrt(size(zallSmooth,1));
    
    % Plot heat map
    subplot(3,1,2);
    imagesc(ts2, 1, zall);
    colormap('jet'); % c1 = colorbar; 
    title(sprintf('Z-Score Heat Map', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14);
    ylabel('Trials', 'FontSize', 12);
    
    % Fill band values for second subplot. Doing here to scale onset bar
    % correctly
    XX = [ts2, fliplr(ts2)];
    YY = [mean(zallSmooth)-zerror, fliplr(mean(zallSmooth)+zerror)];
    
    subplot(3,1,3)
    plot(ts2, mean(zallSmooth), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    h = fill(XX, YY, 'r');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('Z-score', 'FontSize', 12)
    title(sprintf('465 nm Z-Score', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    %c2 = colorbar;
    %%
    figure(2)
    plot(ts2, zallSmooth)
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',12)
    ylabel('Z-score', 'FontSize', 12)
    title(sprintf('465 nm Z-Score', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
    
    %%
    f3 = figure(3);
    subplot(2,3,[1,2,4,5])
    plot(ts2, mean(zallSmooth), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
    line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
    
    h = fill(XX, YY, 'b');
    set(h, 'facealpha',.25,'edgecolor','none')
    
    % Finish up the plot
    axis tight
    xlabel('Time, s','FontSize',18)
    ylabel('Z-score +/- SEM', 'FontSize', 18)
    box off
    
    subplot(2,3,6);
    imagesc(ts2, 1, zallSmooth);
    colormap('jet'); colorbar; 
    title(sprintf('Z-Score/Trial', ...
        numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 16);
    xlabel('Time, s', 'FontSize', 12);
    ylabel('Trial', 'FontSize', 12);
    
    % Fill band values for second subplot. Doing here to scale onset bar
    % correctly
    XX = [ts2, fliplr(ts2)];
    YY = [mean(zallSmooth)-zerror, fliplr(mean(zallSmooth)+zerror)];
    
end

meandFF_smooth = meandFF_smooth';
dFFallSmooth = dFFallSmooth';
zallSmooth = zallSmooth';
meanZallSmooth = meanZall_smooth';
ts1 = ts1';
