clear; clc; close all;
figSavePath = '/Volumes/CUDADRIVE/MielnikCA/20240501_FibreTest/';
saveFig = 1;
setBaseline = 1;
baseAdjust = -2;
BLOCKPATH = '/Volumes/CUDADRIVE/MielnikCA/20240501_FibreTest/1212F_HL_20240501';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
channel = 2;

% data = optoStimEpoch(data,10);
% Creates reward epoc from offset instead of onset
% data.epocs.aReward.onset = data.epocs.aRw_.offset;
% data.epocs.aReward.offset = data.epocs.aRw_.onset;
% data.epocs.aReward.name = 'aReward';
% data.epocs.aReward.data = ones(height(data.epocs.aRw_.offset)) * 10;
% data.epocs.bReward.onset = data.epocs.bRw_.offset;
% data.epocs.bReward.offset = data.epocs.bRw_.onset;
% data.epocs.bReward.name = 'bReward';
% data.epocs.bReward.data = ones(height(data.epocs.bRw_.offset)) * 20;

REF_EPOC = 'St1/'; % Stimulation event to center on

TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
ARANGE = [1 1];
BASELINE_PER = [-2 -1]; % baseline period before stim
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
data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);  
figPath = strcat(BLOCKPATH,'/');
[~,name,~] = fileparts(BLOCKPATH);
brokenID = strsplit(name,'_');
if strcmp(STREAM_STORE1,'x405A')
    ID = brokenID(1);
    task = brokenID(2);
    ROI = 'DLS';
elseif strcmp(STREAM_STORE1,'x405C')
    ID = brokenID(1);
    task = brokenID(2);
    ROI = 'NAc';
else
    disp('Cannot find isosbestic signal. Check the naming and try again.')
end
% remove any "/" from REF_EPOC
REF_EPOC = strrep(REF_EPOC,'/','');
TITLE = strcat(ID,{'-'},task,{'-'},ROI,{'-'},REF_EPOC);
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

zallSmooth = zeros(size(zall));
for k = 1:height(zall)
    zallSmooth(k,:) = smoothdata(zall(k,:),'movmean',50);
end
meanZall_smooth = mean(zallSmooth);

half = floor(height(zall)/2);
if half == 1
    zall_1st = zall(1,:);
    zall_2nd = zall(2,:);
    meanZall_1st = smoothdata(zall_1st,'movmean',50);
    meanZall_2nd = smoothdata(zall_2nd,'movmean',50);
elseif half > 1
    zall_1st = zall(1:half,:);
    zall_2nd = zall(half:end,:);
    meanZall_1st = mean(zall_1st);
    meanZall_2nd = mean(zall_2nd);
    meanZall_1st = smoothdata(meanZall_1st,'movmean',50);
    meanZall_2nd = smoothdata(meanZall_2nd,'movmean',50);
else
    disp('check number of epochs')
end


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
title(TITLE, 'FontSize', 18);
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

if saveFig == 1
    % Save figure 3 to saveFigPath
    file_name1 = char(strcat(figSavePath,TITLE,'.pdf'));
    orient(f3,'landscape');
    print(f3,file_name1,'-dpdf','-vector','-bestfit','');
else
    disp('Figure not saved')
end
