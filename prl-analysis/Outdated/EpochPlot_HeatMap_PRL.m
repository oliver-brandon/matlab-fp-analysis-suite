%% Fiber Photometry Epoch Averaging Example
%
% <html>
% This example goes through fiber photometry analysis using techniques <br>
% such as data smoothing, bleach detrending, and z-score analysis. <br>
% The epoch averaging was done using TDTfilter. <br><br>
% Author Contributions: <br>
% TDT, David Root, and the Morales Lab contributed to the writing and/or conceptualization of the code. <br>
% The signal processing pipeline was inspired by the workflow developed by <a href="https://doi.org/10.1016/j.celrep.2017.10.066">David Barker et al. (2017)</a> for the Morales Lab. <br>
% The data used in the example were provided by David Root. <br><br>
% Author Information: <br>
% David H. Root <br>
% Assistant Professor <br>
% Department of Psychology & Neuroscience <br>
% University of Colorado, Boulder <br>
% Lab Website: <a href="https://www.root-lab.org">https://www.root-lab.org</a> <br>
% david.root@colorado.edu <br><br>
% About the authors: <br>
% The Root lab and Morales lab investigate the neurobiology of reward, aversion, addiction, and depression. <br>
% <br> TDT edits all user submissions in coordination with the contributing
% author(s) prior to publishing.
% </html>

%% Housekeeping
% Clear workspace and close existing figures.
% Add data
clear all; clc; close all;


% [file,path] = uigetfile;
% [~,filename,~] = fileparts(file);
% matfile = fullfile(path,file);
% 
% load(matfile);


BLOCKPATH = '/Volumes/CUDADRIVE/unprocessedTanks/Dual Fiber JZL/465F_Acq2-3_JZL18_Empty_NA_NA';
[~,name,~] = fileparts(BLOCKPATH);
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
box_number = 3; %3 = 465A, 4 = 465C
epoch_number = 3;%1 = cRew, 2 = iNoRew, 3 = cue
training_logic = 2; %1 = training, 2 = Acq/Rev

TITLE1 = name;
if box_number == 3
    STREAM_STORE1 = 'x405A'; % name of the 405 store
    STREAM_STORE2 = 'x465A'; % name of the 465 store
    if epoch_number == 1
        REF_EPOC = 'cRewA'; % Stimulation event to center on
    elseif epoch_number == 2
        REF_EPOC = 'iNoRewA'; % Stimulation event to center on
    elseif epoch_number == 3
        REF_EPOC = 'St1/';
    end
elseif box_number == 4
    STREAM_STORE1 = 'x405C'; % name of the 405 store
    STREAM_STORE2 = 'x465C'; % name of the 465 store
    if epoch_number == 1
        REF_EPOC = 'cRewC'; % Stimulation event to center on
    elseif epoch_number == 2
        REF_EPOC = 'iNoRewC'; % Stimulation event to center on
    elseif epoch_number == 3
        REF_EPOC = 'St2/';
    end
end  
TRANGE = [-5 10]; %window size [start time relative to epoc onset, entire duration]
ARANGE = [1 1];
BASELINE_PER = [-3 -1]; % baseline period before stim
ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
if box_number == 3 && training_logic == 2
    cueTSA = data.epocs.St1_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL1_') == 0
        correct_ts_A = 0;
    else
        correct_ts_A = data.epocs.CL1_.onset;
    end
    
    if isfield(data.epocs, 'IL1_') == 0
        incorrect_ts_A = 0;
    else
        incorrect_ts_A = data.epocs.IL1_.onset;
    end
    pellet_ts_A = data.epocs.Pe1_.onset;
    correct_rewardedA = zeros(double(height(pellet_ts_A)));
    correct_norewardA = zeros(double(height(pellet_ts_A)));
    incorrect_rewardedA = zeros(double(height(pellet_ts_A)));
    incorrect_norewardA = zeros(double(height(pellet_ts_A)));
    for i = 1:height(correct_ts_A)
        var1 = correct_ts_A(i,:);
        var2 = pellet_ts_A;
        [f, ia] = ismember(var1, pellet_ts_A);
        if f == 1
            correct_rewardedA(i,:) = var1;
        elseif f < 1
            correct_norewardA(i,:) = var1;
        end  
    end
    for i = 1:height(incorrect_ts_A)
        var1 = incorrect_ts_A(i,:);
        var2 = pellet_ts_A;
        [f, ia] = ismember(var1, pellet_ts_A);
        if f == 1
            incorrect_rewardedA(i,:) = var1;
        elseif f < 1
            incorrect_norewardA(i,:) = var1;
        end
    end
    
    correct_rewardedA = nonzeros(correct_rewardedA(:,1));
    correct_norewardA = nonzeros(correct_norewardA(:,1));
    incorrect_rewardedA = nonzeros(incorrect_rewardedA(:,1));
    incorrect_norewardA = nonzeros(incorrect_norewardA(:,1));
    
    if isempty(correct_rewardedA) == 1
        correct_rewardedA = 0;
    end
    if isempty(correct_norewardA) == 1
        correct_norewardA = 0;
    end
    if isempty(incorrect_rewardedA) == 1
        incorrect_rewardedA = 0;
    end
    if isempty(incorrect_norewardA) == 1
        incorrect_norewardA = 0;
    end
    data.epocs.cRewA.name = 'cRewA';
    data.epocs.cRewA.onset = correct_rewardedA;
    data.epocs.cRewA.offset = correct_rewardedA+1;
    data.epocs.cRewA.data = ones(height(correct_rewardedA));
    data.epocs.cNoRewA.name = 'cNoRewA';
    data.epocs.cNoRewA.onset = correct_norewardA;
    data.epocs.cNoRewA.offset = correct_norewardA+1;
    data.epocs.cNoRewA.data = ones(height(correct_norewardA))*2;
    data.epocs.iRewA.name = 'iRewA';
    data.epocs.iRewA.onset = incorrect_rewardedA;
    data.epocs.iRewA.offset = incorrect_rewardedA+1;
    data.epocs.iRewA.data = ones(height(incorrect_rewardedA))*3;
    data.epocs.iNoRewA.name = 'iNoRewA';
    data.epocs.iNoRewA.onset = incorrect_norewardA;
    data.epocs.iNoRewA.offset = incorrect_norewardA+1;
    data.epocs.iNoRewA.data = ones(height(incorrect_norewardA))*4;
elseif box_number == 4 && training_logic == 2
    cueTSC = data.epocs.St2_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL2_') == 0
        correct_ts_C = 0;
    else
        correct_ts_C = data.epocs.CL2_.onset;
    end
    
    if isfield(data.epocs, 'IL2_') == 0
        incorrect_ts_C = 0;
    else
        incorrect_ts_C = data.epocs.IL2_.onset;
    end
    pellet_ts_C = data.epocs.Pe2_.onset;
    correct_rewardedC = zeros(double(height(pellet_ts_C)));
    correct_norewardC = zeros(double(height(pellet_ts_C)));
    incorrect_rewardedC = zeros(double(height(pellet_ts_C)));
    incorrect_norewardC = zeros(double(height(pellet_ts_C)));
    for i = 1:height(correct_ts_C)
        var3 = correct_ts_C(i,:);
        var4 = pellet_ts_C;
        [g, ia] = ismember(var3, pellet_ts_C);
        if g == 1
            correct_rewardedC(i,:) = var3;
        elseif g < 1
            correct_norewardC(i,:) = var3;
        end  
    end
    for i = 1:height(incorrect_ts_C)
        var3 = incorrect_ts_C(i,:);
        var4 = pellet_ts_C;
        [g, ia] = ismember(var3, pellet_ts_C);
        if g == 1
            incorrect_rewardedC(i,:) = var3;
        elseif g < 1
            incorrect_norewardC(i,:) = var3;
        end
    end
    correct_rewardedC = nonzeros(correct_rewardedC(:,1));
    correct_norewardC = nonzeros(correct_norewardC(:,1));
    incorrect_rewardedC = nonzeros(incorrect_rewardedC(:,1));
    incorrect_norewardC = nonzeros(incorrect_norewardC(:,1));
    if isempty(correct_rewardedC) == 1
        correct_rewardedC = 0;
    end
    if isempty(correct_norewardC) == 1
        correct_norewardC = 0;
    end
    if isempty(incorrect_rewardedC) == 1
        incorrect_rewardedC = 0;
    end
    if isempty(incorrect_norewardC) == 1
        incorrect_norewardC = 0;
    end
    data.epocs.cRewC.name = 'cRewC';
    data.epocs.cRewC.onset = correct_rewardedC;
    data.epocs.cRewC.offset = correct_rewardedC+1;
    data.epocs.cRewC.data = ones(height(correct_rewardedC));
    data.epocs.cNoRewC.name = 'cNoRewC';
    data.epocs.cNoRewC.onset = correct_norewardC;
    data.epocs.cNoRewC.offset = correct_norewardC+1;
    data.epocs.cNoRewC.data = ones(height(correct_norewardC))*2;
    data.epocs.iRewC.name = 'iRewC';
    data.epocs.iRewC.onset = incorrect_rewardedC;
    data.epocs.iRewC.offset = incorrect_rewardedC+1;
    data.epocs.iRewC.data = ones(height(incorrect_rewardedC))*3;
    data.epocs.iNoRewC.name = 'iNoRewC';
    data.epocs.iNoRewC.onset = incorrect_norewardC;
    data.epocs.iNoRewC.offset = incorrect_norewardC+1;
    data.epocs.iNoRewC.data = ones(height(incorrect_norewardC))*4;
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

