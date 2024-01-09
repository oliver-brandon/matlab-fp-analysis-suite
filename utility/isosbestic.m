clear all; close all;
TANKPATH = '/Users/brandon/My Drive (bloliv95@gmail.com)/self_admin/coc_sa/PrL-aIC/signal_test/tanks/883_PrL_884_PrL';
data = TDTbin2mat(TANKPATH);

REF_EPOC = 'St1/';
STREAM_STORE1 = 'x405A'; % name of the 405 store
STREAM_STORE2 = 'x465A'; % name of the 465 store
TITLE = '883 - PrL GCaMP6f - NAc rCre (Houselight)';
TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-3 -1]; % baseline period before stim
ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
N = 10; % Downsample Nx

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

allSignals405 = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample and average 405 signal
F405 = zeros(size(allSignals405(:,1:N:end-N+1)));
for ii = 1:size(allSignals405,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals405(ii,i:i+N-1)),1:N:length(allSignals405)-N+1);
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

% Standard error of the z-score
meanZall = mean(zall);
zerror = std(zall)/sqrt(size(zall,1));


mean405 = mean(F405);
mean405 = deltaFF(mean405);
mean405 = zScore(mean405);
mean465 = mean(F465);
mean465 = deltaFF(mean465);
mean465 = zScore(mean465);

error405 = std(mean405)/sqrt(size(mean405,1));
error465 = std(mean465)/sqrt(size(mean465,1));

XX = [ts1, fliplr(ts1)];
YY405 = [mean405 + error405, fliplr(mean405 - error405)];
YY465 = [mean465 + error465, fliplr(mean465 - error465)];


f1 = figure;
plot(ts1,mean405,'Color','r','DisplayName', 'Mean Isosbestic');
hold on
plot(ts1,mean465,'Color','b', 'DisplayName', 'Mean GCaMP6f');

h = fill(XX, YY405, 'r', 'DisplayName', 'SEM Isosbestic');
g = fill(XX, YY465, 'b', 'DisplayName', 'SEM GCaMP6f');
set(h, 'facealpha',.20,'edgecolor','none')
set(g, 'facealpha',.20,'edgecolor','none')



xline(0,'LineWidth',2,'Color','black','DisplayName','Epoc Onset')
xlabel('Time (s)');
ylabel('Normalized GCaMP6f');
title(TITLE);
legend('Location','best');

hold off



full405 = data.streams.(STREAM_STORE1).data;
full465 = data.streams.(STREAM_STORE2).data;


full405 = downsample(full405,N);
full405 = detrend(full405);

full465 = downsample(full465,N);
full465 = detrend(full465);

full405 = deltaFF(full405);
full405 = zScore(full405);
full465 = deltaFF(full465);
full465 = zScore(full465);
time = (1:length(full405))/(1017/N);
t = 30;
ind = find(time>t,1);
full405 = full405(ind:end);
full465 = full465(ind:end);
time = time(ind:end);

% [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(data.epocs.aRL_.offset, 20);
ind1 = find(time>300,1);
ind2 = find(time>360,1);

full405 = full405(ind1:ind2);
full465 = full465(ind1:ind2);
time = time(ind1:ind2);
f2 = figure;

plot(time,full405,'Color','r','DisplayName','Isosbestic')
hold on
plot(time,full465,'Color','b','DisplayName','GCaMP6f')

% for kk = 1:height(rewardTimestamps)
%     if rewardTimestamps(kk) < 1000 || rewardTimestamps(kk) > 1200
%         continue
%         xline(rewardTimestamps(kk),'LineWidth',0.5,'Color','black','HandleVisibility','off')
%     end
% end

xlabel('Time (s)')
ylabel('Normalized Signal')
legend('Location','best')
hold off