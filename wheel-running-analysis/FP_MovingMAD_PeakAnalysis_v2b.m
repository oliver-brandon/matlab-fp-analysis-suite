clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = "2b";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_MovingMAD_PeakAnalysis.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% input paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 20; % time threshold below which we will discard
N = 10; %Downsample N times
sampleRate = 1017/N;
session_duration = 3600;
window_size_seconds = 15;
madMultiplier = 2;
minPkWdth = 0.2;
%Snippet Args%
snippet_duration = 30; % Duration of the snippet in seconds
snippet_start_time = 1200; % Start time of the snippet in seconds
%Stream Stores%
DLS_ISOS = 'x405A'; % name of the 405A store
DLS_DA = 'x465A'; % name of the 465A store
NAc_ISOS = 'x405C'; % name of the 405C store
NAc_DA = 'x465C'; % name of the 465C store

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Version: %s\n",VERSION)
myDir = uigetdir('H:\My Drive\wheel-peak-test\tanks');
BLOCKPATH = myDir;
data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});

%time array used for all streams%
time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%

ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);

%downsample streams and time array by N times%
data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
time = downsample(time, N);

%detrend & dFF%
%465A%
bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
detrend_465A = detrend(dFF);
detrend_465A = zscore(detrend_465A);
%465C%
bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
std_dFF2 = std(double(dFF2));
detrend_465C = detrend(dFF2);
detrend_465C = zscore(detrend_465C);

%set up moving MAD window%

% determines if window_size evenly divides into the streams and pads with 
% NaN if not
window_size = ceil(window_size_seconds * sampleRate);
remainderA = ceil((window_size - mod(length(detrend_465A),window_size)));
remainderC = ceil((window_size - mod(length(detrend_465C),window_size)));
if remainderA == window_size
    signalA = detrend_465A;
    signalA_ind = 1:length(signalA);
else
    padding = NaN(1,remainderA);
    signalA = [detrend_465A,padding];
    signalA_ind = 1:length(signalA);
end

if remainderC == window_size
    signalC = detrend_465C;
    signalC_ind = 1:length(signalC);
else
    padding = NaN(1,remainderC);
    signalC = [detrend_465C,padding];
    signalC_ind = 1:length(signalC);
end
time_pad = (1:length(signalA))/sampleRate;

signalA_pks = [];
for i = 1:window_size:length(signalA)
    if (i + window_size - 1) > length(signalA)
        window_end = length(signalA);
    else
        window_end = i + window_size - 1;
    end
    
    sigA_med = median(signalA(i:window_end));
    sigA_mad = mad(signalA(i:window_end), 1);
    sigA_thr = (sigA_med + (madMultiplier * sigA_mad));

    [pksA, locsA] = findpeaks(signalA(i:window_end), ...
        time_pad(i:window_end), 'MinPeakHeight', sigA_thr, ...
        'MinPeakDistance', minPkWdth);

    signalA_pks = [signalA_pks, locsA];
end

signalC_pks = [];
for i = 1:window_size:length(signalC)
    if (i + window_size - 1) > length(signalC)
        window_end = length(signalC);
    else
        window_end = i + window_size - 1;
    end
    
    sigC_med = median(signalC(i:window_end));
    sigC_mad = mad(signalC(i:window_end), 1);
    sigC_thr = (sigC_med + (madMultiplier * sigC_mad));
    % sigC_thr = madMultiplier * sigC_mad;

    [pksC, locsC] = findpeaks(signalC(i:window_end), ...
        time_pad(i:window_end), 'MinPeakHeight', sigC_thr, ...
        'MinPeakDistance', minPkWdth);

    signalC_pks = [signalC_pks, locsC];
end

x1 = ceil(time(1,1));
x2 = ceil(time(1,end));
% Plotting
f1 = figure;
plot(time_pad, signalA)
xlim([x1 x2]);
hold on
peak_indices = ismember(time_pad, signalA_pks);
plot(time_pad(peak_indices), signalA(peak_indices), 'ro')
hold off

f2 = figure;
plot(time_pad, signalC)
xlim([x1 x2]);
hold on
peak_indices = ismember(time_pad, signalC_pks);
plot(time_pad(peak_indices), signalC(peak_indices), 'ro')
hold off


%% Snippet Plotting %%
% Find the indices corresponding to the start and end times of the snippet
[~, snippet_start_index] = min(abs(time_pad - snippet_start_time));
snippet_end_index = snippet_start_index + round(snippet_duration / (time_pad(2) - time_pad(1))) - 1;

% Ensure the end index does not exceed the length of the signal
snippet_end_index = min(snippet_end_index, length(signalA));

% SIGNAL A %
% Plot the snippet of the signal
f3 = figure;
plot(time_pad(snippet_start_index:snippet_end_index), signalA(snippet_start_index:snippet_end_index))
hold on

% Plot the peaks within the snippet
peak_indices_within_snippet = peak_indices & (time_pad >= snippet_start_time) & (time_pad <= (snippet_start_time + snippet_duration));
plot(time_pad(peak_indices_within_snippet), signalA(peak_indices_within_snippet), 'ro')
hold off

% SIGNAL C %
% Plot the snippet of the signal
f4 = figure;
plot(time_pad(snippet_start_index:snippet_end_index), signalC(snippet_start_index:snippet_end_index))
hold on

% Plot the peaks within the snippet
peak_indices_within_snippet = peak_indices & (time_pad >= snippet_start_time) & (time_pad <= (snippet_start_time + snippet_duration));
plot(time_pad(peak_indices_within_snippet), signalC(peak_indices_within_snippet), 'ro')
hold off






















% time2 = (1:length(signalA))/sampleRate;
% snipSt = find(time2>snip(1),1);
% snipEn = find(time2>snip(2),1);
% 
% % fills signal_chunksA and signal_chunksC with the actual signal chunks
% for ii = 1:window_size:length(signalA)
%     x1 = signalA(1,ii:ii+window_size-1);
%     signal_chunksA = [signal_chunksA; x1];
% end
% for ii = 1:window_size:length(signalC)
%     x1 = signalC(1,ii:ii+window_size-1);
%     signal_chunksC = [signal_chunksC; x1];
% end

% % finds peaks using previous window
% signalA_pks = [];
% % signalA_pks(1,size(signal_chunksA,2)) = nan;
% for jj = 2:height(signal_chunksA)
%     signalA_median = median(signal_chunksA(jj-1,:));
%     signalA_MAD = mad(signal_chunksA(jj-1,:),1);
%     signalA_thresh = (signalA_median + (3 * signalA_MAD));
%     % pks = signal_chunksA(jj,:) > signalA_thresh;
%     % signalA_pks(jj,:) = pks;
% 
% end
% 
% % finds peaks using previous window
% signalC_pks = [];
% signalC_pks(1,size(signal_chunksC,2)) = nan;
% for jj = 2:height(signal_chunksC)
%     signalC_median = median(signal_chunksC(jj-1,:));
%     signalC_MAD = mad(signal_chunksC(jj-1,:),1);
%     signalC_thresh = (signalA_median + (3 * signalC_MAD));
%     pks2 = signal_chunksC(jj,:) > signalC_thresh;
% 
%     signalC_pks(jj,:) = pks2;
% end

% % reshapes signal peak matrix into size of signal
% signalA_pks = [signalA_pks(1, :) reshape(signalA_pks(2:end, :)', 1, [])];
% signalC_pks = [signalC_pks(1, :) reshape(signalC_pks(2:end, :)', 1, [])];

        





% 
% % calculates number of peaks and frequency
% totPks_A = length(signalA_pks(signalA_pks == 1));
% transients_A = totPks_A/session_duration;
% totPks_C = length(signalC_pks(signalC_pks == 1));
% transients_C = totPks_C/session_duration;

% % Plots the signal with markers for the peaks
% f1 = figure;
% plot(time2, signalA);
% hold on;
% 
% plot(time2(signalA_pks == 1), signalA(signalA_pks == 1), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
% 
% hold off;
% 
% time_trim = time2(1, snipSt:snipEn);
% signalA_trim = signalA(1,snipSt:snipEn);
% signalA_pks_trim = signalA_pks(1, snipSt:snipEn);
% 
% f2 = figure;
% plot(time_trim, signalA_trim);
% hold on;
% plot(time_trim(signalA_pks_trim == 1), signalA_trim(signalA_pks_trim == 1), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
% hold off;


