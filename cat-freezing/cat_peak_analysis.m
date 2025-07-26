% Updated 4/30/2025
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
N = 100; % downsample
sigHz = 1017/N;

timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % the number of seconds before the onset of a TTL to analyze
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
channel = 1; % 1 = A, 2 = C
t = 60; % remove t seconds from start of recording
smoothFactor = 5;
BLOCKPATH = '/Users/brandon/ucr-drive/collabs/Yaminaka/20250213_FPSF1-05_FPSF1-06';

useSmallMADWindow = 1; % 0 = no, 1 = yes
MADwindow = [t 1000]; % MAD window in seconds [start end]

useExcelTS = 1; % 1 = yes, 0 = no
excelPath = '/Users/brandon/ucr-drive/collabs/Yaminaka/FPSF1-05 Behavior Timestamps.xlsx';
excelSheetName = '-05 Control';
excelRange = 'G2:H9';
epocName = 'ControlFreeze';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Starting analysis...")

data = TDTbin2mat(BLOCKPATH);
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

if channel == 1
    ISOS = 'x405A';
    SIGNAL = 'x465A';
elseif channel == 2
    ISOS = 'x405C';
    SIGNAL = 'x465C';
end

time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
ind = find(time>t,1);
time = time(ind:end);
data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
min_time = min(time);

%downsample
data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data,N);
data.streams.(ISOS).data = downsample(data.streams.(ISOS).data,N);
time = downsample(time,N);

%detrend & dFF%
bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,1);
Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
dFF_signal = detrend(dFF);
dFF_signal_smooth = smoothdata(dFF_signal,'movmean',smoothFactor);

%calculates and plots median absolute deviation for both 465 signals%
if useSmallMADWindow == 0
    MAD1 = mad(dFF_signal_smooth, 1);
elseif useSmallMADWindow == 1
    indSt = find(time>MADwindow(:,1),1);
    indEn = find(time>MADwindow(:,2),1);
    MAD1 = mad(dFF_signal_smooth, 1);
else
    error('Unknown MAD window')
end
%[pks,locs,w,p] = findpeaks(dFF_signal_smooth,time,'MinPeakHeight',MAD1);
findpeaks(dFF_signal_smooth,time,'MinPeakProminence',MAD1);
% findpeaks(dFF_signal_smooth,time,'MinPeakHeight',MAD1);
epoc_ts = [data.epocs.(epocName).onset data.epocs.(epocName).offset];
%peak analysis
for j = 1:height(epoc_ts)
    epocStart = epoc_ts(j,1);
    epocEnd = epoc_ts(j,2);
    indSt = find(time>epocStart,1);
    indEn = find(time>epocEnd,1);
    sig_epoc = dFF_signal_smooth(1,indSt:indEn);
    time_epoc = time(1,indSt:indEn);

    [pks2,locs2,w2,p2] = findpeaks(sig_epoc,time_epoc,'MinPeakHeight',MAD1);
    if isempty(pks2)
        pks2 = 0;
    end
    sig_epoc_numpeak(j,:) = length(pks2);
    sig_epoc_peakfreq(j,:) = (sig_epoc_numpeak(j,:)/(epocEnd-epocStart));
    sig_epoc_max_amp(j,:) = max(pks2);
    sig_epoc_avg_amp(j,:) = mean(pks2);


end
peak_analysis_table = table(epoc_ts(:,1), sig_epoc_numpeak, sig_epoc_peakfreq, ...
    sig_epoc_max_amp, sig_epoc_avg_amp, 'VariableNames',{'TS','NumPeaks',...
    'Peaks/s','Peak Max','Peak Avg'});
ts1 = -baseWindow + (1:epocArrayLen)/data.streams.(SIGNAL).fs*N;
for k = 1:height(epoc_ts)
    epocStart = epoc_ts(k,1) - baseWindow;
    epocEnd = epocStart + timeWindow + baseWindow;
    indSt = find(time>epocStart,1);
    indEn = find(time>epocEnd,1);
    sig_epoc = dFF_signal_smooth(1,indSt:indEn);
    
    if length(sig_epoc) < epocArrayLen
        mn = mean(sig_epoc(1,end-10:end));
        sig_epoc(1,end:epocArrayLen) = mn;
    elseif length(sig_epoc) > epocArrayLen
        op = length(sig_epoc);
        arrayDif = op - epocArrayLen;
        sig_epoc = sig_epoc(1,1:end-arrayDif);
    end
    sig_epoc_all(k,:) = sig_epoc;
    if sig_epoc_all(k,1) < 0
        val = sig_epoc_all(k,1);
        diff = 0-val;
        sig_epoc_all(k,1:epocArrayLen) = sig_epoc_all(k,1:epocArrayLen) + abs(diff);
    elseif sig_epoc_all(k,1) > 0
        val = sig_epoc_all(k,1);
        diff = 0 - val;
        sig_epoc_all(k,1:epocArrayLen) = sig_epoc_all(k, 1:epocArrayLen) - abs(diff);
    end

end
mean_sig_all = mean(sig_epoc_all,1)';
sig_epoc_all = sig_epoc_all';
ts1 = ts1';