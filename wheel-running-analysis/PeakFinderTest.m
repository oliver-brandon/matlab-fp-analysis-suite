% MAD vs RMS peak finder comparison &
clear all; clc; close all;

BLOCKPATH = '/Users/brandon/Desktop/DA_WHEEL/DA80_L_Day7_12_14_2022';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
ISOS1 = 'x405A';
DLS = 'x465A';
ISOS2 = 'x405C';
NAC = 'x465C';
% time window in seconds to analyze peaks (30s window works best)%
X = [1025,1055];
% X = [1230,1300];
peakDist = 0.2;
%time array used for all streams%
time1 = (1:length(data.streams.(DLS).data))/data.streams.(DLS).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 10; % time threshold below which we will discard
ind = find(time1>t,1);% find first index of when time crosses threshold
time1 = time1(ind:end); % reformat vector to only include allowed time
data.streams.(DLS).data = data.streams.(DLS).data(ind:end);
data.streams.(ISOS1).data = data.streams.(ISOS1).data(ind:end);
data.streams.(NAC).data = data.streams.(NAC).data(ind:end);
data.streams.(ISOS2).data = data.streams.(ISOS2).data(ind:end);
% convert raw mV to dFF %
bls = polyfit(data.streams.(ISOS1).data,data.streams.(DLS).data,1);
Y_fit_all = bls(1) .* data.streams.(ISOS1).data + bls(2);
Y_dF_all = data.streams.(DLS).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
ISOS_405A = detrend(data.streams.(ISOS1).data);
DLS_465A = detrend(dFF);

bls2 = polyfit(data.streams.(ISOS2).data,data.streams.(NAC).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(ISOS2).data + bls2(2);
Y_dF_all2 = data.streams.(NAC).data - Y_fit_all2; %dF (units mV) is not dFF
dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
std_dFF2 = std(double(dFF2));
ISOS_405C = detrend(data.streams.(ISOS2).data);
NAC_465C = detrend(dFF2);

%downsample streams and time array by N times%
N = 100;
ISOS_405A = downsample(ISOS_405A, N);
DLS_465A = downsample(DLS_465A, N);
ISOS_405C = downsample(ISOS_405C, N);
NAC_465C = downsample(NAC_465C, N);
time1 = downsample(time1, N);
% DLS_465A = zscore(DLS_465A);
% NAC_465C = zscore(NAC_465C);
MAD1 = mad(DLS_465A,1);
MAD2 = mad(NAC_465C,1);
RMS1 = rms(DLS_465A);
RMS2 = rms(NAC_465C);

% Trimming the stream into X second window %
trimstart = find(time1>X(:,1),1);
trimend = find(time1>X(:,2),1);
DLS_trimmed = DLS_465A(trimstart:trimend);
NAC_trimmed = NAC_465C(trimstart:trimend);
time_trimmed = time1(trimstart:trimend);
DLS_trimmed_smooth = smoothdata(DLS_trimmed,'lowess');
NAC_trimmed_smooth = smoothdata(NAC_trimmed,'lowess');
[pks,locs] = findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD1,"MinPeakDistance",peakDist);
[pks2,locs2] = findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD2,"MinPeakDistance",peakDist);
[pks3,locs3] = findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakHeight",RMS1,"MinPeakDistance",peakDist);
[pks4,locs4] = findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakHeight",RMS2,"MinPeakDistance",peakDist);
[pks5,locs5] = findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD1,"MinPeakDistance",peakDist);
[pks6,locs6] = findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD2,"MinPeakDistance",peakDist);
[pks7,locs7] = findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakProminence",RMS1,"MinPeakDistance",peakDist);
[pks8,locs8] = findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakProminence",RMS2,"MinPeakDistance",peakDist);

% MAD w/ min peak height %
fontSize = 18;
figure;
subplot(4,1,1);
findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD1,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Median Absolute Deviation (Height)", "DLS 465A","FontSize",fontSize);
subplot(4,1,2);
findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakHeight",RMS1,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Root Mean Squared (Height)", "DLS 465A","FontSize",fontSize);
% RMS w/ min peak height %
subplot(4,1,3);
findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD2,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Median Absolute Deviation (Height)", "NAc 465C","FontSize",fontSize);
subplot(4,1,4);
findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakHeight",RMS2,"MinPeakDistance",peakDist);
xlabel("Time (s)","FontSize",fontSize);
ylabel("dFF","FontSize",fontSize);
title("Root Mean Squared (Height)", "NAc 465C","FontSize",fontSize);
% MAD w/ min peak prominence %
figure;
subplot(4,1,1);
findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD1,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Median Absolute Deviation (Prominence)", "DLS 465A","FontSize",fontSize);
subplot(4,1,2);
findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakProminence",RMS1,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Root Mean Squared (Prominence)", "DLS 465A","FontSize",fontSize);
% RMS w/ min peak prominence %
subplot(4,1,3);
findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD2,"MinPeakDistance",peakDist);
ylabel("dFF","FontSize",fontSize);
title("Median Absolute Deviation (Prominence)", "NAc 465C","FontSize",fontSize);
subplot(4,1,4);
findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakProminence",RMS2,"MinPeakDistance",peakDist);
xlabel("Time (s)","FontSize",fontSize);
ylabel("dFF","FontSize",fontSize);
title("Root Mean Squared (Prominence)", "NAc 465C","FontSize",fontSize);



% figure;
% subplot(4,1,1);
% findpeaks(DLS_trimmed,time_trimmed,"MinPeakHeight",MAD1,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Height)", "DLS 465A","FontSize",fontSize);
% subplot(4,1,2);
% findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD1,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Height + Gaussian)", "DLS 465A","FontSize",fontSize);
% subplot(4,1,3);
% findpeaks(NAC_trimmed,time_trimmed,"MinPeakHeight",MAD2,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Height)", "NAc 465C","FontSize",fontSize);
% subplot(4,1,4);
% findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakHeight",MAD2,"MinPeakDistance",peakDist);
% xlabel("Time (s)","FontSize",fontSize);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Height + Gaussian)", "NAc 465C","FontSize",fontSize);
% 
% figure;
% subplot(4,1,1);
% findpeaks(DLS_trimmed,time_trimmed,"MinPeakProminence",MAD1,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Prom)", "DLS 465A","FontSize",fontSize);
% subplot(4,1,2);
% findpeaks(DLS_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD1,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Prom + Gaussian)", "DLS 465A","FontSize",fontSize);
% subplot(4,1,3);
% findpeaks(NAC_trimmed,time_trimmed,"MinPeakProminence",MAD2,"MinPeakDistance",peakDist);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Prom)", "NAc 465C","FontSize",fontSize);
% subplot(4,1,4);
% findpeaks(NAC_trimmed_smooth,time_trimmed,"MinPeakProminence",MAD2,"MinPeakDistance",peakDist);
% xlabel("Time (s)","FontSize",fontSize);
% ylabel("dFF","FontSize",fontSize);
% title("Median Absolute Deviation (Prom + Gaussian)", "NAc 465C","FontSize",fontSize);

figure;
subplot(2,1,1);
plot(time_trimmed,DLS_trimmed_smooth,"Color",'blue');
ylabel("dFF","FontSize",fontSize);
title("DLS","FontSize",fontSize);
axis padded
subplot(2,1,2);
plot(time_trimmed,NAC_trimmed_smooth,"Color",'red');
xlabel("Time (s)","FontSize",fontSize);
ylabel("dFF","FontSize",fontSize);
title("NAc","FontSize",fontSize);
axis padded
disp("MinPeakHeight")
MAD_DLS_Pks = length(pks)
RMS_DLS_Pks = length(pks3)
MAD_NAC_Pks = length(pks2)
RMS_NAC_Pks = length(pks4)
disp("MinPeakProminence")
MAD_DLS_Pks2 = length(pks5)
RMS_DLS_Pks2 = length(pks7)
MAD_NAC_Pks2 = length(pks6)
RMS_NAC_Pks2 = length(pks8)

disp("DONE")