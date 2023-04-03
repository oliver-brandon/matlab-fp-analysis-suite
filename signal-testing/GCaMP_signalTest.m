% GCaMP_signalTest.m created by Brandon L. Oliver, M.A.
% Used to quickly test GCaMP flourescence from TDT fiberphotometry tanks.
% Plots the raw, dFF, and detrended dFF on seperate plots and saves a
% figure to a user designated directory. Can also downsample if necessary
% by changing the value of 'N' below.
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figsavepath = 'Z:\'; % must include backslash at the end of the path
figsavetype = '.pdf'; % can change to '.jpg', '.fig', etc.
t = 5; % first t seconds are discarded to remove laser on artifact
N = 1; % downsample signal N times
ISOS = 'x405C'; % set name of isosbestic signal
GCaMP = 'x465C'; % set name of GCaMP signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gets tank from UI pop-up window
TANK_NAME = uigetdir(pwd, 'Select a tank to plot');
[~,name,~] = fileparts(TANK_NAME);
TITLE = strrep(name,'_',' ');
data = TDTbin2mat(TANK_NAME, 'TYPE', {'streams'});
ISOS_raw = data.streams.(ISOS).data;
GCaMP_raw = data.streams.(GCaMP).data;

% Checks for unequal isosbestic and GCaMP signal length correcting if
% necessary
if length(GCaMP_raw) < length(ISOS_raw)
    disp('Isosbestic signal array is longer than GCaMP signal array')
    ISOS_raw = ISOS_raw(1:length(GCaMP_raw));
    disp('Corrected.')
elseif length(GCaMP_raw) > length(ISOS_raw)
    disp('Isosbestic signal array is shorter than GCaMP signal array')
    GCaMP_raw = GCaMP_raw(1:length(ISOS_raw));
    disp('Corrected.')
else
    disp('Isosbestic and GCaMP signal arrays are of equal size')
    disp('No correction necessary.')
end

% time array
time = (1:length(GCaMP_raw))/data.streams.(GCaMP).fs;

% removes the first (t) seconds of signal
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
GCaMP_raw = GCaMP_raw(ind:end);
ISOS_raw = ISOS_raw(ind:end);

% Downsample streams and time array by N times%
ISOS_raw = downsample(ISOS_raw, N);
GCaMP_raw = downsample(GCaMP_raw, N);
time = downsample(time, N);

% Converts raw mV signal to dFF and detrends to remove photobleaching%
bls = polyfit(ISOS_raw,GCaMP_raw,1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = GCaMP_raw - Y_fit_all; %dF (units mV) is not dFF
GCaMP_dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(GCaMP_dFF));
GCaMP_dFF_detrend = detrend(GCaMP_dFF);

ind1 = find(time > 55,1);
ind2 = find(time > 65,1);
GCaMP_snip = GCaMP_dFF_detrend(ind1:ind2);
time_snip = time(ind1:ind2);
% Creates figure of raw, dFF, and detrended dFF GCaMP signal
f1 = figure;
subplot(4,1,1)
plot(time,GCaMP_raw,'g')
title(TITLE)
ylabel("Raw GCaMP6 (mV)")
subplot(4,1,2)
plot(time,GCaMP_dFF,'g')
ylabel("GCaMP6 (dFF)")
subplot(4,1,3)
plot(time,GCaMP_dFF_detrend,'g')
ylabel('GCaMP6 Detrended (dFF)')
subplot(4,1,4)
plot(time_snip,GCaMP_snip,'g')
title('10s Window')
ylabel('GCaMP6 dFF')
xlabel("Time (s)")

% Saves figure as desired file type to user chosen directory
file_name = strcat(figsavepath,name,figsavetype);
saveas(f1,file_name)