clear;close all;
t = 5; % first t seconds are discarded to remove LED on artifact
N = 100; % downsample signal N times
channel = 1; % 1 = A, 2 = C

BLOCKPATH = '/Users/brandon/personal-drive/prl/GrabDA-eCB/mats/2145_Rev1_NA.mat';
load(BLOCKPATH);
%data = TDTbin2mat(BLOCKPATH);
if channel == 1
    ISOS = 'x405A'; % set name of isosbestic signal
    Grab = 'x560A'; % set name of Grab signal
elseif channel == 2
    ISOS = 'x405C'; % set name of isosbestic signal
    Grab = 'x465C'; % set name of Grab signal
else
    disp('Cannot find isosbestic signal. Check the naming and try again.')
end

ISOS_raw = data.streams.(ISOS).data;
Grab_raw = data.streams.(Grab).data;

% Checks for unequal isosbestic and Grab sensor signal length correcting if
% necessary
if length(Grab_raw) < length(ISOS_raw)
    disp('Isosbestic signal array is longer than Grab signal array')
    ISOS_raw = ISOS_raw(1:length(Grab_raw));
    disp('Corrected.')
elseif length(Grab_raw) > length(ISOS_raw)
    disp('Isosbestic signal array is shorter than Grab signal array')
    Grab_raw = Grab_raw(1:length(ISOS_raw));
    disp('Corrected.')
else
    disp('Isosbestic and Grab signal arrays are of equal size')
    disp('No correction necessary.')
end

% time array
time = (1:length(Grab_raw))/data.streams.(Grab).fs;

% removes the first (t) seconds of signal
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
Grab_raw = Grab_raw(ind:end);
ISOS_raw = ISOS_raw(ind:end);

% Downsample streams and time array by N times%
ISOS_raw = downsample(ISOS_raw, N);
Grab_raw = downsample(Grab_raw, N);
time = downsample(time, N);

% Converts raw mV signal to dFF and detrends to remove photobleaching%
bls = polyfit(ISOS_raw,Grab_raw,1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = Grab_raw - Y_fit_all; %sF (units mV) is not dFF
Grab_dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(Grab_dFF));
Grab_dFF = detrend(Grab_dFF);
Grab_dFF_z = zscore(Grab_dFF);

% Plot the z-scored dFF signal
figure;
plot(time, Grab_dFF_z);
xlabel('Time (s)');
ylabel('Z-scored dFF');
grid on;

% Plot the raw uncorrected signal
figure;
plot(time,Grab_raw)
xlabel('Time (s)')
ylabel('Uncorrected Signal (mV)')
grid on;