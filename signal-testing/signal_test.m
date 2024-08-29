clear all;
close all;
% signal_test.m created by Brandon L. Oliver, M.A.
% Used to quickly test flourescent signals acquired from TDT fiberphotometry tanks.
% Plots the raw, dFF, and detrended dFF on seperate plots and saves a
% figure to a user designated directory. Can also downsample if necessary
% by changing the value of 'N' below.
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figsavepath = '/Volumes/CUDADRIVE/signal-test-tanks/grab-da/20240808-NAc-GrabDA/'; % must include backslash at the end of the path
figsavetype = '.pdf'; % can change to '.jpg', '.fig', etc.
t = 5; % first t seconds are discarded to remove laser on artifact
N = 10; % downsample signal N times
channel = 1; % 1 = mouse on A channel, 2 = mouse on C channel
figSnip = [100, 160];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gets tank from UI pop-up window
TANK_NAME = uigetdir(pwd, 'Select a tank to plot');
[~,name,~] = fileparts(TANK_NAME);
brokenID = strsplit(name,'_');
if channel == 1
    ISOS = 'x405A'; % set name of isosbestic signal
    SIGNAL = 'x465A'; % set name of SIGNAL signal
    animalID = char(brokenID{1});
    region = char(brokenID{2});
elseif channel == 2
    ISOS = 'x405C'; % set name of isosbestic signal
    SIGNAL = 'x465C'; % set name of SIGNAL signal
    animalID = char(brokenID{3});
    region = char(brokenID{4});
end

TITLE = strcat(animalID," ",region);
data = TDTbin2mat(TANK_NAME, 'TYPE', {'streams'});
ISOS_raw = data.streams.(ISOS).data;
SIGNAL_raw = data.streams.(SIGNAL).data;

% Checks for unequal isosbestic and SIGNAL signal length correcting if
% necessary
if length(SIGNAL_raw) < length(ISOS_raw)
    disp('Isosbestic signal array is longer than SIGNAL signal array')
    ISOS_raw = ISOS_raw(1:length(SIGNAL_raw));
    disp('Corrected.')
elseif length(SIGNAL_raw) > length(ISOS_raw)
    disp('Isosbestic signal array is shorter than SIGNAL signal array')
    SIGNAL_raw = SIGNAL_raw(1:length(ISOS_raw));
    disp('Corrected.')
else
    disp('Isosbestic and SIGNAL signal arrays are of equal size')
    disp('No correction necessary.')
end

% time array
time = (1:length(SIGNAL_raw))/data.streams.(SIGNAL).fs;

% removes the first (t) seconds of signal
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
SIGNAL_raw = SIGNAL_raw(ind:end);
ISOS_raw = ISOS_raw(ind:end);

% Downsample streams and time array by N times%
ISOS_raw = downsample(ISOS_raw, N);
SIGNAL_raw = downsample(SIGNAL_raw, N);
time = downsample(time, N);

% Converts raw mV signal to dFF and detrends to remove photobleaching%
bls = polyfit(ISOS_raw,SIGNAL_raw,1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = SIGNAL_raw - Y_fit_all; %dF (units mV) is not dFF
SIGNAL_dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(SIGNAL_dFF));
SIGNAL_dFF_detrend = detrend(SIGNAL_dFF);
SIGNAL_Z = zScore(SIGNAL_dFF_detrend);

ind1 = find(time > figSnip(1),1);
ind2 = find(time > figSnip(2),1);
SIGNAL_snip = SIGNAL_Z(ind1:ind2);
time_snip = time(ind1:ind2);

x1 = ceil(time(1,1));
x2 = ceil(time(1,end));
% Creates figure of raw, dFF, and detrended dFF SIGNAL signal
f1 = figure;
subplot(4,1,1)
plot(time,ISOS_raw,'r')
xlim([x1 x2]);
title(TITLE)
ylabel("Raw Isosbestic (mV)")
subplot(4,1,2)
plot(time,SIGNAL_raw,'b')
xlim([x1 x2]);
ylabel("Signal (mV)")
subplot(4,1,3)
plot(time,SIGNAL_Z,'b')
xlim([x1 x2]);
ylabel('Normalized Signal')
subplot(4,1,4)
plot(time_snip,SIGNAL_snip,'b')
title(sprintf('%d Second Window', (figSnip(2)-figSnip(1))));
ylabel('Normalized Signal')
xlabel("Time (s)")
xlim([figSnip(1) figSnip(2)]);

% Saves figure as desired file type to user chosen directory
file_name = strcat(figsavepath,TITLE,figsavetype);
orient(f1,'landscape');
print(f1,file_name,'-dpdf','-vector','-bestfit','');