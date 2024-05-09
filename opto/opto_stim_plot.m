clear; clc; close all;

t = 10;
N = 10;
channel = 1;


BLOCKPATH = 'E:\Google Drive\optomouse-prime\test-tanks\1034F_20Hz-2mW-5pulse';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});

if channel == 1
    ISOS = 'x405A';
    Grab = 'x465A';
elseif channel == 2
    ISOS = 'x405C';
    Grab = 'x465C';
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



% dF/F
bls = polyfit(ISOS_raw, Grab_raw, 1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = Grab_raw - Y_fit_all; %dF (units mV) is not dFF
Grab_dFF = 100*(Y_dF_all)./Y_fit_all;

% Photobleach correction
ISOS_raw = detrend(ISOS_raw);
Grab_dFF = detrend(Grab_dFF);

% noise reduction using moving median
Grab_filt = smoothdata(Grab_dFF,'movmedian',100);

plot(time, Grab_filt);