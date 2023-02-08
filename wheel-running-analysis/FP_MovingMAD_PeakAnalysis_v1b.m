%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_MovingMAD_PeakAnalysis.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)

VERSION = "1.0b";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');


session_duration = 3600;
disp(VERSION)
myDir = uigetdir;
BLOCKPATH = myDir;
data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});

%Stream Stores%
DLS_ISOS = 'x405A'; % name of the 405A store
DLS_DA = 'x465A'; % name of the 465A store
NAc_ISOS = 'x405C'; % name of the 405C store
NAc_DA = 'x465C'; % name of the 465C store
N = 100; %Downsample N times
%time array used for all streams%
time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 20; % time threshold below which we will discard
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
window_size = 600;
peak_indicies_DLS = [];
peak_indicies_NAc = [];
for ii = 1:length(detrend_465A)
    current_window_DLS = detrend_465A(max(1,ii-window_size):min(length(detrend_465A),...
        ii+window_size));
    DLS_MOVMAD(ii,:) = median(abs(current_window_DLS - median(current_window_DLS)));
    threshold_DLS = DLS_MOVMAD;
    if detrend_465A(ii) > threshold_DLS
        peak_indicies_DLS = [peak_indicies_DLS,ii];
    end
    current_window_NAc = detrend_465C(max(1,ii-window_size):min(length(detrend_465C),...
        ii+window_size));
    NAc_MOVMAD(ii,:) = median(abs(current_window_NAc - median(current_window_NAc)));
    threshold_NAc = NAc_MOVMAD;
    if detrend_465C(ii) > threshold_NAc
        peak_indicies_NAc = [peak_indicies_NAc,ii];
    end
end
