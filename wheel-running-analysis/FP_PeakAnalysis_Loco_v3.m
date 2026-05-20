%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_PeakAnalysis_Loco.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)

VERSION = "3.0";

%Input = Locomotor DA .mat files
%Output = table called "loco_peak_analysis" that includes number of peaks,
%peaks/s, max amp, and avg amp for DLS and NAc signals

%Instructions: Point MATLAB at the desired .mat file folder (data), set
%session time in seconds (session_duration), and click run. Data is output
%as a table called "peak_analysis."


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%House Keeping%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');

%Choose batch or single file analysis%
tanks2analyze = 1; %1 = batch analyze, 2 = single file analysis%
%Peak detection variables
fS = 1017.3;
minPkDist = round(fS/4);
minPkHeight = 1;
minPkProm = 2;
zPkThresh = 2.6;
lLim = 1;
uLim = 1;
%Choose session duration%
session_duration = 3600; %duration of recording in seconds

disp(VERSION);
if tanks2analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~ismember({myFiles.name},{'.','..'}));
    disp("Starting batch analysis...")
    for i = 1:length(myFiles)
        fprintf('Loading tank %d of %d...\n',i,length(myFiles))
        BLOCKPATH = fullfile(myDir, myFiles(i).name);
        data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});
        
        %Stream Stores%
        DLS_ISOS = 'x405A'; % name of the 405A store
        DLS_DA = 'x465A'; % name of the 465A store
        NAc_ISOS = 'x405C'; % name of the 405C store
        NAc_DA = 'x465C'; % name of the 465C store
        N = 0; %Downsample N times
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
        detrend_465A_dFF = detrend(dFF);
        detrend_465A_Z = zscore(detrend_465A_dFF);
        %465C%
        bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
        Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
        Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
        dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
        std_dFF2 = std(double(dFF2));
        detrend_465C_dFF = detrend(dFF2);
        detrend_465C_Z = zscore(detrend_465C_dFF);
        
        % locate peaks using the z-scored signals
        [pks,locs,wid,prm] = findpeaks(detrend_465A_Z, time, 'MinPeakHeight', minPkHeight, 'MinPeakProminence', zPkThresh,...
            "MinPeakDistance", minPkDist);
        [pks2,locs2,wid2,prm2] = findpeaks(detrend_465C_Z, time, 'MinPeakHeight', minPkHeight, 'MinPeakProminence', zPkThresh,...
            "MinPeakDistance", minPkDist);
        
        DLS_pks(i,:) = length(pks);
        DLS_pk_min(i,:) = (DLS_pks(i)/session_duration)*60;
        DLS_amp_max(i,:) = max(pks);
        DLS_amp_avg(i,:) = mean(pks);
        DLS_amp_prom(i,:) = mean(prm);
        NAc_pks(i,:) = length(pks2);
        NAc_pk_min(i,:) = (NAc_pks(i)/session_duration)*60;
        NAc_amp_max(i,:) = max(pks2);
        NAc_amp_avg(i,:) = mean(pks2);
        NAc_amp_prom(i,:) = mean(prm2);

    end
elseif tanks2analyze == 2
    disp("Starting single tank analysis...")
    myDir = uigetdir;
    BLOCKPATH = myDir;
    numFiles = 1;
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
    %calculates and plots median absolute deviation for both 465 signals%
    MAD1 = mad(detrend_465A, 1);
    MAD2 = mad(detrend_465C, 1);
    [pks,locs,wid,prm] = findpeaks(detrend_465A, time, 'MinPeakHeight', MAD1,...
        "MinPeakDistance", minPkDist);
    [pks2,locs2,wid2,prm2] = findpeaks(detrend_465C, time, 'MinPeakHeight', MAD2,...
        "MinPeakDistance", minPkDist);
    
    DLS_pks = length(pks);
    DLS_pk_min = (DLS_pks/session_duration)*60;
    DLS_amp_max = max(pks);
    DLS_amp_avg = mean(pks);
    DLS_amp_prom = mean(prm);
    NAc_pks = length(pks2);
    NAc_pk_min = (NAc_pks/session_duration)*60;
    NAc_amp_max = max(pks2);
    NAc_amp_avg = mean(pks2);
    NAc_amp_prom = mean(prm2);
end


loco_peak_analysis = table(DLS_pks, DLS_pk_min, DLS_amp_max, DLS_amp_avg, ...
DLS_amp_prom, NAc_pks, NAc_pk_min, NAc_amp_max, NAc_amp_avg, NAc_amp_prom,...
'VariableNames', {'DLS Peaks','DLS Peaks/m','DLS Max Amp',...
'DLS Avg Amp','DLS Peak Prom','NAc Peaks','NAc Peaks/m','NAc Max Amp',...
'NAc Avg Amp','NAc Peak Prom'});

%UITable (figure) that displays "peak_analysis" table%
figure;
uitable('Data',loco_peak_analysis{:,:},'ColumnName',loco_peak_analysis.Properties.VariableNames,...
    'RowName',loco_peak_analysis.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

disp("Analysis complete!")
fprintf('Successfully analyzed %d tank(s)!',numFiles)