%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%DA_Wheel_PeakAnalysis_Loco.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = 2.0;

%Input = TDT fiber photometry tank with epocs stored in the "Cam1" notes 
%Output = uitable called "master_peak_analysis" that includes number of peaks,
%peaks/s for DLS and NAc signals during running bouts and during
%non-running bouts.

%Instructions: Point MATLAB at the desired TDT tank folder (data) and click run. 
%Data is output as a table called "master_peak_analysis." There is also a 
%"session_avg" table that displays the average peaks, peaks/s, max peak amp,
% and avg peak amp for the session. 
%This script uses manual TTLs from the Cam1 epoch to create "runBout" time
%stamps. MAD1 = MAD value for the 465A and MAD2 = MAD value for the 465C


myDir = uigetdir; %gets directory%
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~ismember({myFiles.name},{'.','..'}));
numFiles = length(myFiles);
master_bout_analysis = cell(1,length(myFiles));
session_avg = cell(1,length(myFiles));
fprintf("VERSION: %d\n",VERSION)
disp("Starting batch analysis...")
for i = 1:length(myFiles)
    BLOCKPATH = fullfile(myDir, myFiles(i).name);
    [tankpath,tankname,ext] = fileparts(BLOCKPATH);
    fprintf('Loading tank %d of %d...\n',i,length(myFiles))
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
    
    %% makes new wheel data epocs %%
    %% epocs created in TDT OpenScope save in the notes of Cam1 %% 
    %combine index with timestamp data from Cam1 notes%
    ind = double(data.epocs.Cam1.notes.index);
    ts = data.epocs.Cam1.notes.ts;
    var1 = [ind ts];
    %separate by index to make separate epocs%
    % onWheel = var1(ismember(var1(:,1),[1]),:);
    % offWheel = var1(ismember(var1(:,1),[2]),:);
    runStart = var1(ismember(var1(:,1),[3]),:);
    runStop = var1(ismember(var1(:,1),[4]),:);
    %extract time stamps%
    % onWheelTs = onWheel(:,2);
    % offWheelTs = offWheel(:,2);
    runStartTs = runStart(:,2);
    runStopTs = runStop(:,2);
    %extract indicies%
    % onWheelInd = onWheel(:,1);
    % offWheelInd = offWheel(:,1);
    runStartInd = runStart(:,1);
    runStopInd = runStop(:,1);
    %runStart%
    data.epocs.runStart.onset = runStartTs;
    data.epocs.runStart.offset = runStartTs + 0.01;
    data.epocs.runStart.name = 'runStart';
    data.epocs.runStart.data = runStartInd;
    %runStop%
    data.epocs.runStop.onset = runStopTs;
    data.epocs.runStop.offset = runStopTs + 0.01;
    data.epocs.runStop.name = 'runStop';
    data.epocs.runStop.data = runStopInd;

    
    %creates running bout epoc from runStart and runStop onsets%
    data.epocs.runBout.onset = runStartTs;
    data.epocs.runBout.offset = runStopTs;
    data.epocs.runBout.name = 'runBout';
    data.epocs.runBout.typeStr = data.epocs.Cam1.typeStr;
    data.epocs.runBout.data = ones(length(data.epocs.runBout.onset),1)*5;
    %Finds the running bout time stamps%
    runBout_ts = [data.epocs.runBout.onset data.epocs.runBout.offset]; 
    end_ts = height(runBout_ts);
    
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
    min_time = min(time);
    
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
    
    %465C%
    bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
    Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
    Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
    dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
    std_dFF2 = std(double(dFF2));
    detrend_465C = detrend(dFF2);
    
    %calculates and plots median absolute deviation for both 465 signals%
    MAD1 = mad(detrend_465A, 1);
    MAD2 = mad(detrend_465C, 1);
    [pks,locs,w,p] = findpeaks(detrend_465A, time, 'MinPeakHeight', MAD1);
    [pks2,locs2,w2,p2] = findpeaks(detrend_465C, time, 'MinPeakHeight', MAD2);
    
    %plots entire streams with peak indicators% 
%     subplot(2,1,1);
%     findpeaks(detrend_465A, time, 'MinPeakHeight', MAD1);
%     subplot(2,1,2)
%     findpeaks(detrend_465C, time, 'MinPeakHeight', MAD2);
    
    %establish time stamps for non-running "bouts"%
    noRun_ts = [min_time runBout_ts(1,1)];
    for h = 1:end_ts
        h1 = runBout_ts(h,2);
        h2 = runBout_ts(h+1,1);
        noRun_ts(h+1,:) = [h1 h2];
        if h+1 == end_ts
            break
        end
    end
    
    %Peak analysis for non-running bouts%
    for y = 1:end_ts
        y1 = noRun_ts(y,1);
        y2 = noRun_ts(y,2);
        [c, index] = min(abs(time-y1));
        [c2, index2] = min(abs(time-y2)); 
        DOPE1_sig = detrend_465A(1,index:index2);
        DOPE2_sig = detrend_465C(1,index:index2);
        DOPE_time = time(1,index:index2);
      
        [pks3,locs3,w3,p3] = findpeaks(DOPE1_sig, DOPE_time, 'MinPeakHeight', MAD1);
        [pks4,locs4,w4,p4] = findpeaks(DOPE2_sig, DOPE_time, 'MinPeakHeight', MAD2);
        if isempty(pks3)
            pks3 = 0;
        end
        if isempty(pks4)
            pks4 = 0;
        end
        DLS_numpeak = length(pks3);
        DLS_pk_min = (DLS_numpeak/(y2-y1)*60);
        DLS_max_amp = max(pks3);
        DLS_avg_amp = mean(pks3);
        NAc_numpeak = length(pks4);
        NAc_pk_min = (NAc_numpeak/(y2-y1)*60);
        NAc_max_amp = max(pks4);
        NAc_avg_amp = mean(pks4);
        noRun_peak_analysis(y,:) = table(y1, DLS_numpeak, DLS_pk_min, DLS_max_amp, DLS_avg_amp, NAc_numpeak, ...
            NAc_pk_min, NAc_max_amp, NAc_avg_amp, 'VariableNames', {'No Run Ts', 'DLS Peaks Stop', 'DLS Peaks/m Stop', ...
            'DLS Max Amp Stop', 'DLS Avg Amp Stop', 'NAc Peaks Stop', 'NAc Peaks/m Stop', 'NAc Max Amp Stop', ...
            'NAc Avg Amp Stop'});
    end
    
    %Peak analysis for running bouts%
    for x = 1:end_ts
        x1 = runBout_ts(x,1);
        x2 = runBout_ts(x,2);
        [c, index] = min(abs(time-x1));
        [c2, index2] = min(abs(time-x2)); 
        DOPE1_sig = detrend_465A(1,index:index2);
        DOPE2_sig = detrend_465C(1,index:index2);
        DOPE_time = time(1,index:index2);
      
        [pks3,locs3,w3,p3] = findpeaks(DOPE1_sig, DOPE_time, 'MinPeakHeight', MAD1);
        [pks4,locs4,w4,p4] = findpeaks(DOPE2_sig, DOPE_time, 'MinPeakHeight', MAD2);
        if isempty(pks3)
            pks3 = 0;
        end
        if isempty(pks4)
            pks4 = 0;
        end
        DLS_numpeak = length(pks3);
        DLS_pk_min = (DLS_numpeak/(x2-x1)*60);
        DLS_max_amp = max(pks3);
        DLS_avg_amp = mean(pks3);
        NAc_numpeak = length(pks4);
        NAc_pk_min = (NAc_numpeak/(x2-x1)*60);
        NAc_max_amp = max(pks4);
        NAc_avg_amp = mean(pks4);
        run_peak_analysis(x,:) = table(x1, DLS_numpeak, DLS_pk_min, DLS_max_amp, DLS_avg_amp, NAc_numpeak, ...
            NAc_pk_min, NAc_max_amp, NAc_avg_amp, 'VariableNames', {'Run Start Ts', 'DLS Peaks Run', 'DLS Peaks/m Run', ...
            'DLS Max Amp Run', 'DLS Avg Amp Run', 'NAc Peaks Run', 'NAc Peaks/m Run', ...
            'NAc Max Amp Run', 'NAc Avg Amp Run'});
    end
    %combines peak analysis for run bouts and non run bouts into one table%
        master_bout_analysis{i} = [run_peak_analysis noRun_peak_analysis];
        session_avg{i} = varfun(@mean, master_bout_analysis{i});
end

        



%UITable (figure) that displays "master_peak_analysis" table%
% g = figure;
% v = uipanel(g);
% panel_loc = 0;
% for ii = 1:numel(master_bout_analysis)
%     uitable(v,'Data',master_bout_analysis{1,ii},'ColumnName',master_bout_analysis{1,ii}.Properties.VariableNames,...
%          'RowName',master_bout_analysis{1,ii}.Properties.RowNames,'Units', 'Normalized', 'Position',[0, panel_loc, 1, 1]);
%     panel_loc = panel_loc + 1;
% end



% %gathers session averages for master_peak_analysis (mean of each column)%
% 
% figure;
% uitable('Data',session_avg{:,:},'ColumnName',session_avg.Properties.VariableNames,...
%     'RowName',session_avg.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    




