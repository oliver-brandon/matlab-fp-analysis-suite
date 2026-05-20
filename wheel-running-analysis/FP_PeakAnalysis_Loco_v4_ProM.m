%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_PeakAnalysis_Loco.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)

VERSION = '4.0';

%Input = Locomotor DA .mat files
%Output = table called "loco_peak_analysis" that includes number of peaks,
%peaks/s, max amp, avg amp, avg prominence, and avg width for DLS and NAc
%signals

%Instructions: Point MATLAB at the desired .mat file folder (data), set
%session time in seconds (session_duration), and click run. Data is output
%as a table called "loco_peak_analysis."

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%House Keeping%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');

%Choose batch or single file analysis%
tanks2analyze = 2; %1 = batch analyze, 2 = single file analysis%
%Peak detection variables
fS = 1017.3;
minPkDist = round(fS/4);
minPkHeight = 1;
minPkProm = 2;
zPkThresh = 1;
lLim = 1;
uLim = 1;
movAvgWinSec = 50;
%Choose session duration%
session_duration = 3600; %duration of recording in seconds

if tanks2analyze == 1
    myDir = uigetdir(pwd, 'Select folder containing .mat files'); %gets directory%
    if isequal(myDir,0)
        error('No folder selected.');
    end
    myFiles = dir(fullfile(myDir, '*.mat')); %gets all .mat files in directory%
    numFiles = length(myFiles);
    if numFiles == 0
        error('No .mat files found in the selected directory.');
    end
    disp("Starting batch analysis...")
    
    DLS_pks = nan(numFiles,1);
    DLS_pk_min = nan(numFiles,1);
    DLS_amp_max = nan(numFiles,1);
    DLS_amp_avg = nan(numFiles,1);
    DLS_amp_prom = nan(numFiles,1);
    DLS_amp_width = nan(numFiles,1);
    NAc_pks = nan(numFiles,1);
    NAc_pk_min = nan(numFiles,1);
    NAc_amp_max = nan(numFiles,1);
    NAc_amp_avg = nan(numFiles,1);
    NAc_amp_prom = nan(numFiles,1);
    NAc_amp_width = nan(numFiles,1);
    file_name = strings(numFiles,1);
    
    for i = 1:numFiles
        fprintf('Loading file %d of %d...\n',i,numFiles)
        BLOCKPATH = fullfile(myDir, myFiles(i).name);
        file_name(i) = string(myFiles(i).name);
        load(BLOCKPATH);
        
        %Stream Stores%
        DLS_ISOS = 'x405A'; % name of the 405A store
        DLS_DA = 'x465A'; % name of the 465A store
        NAc_ISOS = 'x405C'; % name of the 405C store
        NAc_DA = 'x465C'; % name of the 465C store
        N = 1; %Downsample factor (1 = no downsampling)
        
        [time, DLS_out, NAc_out] = process_tank(data, DLS_ISOS, DLS_DA, NAc_ISOS, NAc_DA, ...
            session_duration, N, minPkDist, minPkHeight, minPkProm, zPkThresh, movAvgWinSec);
        
        DLS_pks(i) = DLS_out.numPeaks;
        DLS_pk_min(i) = DLS_out.peaksPerMin;
        DLS_amp_max(i) = DLS_out.maxAmp;
        DLS_amp_avg(i) = DLS_out.avgAmp;
        DLS_amp_prom(i) = DLS_out.avgProm;
        DLS_amp_width(i) = DLS_out.avgWidth;
        
        NAc_pks(i) = NAc_out.numPeaks;
        NAc_pk_min(i) = NAc_out.peaksPerMin;
        NAc_amp_max(i) = NAc_out.maxAmp;
        NAc_amp_avg(i) = NAc_out.avgAmp;
        NAc_amp_prom(i) = NAc_out.avgProm;
        NAc_amp_width(i) = NAc_out.avgWidth;
        
        plot_session_traces(time, DLS_out, NAc_out, file_name(i), zPkThresh, lLim, uLim);
    end
elseif tanks2analyze == 2
    disp("Starting single file analysis...")
    [myFile, myPath] = uigetfile('*.mat', 'Select .mat file');
    if isequal(myFile,0)
        error('No file selected.');
    end
    BLOCKPATH = fullfile(myPath, myFile);
    numFiles = 1;
    load(BLOCKPATH);
    
    %Stream Stores%
    DLS_ISOS = 'x405A'; % name of the 405A store
    DLS_DA = 'x465A'; % name of the 465A store
    NAc_ISOS = 'x405C'; % name of the 405C store
    NAc_DA = 'x465C'; % name of the 465C store
    N = 1; %Downsample factor%
    
    [time, DLS_out, NAc_out] = process_tank(data, DLS_ISOS, DLS_DA, NAc_ISOS, NAc_DA, ...
        session_duration, N, minPkDist, minPkHeight, minPkProm, zPkThresh, movAvgWinSec);
    
    DLS_pks = DLS_out.numPeaks;
    DLS_pk_min = DLS_out.peaksPerMin;
    DLS_amp_max = DLS_out.maxAmp;
    DLS_amp_avg = DLS_out.avgAmp;
    DLS_amp_prom = DLS_out.avgProm;
    DLS_amp_width = DLS_out.avgWidth;
    
    NAc_pks = NAc_out.numPeaks;
    NAc_pk_min = NAc_out.peaksPerMin;
    NAc_amp_max = NAc_out.maxAmp;
    NAc_amp_avg = NAc_out.avgAmp;
    NAc_amp_prom = NAc_out.avgProm;
    NAc_amp_width = NAc_out.avgWidth;
    
    DLS_peak_trace = DLS_out.zFilt;
    NAc_peak_trace = NAc_out.zFilt;
    Time_trace = time;
    DLS_peak_timestamps = DLS_out.finalPeakTimes;
    NAc_peak_timestamps = NAc_out.finalPeakTimes;
    
    plot_session_traces(time, DLS_out, NAc_out, string(myFile), zPkThresh, lLim, uLim);
end

if tanks2analyze == 1
    loco_peak_analysis = table(file_name, DLS_pks, DLS_pk_min, DLS_amp_max, DLS_amp_avg, ...
        DLS_amp_prom, DLS_amp_width, NAc_pks, NAc_pk_min, NAc_amp_max, NAc_amp_avg, ...
        NAc_amp_prom, NAc_amp_width, ...
        'VariableNames', {'File','DLS Peaks','DLS Peaks/m','DLS Max Amp',...
        'DLS Avg Amp','DLS Peak Prom','DLS Peak Width','NAc Peaks','NAc Peaks/m','NAc Max Amp',...
        'NAc Avg Amp','NAc Peak Prom','NAc Peak Width'});
else
    loco_peak_analysis = table(DLS_pks, DLS_pk_min, DLS_amp_max, DLS_amp_avg, ...
        DLS_amp_prom, DLS_amp_width, NAc_pks, NAc_pk_min, NAc_amp_max, NAc_amp_avg, ...
        NAc_amp_prom, NAc_amp_width, ...
        'VariableNames', {'DLS Peaks','DLS Peaks/m','DLS Max Amp',...
        'DLS Avg Amp','DLS Peak Prom','DLS Peak Width','NAc Peaks','NAc Peaks/m','NAc Max Amp',...
        'NAc Avg Amp','NAc Peak Prom','NAc Peak Width'});
end

%UITable (figure) that displays "loco_peak_analysis" table%
figure;
uitable('Data',loco_peak_analysis{:,:},'ColumnName',loco_peak_analysis.Properties.VariableNames,...
    'RowName',loco_peak_analysis.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

disp("Analysis complete!")
fprintf('Successfully analyzed %d file(s)!',numFiles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time, chanA, chanC] = process_tank(data, DLS_ISOS, DLS_DA, NAc_ISOS, NAc_DA, ...
    session_duration, N, minPkDist, minPkHeight, minPkProm, zPkThresh, movAvgWinSec)

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
    
    %downsample streams and time array by N factor%
    if N > 1
        data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
        data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
        data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
        data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
        time = downsample(time, N);
    end
    
    fs_ds = data.streams.(DLS_DA).fs / N;
    minPkDist_ds = max(1, round(minPkDist / N));
    movAvgWinSamp = max(1, round(movAvgWinSec * fs_ds));
    
    %465A processing%
    chanA = process_channel(data.streams.(DLS_ISOS).data, data.streams.(DLS_DA).data, ...
        time, fs_ds, session_duration, movAvgWinSamp, minPkDist_ds, minPkHeight, minPkProm, zPkThresh);
    
    %465C processing%
    chanC = process_channel(data.streams.(NAc_ISOS).data, data.streams.(NAc_DA).data, ...
        time, fs_ds, session_duration, movAvgWinSamp, minPkDist_ds, minPkHeight, minPkProm, zPkThresh);
end

function out = process_channel(isosData, expData, time, fs_ds, session_duration, movAvgWinSamp, ...
    minPkDist_ds, minPkHeight, minPkProm, zPkThresh)

    % Step 1: correct experimental channel using isosbestic channel, then
    % detrend the resulting %dF/F trace
    bls = polyfit(isosData, expData, 1);
    Y_fit_all = bls(1) .* isosData + bls(2);
    Y_dF_all = expData - Y_fit_all;
    raw_dFF = 100 * (Y_dF_all) ./ Y_fit_all;
    raw_dFF = detrend(raw_dFF);
    
    % Step 2: z-score normalization
    zTrace = zscore(raw_dFF);
    
    % Step 3: 100-s moving average filter
    zFilt = movmean(zTrace, movAvgWinSamp);
    
    % Step 4: detect significant peaks on filtered z-scored trace
    [candPks, candLocs, candWid, candProm] = findpeaks(zFilt, ...
        'MinPeakHeight', zPkThresh, 'MinPeakDistance', minPkDist_ds);
    candTimes = time(candLocs);
    
    % Step 5: timestamp candidate peaks
    
    % Steps 6-8: revisit raw %dF/F trace and retain only peaks that also
    % satisfy prominence, height, and distance criteria in the raw trace
    [rawPks, rawLocs, rawWid, rawProm] = findpeaks(raw_dFF, ...
        'MinPeakHeight', minPkHeight, 'MinPeakProminence', minPkProm, ...
        'MinPeakDistance', minPkDist_ds);
    rawTimes = time(rawLocs);
    
    finalMask = false(size(candLocs));
    finalRawIdx = nan(size(candLocs));
    matchTol = max(1, floor(minPkDist_ds/2));
    
    for j = 1:numel(candLocs)
        nearby = find(abs(rawLocs - candLocs(j)) <= matchTol);
        if ~isempty(nearby)
            [~, bestIdx] = max(rawProm(nearby));
            finalMask(j) = true;
            finalRawIdx(j) = nearby(bestIdx);
        end
    end
    
    finalRawIdx = unique(finalRawIdx(finalMask));
    finalRawIdx = finalRawIdx(~isnan(finalRawIdx));
    
    finalPeakVals = rawPks(finalRawIdx);
    finalPeakLocs = rawLocs(finalRawIdx);
    finalPeakTimes = rawTimes(finalRawIdx);
    finalPeakProm = rawProm(finalRawIdx);
    finalPeakWidSamples = rawWid(finalRawIdx);
    finalPeakWidSec = finalPeakWidSamples ./ fs_ds;
    
    if isempty(finalPeakVals)
        numPeaks = 0;
        peaksPerMin = 0;
        maxAmp = NaN;
        avgAmp = NaN;
        avgProm = NaN;
        avgWidth = NaN;
    else
        numPeaks = numel(finalPeakVals);
        peaksPerMin = (numPeaks / session_duration) * 60;
        maxAmp = max(finalPeakVals);
        avgAmp = mean(finalPeakVals);
        avgProm = mean(finalPeakProm);
        avgWidth = mean(finalPeakWidSec);
    end
    
    out = struct();
    out.rawDFF = raw_dFF;
    out.zTrace = zTrace;
    out.zFilt = zFilt;
    out.candPks = candPks;
    out.candLocs = candLocs;
    out.candTimes = candTimes;
    out.candWid = candWid;
    out.candProm = candProm;
    out.rawPks = rawPks;
    out.rawLocs = rawLocs;
    out.rawTimes = rawTimes;
    out.rawWid = rawWid;
    out.rawProm = rawProm;
    out.finalPeakVals = finalPeakVals;
    out.finalPeakLocs = finalPeakLocs;
    out.finalPeakTimes = finalPeakTimes;
    out.finalPeakProm = finalPeakProm;
    out.finalPeakWidSec = finalPeakWidSec;
    out.numPeaks = numPeaks;
    out.peaksPerMin = peaksPerMin;
    out.maxAmp = maxAmp;
    out.avgAmp = avgAmp;
    out.avgProm = avgProm;
    out.avgWidth = avgWidth;
end

function plot_session_traces(time, DLS_out, NAc_out, plotLabel, zPkThresh, lLim, uLim)
    figure('Name', char(plotLabel), 'Color', 'w');
    tiledlayout(2,1);
    
    % DLS plot
    nexttile;
    plot(time, DLS_out.zFilt, 'b'); hold on;
    if ~isempty(DLS_out.candLocs)
        plot(DLS_out.candTimes, DLS_out.candPks, 'rv', 'MarkerFaceColor', 'r');
    end
    if ~isempty(DLS_out.finalPeakLocs)
        plot(DLS_out.finalPeakTimes, DLS_out.zFilt(DLS_out.finalPeakLocs), 'ko', 'MarkerSize', 5);
    end
    yline(zPkThresh, '--k');
    xlim([time(1) time(end)]);
    if lLim ~= 1 || uLim ~= 1
        ylim([lLim uLim]);
    end
    xlabel('Time (s)');
    ylabel('z-scored 465A');
    title(['DLS (465A) - ' char(plotLabel)]);
    txtA = sprintf('Final peaks/min: %.3f\nAvg prominence: %.3f\nAvg width (s): %.3f', ...
        DLS_out.peaksPerMin, DLS_out.avgProm, DLS_out.avgWidth);
    text(0.01, 0.98, txtA, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold off;
    
    % NAc plot
    nexttile;
    plot(time, NAc_out.zFilt, 'm'); hold on;
    if ~isempty(NAc_out.candLocs)
        plot(NAc_out.candTimes, NAc_out.candPks, 'rv', 'MarkerFaceColor', 'r');
    end
    if ~isempty(NAc_out.finalPeakLocs)
        plot(NAc_out.finalPeakTimes, NAc_out.zFilt(NAc_out.finalPeakLocs), 'ko', 'MarkerSize', 5);
    end
    yline(zPkThresh, '--k');
    xlim([time(1) time(end)]);
    if lLim ~= 1 || uLim ~= 1
        ylim([lLim uLim]);
    end
    xlabel('Time (s)');
    ylabel('z-scored 465C');
    title(['NAc (465C) - ' char(plotLabel)]);
    txtC = sprintf('Final peaks/min: %.3f\nAvg prominence: %.3f\nAvg width (s): %.3f', ...
        NAc_out.peaksPerMin, NAc_out.avgProm, NAc_out.avgWidth);
    text(0.01, 0.98, txtC, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold off;
end
