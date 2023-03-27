clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ = [9 1];
N = 100; %Downsample N times
minArrayLen = 72; %timeWindow = 5 - 72, timeWindow = 10 - 123
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
zMin = -10;
zMax = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prompt the user to select a .mat file using uigetfile
[file, path] = uigetfile('*.mat', 'Select .mat file to load');
if file == 0
    disp('Select a .mat file to start')
    return
end
% If the user didn't click cancel, load the selected file into the workspace
if ischar(file)
    % Combine the file and path to create the full file path
    fullFilePath = fullfile(path, file);
    
    % Load the .mat file into the MATLAB workspace
    load(fullFilePath);
    
    % Display a message to confirm that the file was loaded
    disp(['Loaded ' file]);
end

if isfield(data.streams, 'x405A')
    ISOS = 'x405A';
    GRABDA = 'x465A';
    %time array used for all streams%
    time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    t = 0; % time threshold below which we will discard
    ind = find(time1>t,1);% find first index of when time crosses threshold
    time1 = time1(ind:end); % reformat vector to only include allowed time
    data.streams.(GRABDA).data = data.streams.(GRABDA).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
    
    %downsample streams and time array by N times%
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(GRABDA).data = downsample(data.streams.(GRABDA).data, N);
    minLength = min(length(data.streams.(ISOS).data),length(data.streams.(GRABDA).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
    data.streams.(GRABDA).data = data.streams.(GRABDA).data(1:minLength);
    time1 = downsample(time1, N);
    ts1 = -baseline + (1:minArrayLen) / data.streams.(GRABDA).fs*N;
    
    %detrend & dFF%
    bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
    Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
    Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    detrend_465 = detrend(dFF);
    z465 = zscore(data.streams.(GRABDA).data);
    cueTSA = data.epocs.St1_.onset;
    correct_rewardedA = data.epocs.cRewA.onset;
    correct_norewardA = data.epocs.cNoRewA.onset;
    incorrect_rewardedA = data.epocs.iRewA.onset;
    incorrect_norewardA = data.epocs.iNoRewA.onset;
    [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cueTSA,correct_rewardedA,...
            correct_norewardA,incorrect_rewardedA,incorrect_norewardA);
    for e = 1:height(session_ts)
        e1 = session_ts(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        e3 = session_ts(e,1)-baselineZ(:,2);
        e4 = session_ts(e,1)-baselineZ(:,1);
        [~,ind1] = min(abs(time1-e1));
        [~,ind2] = min(abs(time1-e2));
        [~,ind3] = min(abs(time1-e3));
        [~,ind4] = min(abs(time1-e4));
        GRABDA_zBase = detrend_465(1,ind4:ind3);
        GRABDA_signal = detrend_465(1,ind1:ind2);
        GRABDA_time = time1(1,ind1:ind2);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_signal - zb)/zsd;
        if length(zfinal) < minArrayLen
            continue
        end
        sessionSTREAM(e,:) = zfinal(1:minArrayLen);
        sessionAMP(e,:) = max(zfinal);
        sessionAUC(e,:) = trapz(GRABDA_time,zfinal);
    end
    for e = 1:height(cueTSA)
            e1 = cueTSA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = cueTSA(e,1)-baselineZ(:,2);
            e4 = cueTSA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            GRABDA_zBase = detrend_465(1,ind4:ind3);
            GRABDA_signal = detrend_465(1,ind1:ind2);
            GRABDA_time = time1(1,ind1:ind2);
            zb = mean(GRABDA_zBase);
            zsd = std(GRABDA_zBase);
            zfinal = (GRABDA_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                break
            end
            cueSTREAM(e,:) = zfinal(1:minArrayLen);
            cueAMP(e,:) = max(zfinal);
            cueAUC(e,:) = trapz(GRABDA_time,zfinal);
    end
else
    ISOS = 'x405C';
    GRABDA = 'x465C';
    %time array used for all streams%
    time1 = (1:length(data.streams.(GRABDA).data))/data.streams.(GRABDA).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    t = 0;% time threshold below which we will discard
    ind = find(time1>t,1);% find first index of when time crosses threshold
    time1 = time1(ind:end); % reformat vector to only include allowed time
    data.streams.(GRABDA).data = data.streams.(GRABDA).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
    minLength = min(length(data.streams.(ISOS).data),length(data.streams.(GRABDA).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
    data.streams.(GRABDA).data = data.streams.(GRABDA).data(1:minLength);
    
    %downsample streams and time array by N times%
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(GRABDA).data = downsample(data.streams.(GRABDA).data, N);
    time1 = downsample(time1, N);
    ts1 = -baseline + (1:minArrayLen) / data.streams.(GRABDA).fs*N;
    
    %detrend & dFF%
    bls = polyfit(data.streams.(ISOS).data,data.streams.(GRABDA).data,1);
    Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
    Y_dF_all = data.streams.(GRABDA).data - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    detrend_465 = detrend(dFF);
    z465 = zscore(data.streams.(GRABDA).data);
    cueTSC = data.epocs.St2_.onset;
    correct_rewardedC = data.epocs.cRewC.onset;
    correct_norewardC = data.epocs.cNoRewC.onset;
    incorrect_rewardedC = data.epocs.iRewC.onset;
    incorrect_norewardC = data.epocs.iNoRewC.onset;
    [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cueTSC,correct_rewardedC,...
            correct_norewardC,incorrect_rewardedC,incorrect_norewardC);
    
    for e = 1:height(session_ts)
        e1 = session_ts(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        e3 = session_ts(e,1)-baselineZ(:,2);
        e4 = session_ts(e,1)-baselineZ(:,1);
        [~,ind1] = min(abs(time1-e1));
        [~,ind2] = min(abs(time1-e2));
        [~,ind3] = min(abs(time1-e3));
        [~,ind4] = min(abs(time1-e4));
        GRABDA_zBase = detrend_465(1,ind4:ind3);
        GRABDA_signal = detrend_465(1,ind1:ind2);
        GRABDA_time = time1(1,ind1:ind2);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_signal - zb)/zsd;
        if length(zfinal) < minArrayLen
            continue
        end
        sessionSTREAM(e,:) = zfinal(1:minArrayLen);
        sessionAMP(e,:) = max(zfinal);
        sessionAUC(e,:) = trapz(GRABDA_time,zfinal);
    end
    for e = 1:height(cueTSC)
        e1 = cueTSC(e,1)-baseline;
        e2 = e1+timeWindow+baseline;
        e3 = cueTSC(e,1)-baselineZ(:,2);
        e4 = cueTSC(e,1)-baselineZ(:,1);
        [c1,ind1] = min(abs(time1-e1));
        [c2,ind2] = min(abs(time1-e2));
        [c3,ind3] = min(abs(time1-e3));
        [c4,ind4] = min(abs(time1-e4));
        GRABDA_zBase = detrend_465(1,ind4:ind3);
        GRABDA_signal = detrend_465(1,ind1:ind2);
        GRABDA_time = time1(1,ind1:ind2);
        zb = mean(GRABDA_zBase);
        zsd = std(GRABDA_zBase);
        zfinal = (GRABDA_signal - zb)/zsd;
        if length(zfinal) < minArrayLen
            break
        end
        cueSTREAM(e,:) = zfinal(1:minArrayLen);
        cueAMP(e,:) = max(zfinal);
        cueAUC(e,:) = trapz(GRABDA_time,zfinal);
    end
end
% Generate a colormap for the heatmap
myColorMap = turbo;
% Set a master font size
myFontSize = 14;
session_time = ts1(:,1:minArrayLen);
% Define the number of signals and time steps
numSignals = size(cueSTREAM, 1);


% Find the maximum value of each signal
maxValues = max(abs(cueSTREAM), [], 2);

% Compute the y-axis locations for each signal based on its max value
yLocs = cumsum([0; maxValues+2]);
yLocs = yLocs(1:end-1);

% Create a figure
f1 = figure;
subplot(1,2,1);
% Loop through each signal and plot it on the graph
for i = 0:numSignals-1
    % Plot the signal on the graph at its computed y-axis location
    plot(session_time, cueSTREAM(numSignals-i,:) + yLocs(numSignals-i), 'LineWidth', 1, 'Color', 'b', 'LineStyle', '-');
    hold on;
    
    % Add a label for this signal
    text(timeWindow+0.1, yLocs(numSignals-i) + maxValues(numSignals-i)/2, trial_name(numSignals-i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', myFontSize-1);
    
end

% Set the x-axis limits and label
xlim([session_time(:,1) session_time(:,minArrayLen)]);
xlabel('Time (s)','FontSize',myFontSize+2);
xline(0,'LineWidth',2,'Color','black')
% Set the y-axis limits and label
ylim([0 max(yLocs)+maxValues(end)]);
ylabel('GrabDA (Z-Score)','FontSize',myFontSize+2);
set(gca,'FontSize',myFontSize);
% Add a title to the graph
[~,name,~] = fileparts(file);
name = strrep(name,'_',' ');
title(name,'FontSize',myFontSize+3);


% Create a figure
subplot(1,2,2);

signals = cueSTREAM; 
x_values = lever_ts; 

% Define the time array
timeArray = session_time;

% Plot the heatmap
imagesc(timeArray, 1:size(signals, 1), cueSTREAM);
set(gca,'YDir','normal','FontSize',myFontSize);
% Set the colormap for the heatmap
colormap(myColorMap);

% Add a colorbar to the heatmap
colorbar('FontSize',myFontSize);
caxis([zMin, zMax]);
% Loop through each row of signals
% Get the y-limits and y-axis unit size of the heatmap
ylims = ylim;
y_unit = diff(ylims) / size(signals, 1);
for k = 1:size(signals, 1)
    % Calculate the y position of the indicator
    y_pos = ylims(1) + (k - 1) * y_unit;
    
    % Get the x and y data for the indicator
    x_data = x_values(k) * ones(2, 1);
    y_data = [y_pos, y_pos + y_unit];
    
    % Add the indicator using the line function
    line(x_data, y_data, 'Color', 'white', 'LineWidth', 3);
end
% Set the x-axis label
xlabel('Time (s)','FontSize',myFontSize+2);
xline(0,'LineWidth',2,'Color','black')
% Set the y-axis label
ylabel('Trial', 'FontSize',myFontSize+2);

% Add a title to the heatmap
title('Session Heatmap','FontSize',myFontSize+3);
set(f1,'Position',[200 200 1600 800]);