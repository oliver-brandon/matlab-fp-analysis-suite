clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 15; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [3 1];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
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
    SIGNAL = 'x465A';
    %time array used for all streams%
    time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    ind = find(time>t,1);% find first index of when time crosses threshold
    time = time(ind:end); % reformat vector to only include allowed time
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
    
    %downsample streams and time array by N times%
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
    minLength = min(length(data.streams.(ISOS).data),length(data.streams.(SIGNAL).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(1:minLength);
    time = downsample(time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;
    
    ISOS_raw = data.streams.(ISOS).data;
    SIGNAL_raw = data.streams.(SIGNAL).data;
    
    %detrend & dFF%
    bls = polyfit(ISOS_raw,SIGNAL_raw,1);
    Y_fit_all = bls(1) .* ISOS_raw + bls(2);
    Y_dF_all = SIGNAL_raw - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    dFF = detrend(dFF);
    z465 = zscore(dFF);
    
    cueTSA = data.epocs.St1_.onset;
    correct_rewardedA = data.epocs.cRewA.onset;
    correct_norewardA = data.epocs.cNoRewA.onset;
    incorrect_rewardedA = data.epocs.iRewA.onset;
    incorrect_norewardA = data.epocs.iNoRewA.onset;
    [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cueTSA,correct_rewardedA,...
            correct_norewardA,incorrect_rewardedA,incorrect_norewardA);
    for e = 1:height(session_identifiers(:,1))
        e1 = session_identifiers(e,1)-baseWindow;
        e2 = e1+timeWindow+baseWindow;
        e3 = session_identifiers(e,1)-baseline(:,2);
        e4 = session_identifiers(e,1)-baseline(:,1);
        [~,ind1] = min(abs(time-e1));
        [~,ind2] = min(abs(time-e2));
        [~,ind3] = min(abs(time-e3));
        [~,ind4] = min(abs(time-e4));
        SIGNAL_meanBase = mean(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_std = std(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_window = SIGNAL_raw(1,ind1:ind2);
        SIGNAL_dFF = SIGNAL_window - SIGNAL_meanBase;
        SIGNAL_dFF = 100*(SIGNAL_dFF / SIGNAL_meanBase);
        SIGNAL_time = time(1,ind1:ind2);
        SIGNAL_z = (SIGNAL_window - SIGNAL_meanBase)/SIGNAL_std;
        if length(SIGNAL_z) < epocArrayLen
            continue
        end
        sessionSTREAM(e,:) = SIGNAL_z(1:epocArrayLen);
        sessionAMP(e,:) = max(SIGNAL_z);
        sessionAUC(e,:) = trapz(SIGNAL_time,SIGNAL_z);
    end
    for e = 1:height(cueTSA)
        e1 = cueTSA(e,1)-baseWindow;
        e2 = e1+timeWindow+baseWindow;
        e3 = cueTSA(e,1)-baseline(:,2);
        e4 = cueTSA(e,1)-baseline(:,1);
        [c1,ind1] = min(abs(time-e1));
        [c2,ind2] = min(abs(time-e2));
        [c3,ind3] = min(abs(time-e3));
        [c4,ind4] = min(abs(time-e4));
        SIGNAL_meanBase = mean(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_std = std(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_window = SIGNAL_raw(1,ind1:ind2);
        SIGNAL_dFF = SIGNAL_window - SIGNAL_meanBase;
        SIGNAL_dFF = 100*(SIGNAL_dFF / SIGNAL_meanBase);
        SIGNAL_time = time(1,ind1:ind2);
        SIGNAL_z = (SIGNAL_window - SIGNAL_meanBase)/SIGNAL_std;
        if length(SIGNAL_z) < epocArrayLen
            mn = mean(SIGNAL_z(1,end-10:end));
            SIGNAL_z(1,end:epocArrayLen) = mn;
        elseif length(SIGNAL_z) > epocArrayLen
            op = length(SIGNAL_z);
            arrayDif = op - epocArrayLen;
            SIGNAL_z = SIGNAL_z(1,1:end-arrayDif);
        end
        cueSTREAM(e,:) = SIGNAL_z(1:epocArrayLen);
        cueAMP(e,:) = max(SIGNAL_z);
        % cueAUC(e,:) = trapz(SIGNAL_time,SIGNAL_z);
    end
else
    ISOS = 'x405C';
    SIGNAL = 'x465C';
    %time array used for all streams%
    time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    ind = find(time>t,1);% find first index of when time crosses threshold
    time = time(ind:end); % reformat vector to only include allowed time
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
    
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
    minLength = min(length(data.streams.(ISOS).data),length(data.streams.(SIGNAL).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(1:minLength);
    time = downsample(time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;
    
    ISOS_raw = data.streams.(ISOS).data;
    SIGNAL_raw = data.streams.(SIGNAL).data;

    %detrend & dFF%
    bls = polyfit(ISOS_raw,SIGNAL_raw,1);
    Y_fit_all = bls(1) .* ISOS_raw + bls(2);
    Y_dF_all = SIGNAL_raw - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    dFF = detrend(dFF);
    z465 = zscore(dFF);

    cueTSC = data.epocs.St2_.onset;
    correct_rewardedC = data.epocs.cRewC.onset;
    correct_norewardC = data.epocs.cNoRewC.onset;
    incorrect_rewardedC = data.epocs.iRewC.onset;
    incorrect_norewardC = data.epocs.iNoRewC.onset;
    [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cueTSC,correct_rewardedC,...
            correct_norewardC,incorrect_rewardedC,incorrect_norewardC);
    
    for e = 1:height(session_identifiers(:,1))
        e1 = session_identifiers(e,1)-baseWindow;
        e2 = e1+timeWindow+baseWindow;
        e3 = session_identifiers(e,1)-baseline(:,2);
        e4 = session_identifiers(e,1)-baseline(:,1);
        [~,ind1] = min(abs(time-e1));
        [~,ind2] = min(abs(time-e2));
        [~,ind3] = min(abs(time-e3));
        [~,ind4] = min(abs(time-e4));
        SIGNAL_meanBase = mean(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_std = std(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_window = SIGNAL_raw(1,ind1:ind2);
        SIGNAL_dFF = SIGNAL_window - SIGNAL_meanBase;
        SIGNAL_dFF = 100*(SIGNAL_dFF / SIGNAL_meanBase);
        SIGNAL_time = time(1,ind1:ind2);
        SIGNAL_z = (SIGNAL_window - SIGNAL_meanBase)/SIGNAL_std;
        if length(SIGNAL_z) < epocArrayLen
            continue
        end
        sessionSTREAM(e,:) = SIGNAL_z(1:epocArrayLen);
        sessionAMP(e,:) = max(SIGNAL_z);
        sessionAUC(e,:) = trapz(SIGNAL_time,SIGNAL_z);
    end
    for e = 1:height(cueTSC)
        e1 = cueTSC(e,1)-baseWindow;
        e2 = e1+timeWindow+baseWindow;
        e3 = cueTSC(e,1)-baseline(:,2);
        e4 = cueTSC(e,1)-baseline(:,1);
        [c1,ind1] = min(abs(time-e1));
        [c2,ind2] = min(abs(time-e2));
        [c3,ind3] = min(abs(time-e3));
        [c4,ind4] = min(abs(time-e4));
        SIGNAL_meanBase = mean(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_std = std(SIGNAL_raw(1,ind4:ind3));
        SIGNAL_window = SIGNAL_raw(1,ind1:ind2);
        SIGNAL_dFF = SIGNAL_window - SIGNAL_meanBase;
        SIGNAL_dFF = 100*(SIGNAL_dFF / SIGNAL_meanBase);
        SIGNAL_time = time(1,ind1:ind2);
        SIGNAL_z = (SIGNAL_window - SIGNAL_meanBase)/SIGNAL_std;
        % if length(SIGNAL_z) < epocArrayLen
        %     break
        % end
        cueSTREAM(e,:) = SIGNAL_z(1:epocArrayLen);
        cueAMP(e,:) = max(SIGNAL_z);
        cueAUC(e,:) = trapz(SIGNAL_time,SIGNAL_z);
    end
end
% Generate a colormap for the heatmap
myColorMap = turbo;
% Set a master font size
myFontSize = 14;
session_time = ts1(:,1:epocArrayLen);
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
xlim([session_time(:,1) session_time(:,epocArrayLen)]);
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
x_values = session_identifiers(2:2:end,1); 

% Define the time array
timeArray = session_time;

% Plot the heatmap
imagesc(timeArray, 1:size(signals, 1), cueSTREAM);
set(gca,'YDir','normal','FontSize',myFontSize);
% Set the colormap for the heatmap
colormap(myColorMap);

% Add a colorbar to the heatmap
colorbar('FontSize',myFontSize);
clim([zMin, zMax]);
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