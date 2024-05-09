clear;
close all;
clc;
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dType = 1; % 1 = TDT, 2 = mat
t = 10; % first t seconds are discarded to remove LED on artifact
N = 100; % downsample signal N times
channel = 1; % 1 = A, 2 = C
ARTIFACT465 = 100; % eliminates parts of the signal with values greater than this threshold
saveArtifact = 0; % 0 = do not save, 1 = save, 2 = overwrite; *only works for mat files*
figSave = 1; % saves figure to folder of data file
session_duration = 5400; % duration of session in seconds, used for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thresh_base1 and thresh_base2 are used as the windows to find peaks in. 
% If useOtherThresh is set to 1, then the peak threshold will be calculated 
% from the other_thresh1 and other_thresh2 windows but the peaks will be searched for in 
% thresh_base1 and thresh_base2. This is useful for when you want to use the signal 
% from the habituation period to calculate the peak threshold for the amphetamine period.
% If useAmphEpoc is set to 1, then the amphetamine time stamp in the file will be used 
% as the starting point of the amphetamine period rather than 1800 seconds.
useAmphEpoc = 1; % 0 = use time, 1 = use epoc
thresh_base1 = [0, 1800]; % habituation window
thresh_base2 = [1800, 5400]; % treatment window
useOtherThresh = 1; % 0 = use thresh_base, 1 = use other_thresh
other_thresh1 = [1200, 1800]; % other habituation window
other_thresh2 = [1200, 1800]; % other amphetamine window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dType == 1
    % Gets tank from UI pop-up window
    TANK_NAME = uigetdir('H:\\My Drive\\wheel-peak-test\\tanks', 'Select a tank to plot');
    if TANK_NAME == 0
        disp('Select a tank to start!')
        return
    end
    disp('Data Type: TDT Tank.')
    data = TDTbin2mat(TANK_NAME, 'T2', session_duration+t, 'TYPE', {'streams','epocs'});
elseif dType == 2
    [TANK_NAME, pathname] = uigetfile('E:\\Google Drive\\hot-wheels\\within-session-wheel-mats');
    
    if TANK_NAME == 0
        disp('Select a mat to start!')
        return
    end
    [~,name,~] = fileparts(TANK_NAME);
    TANK_NAME = strcat(pathname,TANK_NAME);
    load(TANK_NAME)
else
    disp('Select a data type')
end
if channel == 1
    ISOS = 'x405A';
    Grab = 'x465A';
elseif channel == 2
    ISOS = 'x405C';
    Grab = 'x465C';
else
    disp('Cannot resolve channel.')
    return
end

if strcmp(Grab,'x465A')
    ROI = 'DLS';
    thresh_mult = 0.5; % finds peaks that are 50% of the max peak
    prom = 3; % prominence threshold for peak detection
    if saveArtifact == 1 && ~isfield(data, "DLSartifact")
        data.DLSartifact = ARTIFACT465;
        disp('New artifact level saved')
    elseif saveArtifact == 2
        data.DLSartifact = ARTIFACT465;
        disp('Overwriting artifact level')
    elseif saveArtifact == 0 && isfield(data, "DLSartifact")
        ARTIFACT465 = data.DLSartifact;
        disp('Artifact level loaded')
    else
        disp('Artifact level not saved')
    end
elseif strcmp(Grab,'x465C') 
    ROI = 'NAc';
    thresh_mult = 0.3; % finds peaks that are 30% of the max peak
    prom = 3; % prominence threshold for peak detection
    if saveArtifact == 1 && ~isfield(data, "NACartifact")
        data.NACartifact = ARTIFACT465;
        disp('New artifact level saved')
    elseif saveArtifact == 2
        data.NACartifact = ARTIFACT465;
        disp('Overwriting artifact level')
    elseif saveArtifact == 0 && isfield(data, "NACartifact")
        ARTIFACT465 = data.NACartifact;
        disp('Artifact level loaded')
    else
        disp('Artifact level not saved')
    end
else
    disp('Cannot resolve ROI.')
    return
end

if saveArtifact == 1 || saveArtifact == 2
    disp('Saved.')
    save(TANK_NAME, 'data')
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

% Display max dFF
fprintf('Max dFF: %.2f\n', max(Grab_dFF))

% noise reduction using moving median
Grab_filt = smoothdata(Grab_dFF,'movmedian',10);


% artifact removal
% Find indices where the signal is greater than the artifact threshold
art1_indices = find(any(Grab_filt > ARTIFACT465, 1));
% Find indices where the signal is less than the negative artifact threshold
art2_indices = find(any(Grab_filt < -ARTIFACT465, 1));
% Combine the indices of artifacts
artifact_indices = unique([art1_indices, art2_indices]);
% Create a logical array indicating good (non-artifact) signals
good_signals = true(1, size(Grab_filt, 2));
good_signals(artifact_indices) = false;
% Filter out the artifact signals
Grab_filt = Grab_filt(:, good_signals);
time_filt = time(:, good_signals);

f1 = figure(1);
% thresh_base1
subplot(4,1,1)
plot(time, Grab_dFF)
title('Without artifact removal')
ylabel('dFF')
xlim([t, session_duration])

subplot(4,1,2)
plot(time_filt, Grab_filt)
title('With artifact removal')
ylabel('dFF')
xlim([t, session_duration])

idx1 = find(thresh_base1(1) < time_filt & time_filt < thresh_base1(2));
if useOtherThresh == 0
    threshold1 = max(Grab_filt(:,idx1))*thresh_mult;
elseif useOtherThresh == 1
    alt_idx1 = find(other_thresh1(1) < time_filt & time_filt < other_thresh1(2));
    threshold1 = max(Grab_filt(:,alt_idx1))*thresh_mult;
end
[pks,locs,w,p] = findpeaks(Grab_filt(:,idx1), 'MinPeakHeight',threshold1, 'MinPeakProminence',prom);
numPeaks = length(locs);
peakFq = num2str(((numPeaks/(thresh_base1(2) - thresh_base1(1) - t))*60), '%.4f');

disp('%%%%%%%%%%%%%% HABITUATION WINDOW %%%%%%%%%%%%%%')
fprintf('Average Peak Width (HAB): %.2f\n',mean(w))
fprintf('Average Peak Prominence (HAB): %.2f\n',mean(p))
fprintf('AUC (HAB): %.2f\n',trapz(time_filt(:,idx1), Grab_filt(:,idx1)))

subplot(4,1,3)
plot(time_filt(:, idx1), Grab_filt(:, idx1))
title(sprintf('ROI: %s, Window duration: %d, Prom: %d, ThreshMax: %.2f, Peaks: %d, Pk/m: %s', ROI, (thresh_base1(2) - thresh_base1(1)), prom, thresh_mult, numPeaks, peakFq))
ylabel('dFF')


xlim([thresh_base1(1), thresh_base1(2)]);

hold on
plot(time_filt(idx1(1) + locs - 1), pks, 'ro'); % Adjust locs by adding the starting index of the chunk
hold off


% thresh_base2

if useOtherThresh == 0 && useAmphEpoc == 1
    idx2 = find(data.epocs.Coc_.onset < time_filt & time_filt < thresh_base2(2));
    threshold2 = max(Grab_filt(:,idx2))*thresh_mult;
elseif useOtherThresh == 0 && useAmphEpoc == 0
    idx2 = find(thresh_base2(1) < time_filt & time_filt < thresh_base2(2));
    threshold2 = max(Grab_filt(:,idx2))*thresh_mult;
elseif useOtherThresh == 1 && useAmphEpoc == 0
    idx2 = find(thresh_base2(1) < time_filt & time_filt < thresh_base2(2));
    other_idx2 = find(other_thresh2(1) < time_filt & time_filt < other_thresh2(2));
    threshold2 = max(Grab_filt(:,other_idx2))*thresh_mult;
elseif useOtherThresh == 1 && useAmphEpoc == 1
    idx2 = find(data.epocs.Coc_.onset < time_filt & time_filt < thresh_base2(2));
    other_idx2 = find(other_thresh2(1) < time_filt & time_filt < other_thresh2(2));
    threshold2 = max(Grab_filt(:,other_idx2))*thresh_mult;
else
    disp('Cannot resolve threshold.')
    return
end
[pks2,locs2,w2,p2] = findpeaks(Grab_filt(:,idx2), 'MinPeakHeight',threshold2, 'MinPeakProminence',prom);
numPeaks2 = length(locs2);
if useAmphEpoc == 1
    peakFq2 = num2str((numPeaks2/((thresh_base2(2) - data.epocs.Coc_.onset))*60), '%.4f');
else
    peakFq2 = num2str((numPeaks2/((thresh_base2(2) - thresh_base2(1)))*60), '%.4f');
end


% creates a new signal that is zeroed to the amphetamine onset
auc_idx = find(time_filt > data.epocs.Coc_.onset,1);
if Grab_filt(:,auc_idx) < 0
    diff = 0 - Grab_filt(:,auc_idx:end); 
    Grab_auc = Grab_filt(:,auc_idx:end) + abs(diff);
elseif Grab_filt(:,auc_idx1) > 0
    diff = 0 - Grab_filt(:,auc_idx:end); 
    Grab_auc = Grab_filt(:,auc_idx:end) - abs(diff);
end
disp('%%%%%%%%%%%%%% AMPHETAMINE WINDOW %%%%%%%%%%%%%%')
fprintf('Average Peak Width (AMP): %.2f\n',mean(w2))
fprintf('Average Peak Prominence (AMP): %.2f\n',mean(p2))
fprintf('AUC (AMP): %.2f\n',trapz(time_filt(:,auc_idx:end), Grab_auc))

subplot(4,1,4)
plot(time_filt(:,idx2), Grab_filt(:,idx2))
title(sprintf('ROI: %s, Window duration: %d, Prom: %d, ThreshMax: %.2f, Peaks: %d, Pk/m: %s', ROI, (thresh_base2(2) - thresh_base2(1)), prom, thresh_mult, numPeaks2, peakFq2))
ylabel('dFF')


xlim([thresh_base2(1), thresh_base2(2)]);

hold on
plot(time_filt(idx2(1) + locs2 - 1), pks2, 'ro'); % Adjust locs by adding the starting index of the chunk
hold off



% save figure logic
if figSave == 1 && dType == 2
    disp('Saving figure.')
    % Save figure 1 to pathname
    file_name = char(strcat(pathname,name,{'_'},ROI,'.pdf'));
    orient(f1,'landscape');
    print(f1,file_name,'-dpdf','-vector','-bestfit','');
elseif figSave == 1 && dType == 1
    disp('Saving figure.')
    % Save figure 1 to pathname
    file_name = char(strcat(TANK_NAME,{'_'},ROI,'.pdf'));
    orient(f1,'landscape');
    print(f1,file_name,'-dpdf','-vector','-bestfit','');
elseif figSave == 0
    disp('Figure not saved')
end

% create data table
data_table = [mean(w), mean(w2), mean(p), mean(p2), trapz(time_filt(:,idx1), Grab_filt(:,idx1)), trapz(time_filt(:,auc_idx:end), Grab_auc)];
data_table = array2table(data_table, 'VariableNames', {'AvgPeakWidth_HAB', 'AvgPeakWidth_AMP', 'AvgPeakProm_HAB', 'AvgPeakProm_AMP', 'AUC_HAB', 'AUC_AMP'});
