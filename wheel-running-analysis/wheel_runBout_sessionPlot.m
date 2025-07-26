clear; close all

channel = 2;
camNumber = 2;
smoothFactor = 50;
fit = 1; % fit for polynomial
dataType = 1; % 1 = mat, 2 = tank
manualTTL = 2; % 1 = yes, 0 = no, 2 = TTLs in mat file
N = 51; %Downsample N times
BLOCKPATH = '/Users/brandon/personal-drive/hot-wheels/data/wheel-running-mats/unlocked/DA40_UL_Day7_10_3_2022.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dataType == 1
    load(BLOCKPATH)
elseif dataType == 2
    data = TDTbin2mat(BLOCKPATH);
end
if manualTTL == 0
    if camNumber == 1
        ind = double(data.epocs.Cam1.notes.index);
        ts = data.epocs.Cam1.notes.ts;
    elseif camNumber == 2
        ind = double(data.epocs.Cam2.notes.index);
        ts = data.epocs.Cam2.notes.ts;
    else
        disp('Invalid camera number')
    end
    var1 = [ind ts];
    %separate by index to make separate epocs%
    runStart = var1(ismember(var1(:,1),[1]),:);
    runStop = var1(ismember(var1(:,1),[2]),:);
    %extract time stamps%
    runStartTs = runStart(:,2);
    runStopTs = runStop(:,2);
    %extract indicies%
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
    if camNumber == 1
        data.epocs.runBout.typeStr = data.epocs.Cam1.typeStr;
    elseif camNumber == 2
        data.epocs.runBout.typeStr = data.epocs.Cam2.typeStr;
    end
    data.epocs.runBout.data = ones(length(data.epocs.runBout.onset),1)*5;    
    %Finds the running bout time stamps%
    runBout_ts = [data.epocs.runBout.onset data.epocs.runBout.offset]; 
    end_ts = height(runBout_ts);
elseif manualTTL == 1
    %Finds the running bout time stamps%
    runBout_ts = [data.epocs.runBout.onset data.epocs.runBout.offset]; 
    end_ts = height(runBout_ts);
elseif manualTTL == 2
    runBout_ts = [data.epocs.runStart.onset data.epocs.runStop.onset];
    end_ts = height(runBout_ts);
end 

%Stream Stores%
if channel == 1
    ISOS = 'x405A'; % name of the 405A store
    SIGNAL = 'x465A'; % name of the 465A store
elseif channel == 2
    ISOS = 'x405C';
    SIGNAL = 'x465C';
else
    error('Unknown channel number')
end


%time array used for all streams%
time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 5; % time threshold below which we will discard
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
min_time = min(time);

%downsample streams and time array by N times%
data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
time = downsample(time, N);

%detrend & dFF%
bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,fit);
Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
detrend_signal = detrend(dFF);
detrend_signal = smoothdata(detrend_signal,'movmean',smoothFactor);

% Plot signal
figure;
plot(time, detrend_signal, 'b'); hold on;
xlabel('Time (s)');
ylabel('zScore');
title('eCB Signal During Wheel Running');
ylim padded;
y_limits = ylim;
margin = 0.3 * diff(y_limits);
% Add markers for each onset-offset pair
for i = 1:length(runBout_ts)
    if runBout_ts(i,2) - runBout_ts(i,1) >= 20
        % Shade region
        fill([runBout_ts(i,1) runBout_ts(i,2) runBout_ts(i,2) runBout_ts(i,1)],...
            [y_limits(1)-margin y_limits(1)-margin y_limits(2)+margin...
            y_limits(2)+margin], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
end

legend('eCB Signal', 'Running Bout');