clear all; clc; close all;
warning off;


savetype = ".jpg";
savepath = '/Volumes/CUDADRIVE/DA_Wheel/Wheel_Running_Figs/Wheel_Whole_Session_Plots/';
myDir = uigetdir; %gets directory%
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv','.csv'}));
fig = cell(length(myFiles),1);
FILENAME = cell(length(myFiles),1);
FIG_TITLE = cell(length(myFiles),1);
disp("Starting batch plot...")
for batch = 1:length(myFiles)
    fprintf('Loading tank %d of %d...\n',batch,length(myFiles))
    BLOCKPATH = fullfile(myDir, myFiles(batch).name);
    [~,name,~] = fileparts(BLOCKPATH);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    animalID = char(brokenID{1});
    wheelAccess = char(brokenID{2});
    day = char(brokenID{3});
    
   
    DLS_ISOS = 'x405A'; % name of the 405 store
    DLS_DA = 'x465A'; % name of the 465 store
    NAc_ISOS = 'x405C'; % name of the 405 store
    NAc_DA = 'x465C'; % name of the 465 store
    FILENAME{batch} = strcat(animalID,'_', wheelAccess,'_', day);
    FIG_TITLE{batch} = strcat(animalID," ",wheelAccess," ",day);
    
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs','streams'});
    N = 1000; %Downsample N times
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

    disp("Plotting signals...")
    fig{batch} = figure;
    subplot(2,1,1)
    plot(time,detrend_465A,'r');
    ylabel("Z-Score")
    title("DLS 465A")
    sgtitle(FIG_TITLE{batch});
    subplot(2,1,2)
    plot(time,detrend_465C,'b');
    title("NAc 465C")
    xlabel("Time (s)")
    ylabel("Z-Score")
    
    
    

    savdisp = strcat("Saving figure for: ",FILENAME{batch});
    disp(savdisp)
    filename = strcat(savepath,FILENAME{batch},savetype);
    saveas(fig{batch},filename)
    clf reset   
end
close all
fprintf("Finished plotting and saving %d figures...\n",length(myFiles))
