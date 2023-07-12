clear all; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 10; % the number of seconds after the onset of a TTL to analyze
baseWindow = 2; % the number of seconds before the onset of a TTL to use for stream extraction
baseline = [5 1]; % baseline signal for dFF/zscore
amp_window = [-1 2];
N = 10; %Downsample N times
minArrayLen = 1221; 
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/DA_PREF/AM251_Collab','Choose the .mat files you want to analyze.'); %gets directory%
if myDir == 0
    disp("Select a .mat file to start")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
IDs = {};
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    tempID = cellstr(name);
    IDs = vertcat(IDs,tempID);
    brokenID = strsplit(name,'_');
    animalID = char(brokenID(1));
    treatment = char(brokenID(2));
    wd_phase = char(brokenID(3));
    load(filename)
    
    if isfield(data.streams, 'x405A')
        ISOS = 'x405A';
        SIGNAL = 'x465A';
    elseif isfield(data.streams,'x405C')
        ISOS = 'x405C';
        SIGNAL = 'x465C';
    end
        %time array used for all streams%
        session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
        %removes the first (t) seconds where the data is wild due to turning on LEDs%
        t = 5; % time threshold below which we will discard
        ind = find(session_time>t,1);% find first index of when time crosses threshold
        session_time = session_time(ind:end); % reformat vector to only include allowed time
        data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
        data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
       
        %downsample streams and time array by N times%
        data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
        data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
        minLength = min(length(data.streams.(ISOS).data),length(data.streams.(SIGNAL).data));
        data.streams.(ISOS).data = data.streams.(ISOS).data(1:minLength);
        data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(1:minLength);
        session_time = downsample(session_time, N);
        ts1 = -baseWindow + (1:minArrayLen) / data.streams.(SIGNAL).fs*N;
        ISOS_raw = data.streams.(ISOS).data;
        SIGNAL_raw = data.streams.(SIGNAL).data;
        
        %detrend & dFF%
        bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,1);
        Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
        Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
        dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(dFF));
        detrend_465 = detrend(dFF);
        z465 = zscore(detrend_465);

        leftApp = data.epocs.lApp.onset;
        rightApp = data.epocs.rApp.onset;
        if isempty(leftApp) && isempty(rightApp)
            continue
        end

        for e = 1:height(leftApp)
            if isempty(leftApp)
                leftCupSTREAMdff = zeros(1:minArrayLen);
                leftCupAMPdff = zeros(1);
                leftCupSTREAMz = zeros(1:minArrayLen);
                leftCupAMPz = zeros(1);
                break
            end
            e1 = leftApp(e,1)-baseWindow;
            e2 = e1+timeWindow+baseWindow;
            e3 = leftApp(e,1)-baseline(2);
            e4 = leftApp(e,1)-baseline(1);
            [~,ind1] = min(abs(session_time-e1));
            [~,ind2] = min(abs(session_time-e2));
            [~,ind3] = min(abs(session_time-e3));
            [~,ind4] = min(abs(session_time-e4));
            SIGNAL_dffBase = mean(SIGNAL_raw(1,ind4:ind3));
            SIGNAL_dff = 100*(SIGNAL_raw(1,ind1:ind2) - SIGNAL_dffBase)/SIGNAL_dffBase;
            SIGNAL_zBase = SIGNAL_raw(1,ind4:ind3);
            SIGNAL_z = SIGNAL_raw(1,ind1:ind2);
            epoc_time = session_time(1,ind1:ind2);
            zb = mean(SIGNAL_zBase);
            zsd = std(SIGNAL_dffBase);
            zfinal = (SIGNAL_z - SIGNAL_dffBase)/zsd;
            if length(zfinal) < minArrayLen || length(SIGNAL_dff) < minArrayLen
                continue
            end
            leftCupSTREAMdff(e,:) = SIGNAL_dff(1:minArrayLen);
            leftCupAMPdff(e,:) = max(SIGNAL_dff(1:minArrayLen));
            leftCupSTREAMz(e,:) = zfinal(1:minArrayLen);
            leftCupAMPz(e,:) = max(zfinal(1:minArrayLen));
        end
        for e = 1:height(rightApp)
            if isempty(rightApp)
                rightCupSTREAMdff = zeros(1:minArrayLen);
                rightCupAMPdff = zeros(1);
                rightCupSTREAMz = zeros(1:minArrayLen);
                rightCupAMPz = zeros(1);
                break
            end
            e1 = rightApp(e,1)-baseWindow;
            e2 = e1+timeWindow+baseWindow;
            e3 = rightApp(e,1)-baseline(2);
            e4 = rightApp(e,1)-baseline(1);
            [~,ind1] = min(abs(session_time-e1));
            [~,ind2] = min(abs(session_time-e2));
            [~,ind3] = min(abs(session_time-e3));
            [~,ind4] = min(abs(session_time-e4));
            SIGNAL_dffBase = mean(SIGNAL_raw(1,ind4:ind3));
            SIGNAL_dff = 100*(SIGNAL_raw(1,ind1:ind2) - SIGNAL_dffBase)/SIGNAL_dffBase;
            SIGNAL_zBase = SIGNAL_raw(1,ind4:ind3);
            SIGNAL_z = SIGNAL_raw(1,ind1:ind2);
            epoc_time = session_time(1,ind1:ind2);
            zb = mean(SIGNAL_zBase);
            zsd = std(SIGNAL_dffBase);
            zfinal = (SIGNAL_z - SIGNAL_dffBase)/zsd;
            if length(zfinal) < minArrayLen || length(SIGNAL_dff) < minArrayLen
                continue
            end
            rightCupSTREAMdff(e,:) = SIGNAL_dff(1:minArrayLen);
            rightCupAMPdff(e,:) = max(SIGNAL_dff(1:minArrayLen));
            rightCupSTREAMz(e,:) = zfinal(1:minArrayLen);
            rightCupAMPz(e,:) = max(zfinal(1:minArrayLen));
        end
        

    
        
    
    master_left_STREAMdff(i,:) = mean(leftCupSTREAMdff,1);
    master_right_STREAMdff(i,:) = mean(rightCupSTREAMdff,1);
    master_left_STREAMz(i,:) = mean(leftCupSTREAMz,1);
    master_right_STREAMz(i,:) = mean(rightCupSTREAMz,1);
    
    leftCupAMPdff(leftCupAMPdff == 0) = nan;
    rightCupAMPdff(rightCupAMPdff == 0) = nan;
    AMPdff_analysis(i,1:2) = {mean(leftCupAMPdff,1) mean(rightCupAMPdff,1)};
    leftCupAMPz(leftCupAMPz == 0) = nan;
    rightCupAMPz(rightCupAMPz == 0) = nan;
    AMPz_analysis(i,1:2) = {mean(leftCupAMPz,1) mean(rightCupAMPz,1)};

    
        
end
% IDs = cellstr(IDs);
IDlist = cell2table(IDs,'VariableNames',{'ID'});
AMPdff_analysis(cellfun(@(x) x==0,AMPdff_analysis)) = {NaN};
AMPdff_analysis_table = cell2table(AMPdff_analysis,'VariableNames',{'Left Approach','Right Approach'});
AMPz_analysis(cellfun(@(x) x==0,AMPz_analysis)) = {NaN};
AMPz_analysis_table = cell2table(AMPz_analysis,'VariableNames',{'Left Approach','Right Approach'});
AMPdff_analysis_table = horzcat(IDlist,AMPdff_analysis_table);
toc

disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);