clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore
% amp_window = [-1 2];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths
%as of this version, the easiest way to determine minArrayLength is to set
%it low, run the script, and see how large the arrays are coming out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/DA_PRL','Choose the .mat files you want to analyze.'); %gets directory%
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
treatList = {};

AMPdFF_analysis = cell(numFiles,5);
AMPz_analysis = cell(numFiles,5);
AUCdFF_analysis = cell(numFiles,5);
AUCz_analysis = cell(numFiles,5);
master_cue_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cNoRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_iRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_iNoRew_STREAMdFF = zeros(numFiles,epocArrayLen);
master_cue_STREAMz = zeros(numFiles,epocArrayLen);
master_cRew_STREAMz = zeros(numFiles,epocArrayLen);
master_cNoRew_STREAMz = zeros(numFiles,epocArrayLen);
master_iRew_STREAMz = zeros(numFiles,epocArrayLen);
master_iNoRew_STREAMz = zeros(numFiles,epocArrayLen);

for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    [~,treatment,~] = fileparts(myDir);
    tempID = cellstr(name);
    tempTreat = cellstr(treatment);
    IDs = vertcat(IDs,tempID);
    treatList = vertcat(treatList,tempTreat);
    brokenID = strsplit(name,'_');
    prl_phase = char(brokenID(2));
    load(filename)

    if isfield(data.streams, 'x405A')
        ISOS = 'x405A';
        SIGNAL = 'x465A';

        cue = data.epocs.St1_.onset;
        cRew = data.epocs.cRewA.onset;
        cNoRew = data.epocs.cNoRewA.onset;
        iRew = data.epocs.iRewA.onset;
        iNoRew = data.epocs.iNoRewA.onset;

        cueSTREAMraw = zeros(height(cue),epocArrayLen);
        cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);

        cueSTREAMdFF = zeros(height(cue),epocArrayLen);
        cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
        
        cueSTREAMz = zeros(height(cue),epocArrayLen);
        cRewSTREAMz = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMz = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);

        cueAMPdFF = zeros(height(cue),1);
        cRewAMPdFF = zeros(height(cRew),1);
        cNoRewAMPdFF = zeros(height(cNoRew),1);
        iRewAMPdFF = zeros(height(iRew),1);
        iNoRewAMPdFF = zeros(height(iNoRew),1);

        cueAMPz = zeros(height(cue),1);
        cRewAMPz = zeros(height(cRew),1);
        cNoRewAMPz = zeros(height(cNoRew),1);
        iRewAMPz = zeros(height(iRew),1);
        iNoRewAMPz = zeros(height(iNoRew),1);

        cueAUCdFF = zeros(height(cue),1);
        cRewAUCdFF = zeros(height(cRew),1);
        cNoRewAUCdFF = zeros(height(cNoRew),1);
        iRewAUCdFF = zeros(height(iRew),1);
        iNoRewAUCdFF = zeros(height(iNoRew),1);
        
        cueAUCz = zeros(height(cue),1);
        cRewAUCz = zeros(height(cRew),1);
        cNoRewAUCz = zeros(height(cNoRew),1);
        iRewAUCz = zeros(height(iRew),1);
        iNoRewAUCz = zeros(height(iNoRew),1);



        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
        epocList = {cue;cRew;cNoRew;iRew;iNoRew};
        outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
            iRewSTREAMraw;iNoRewSTREAMraw};
        outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
            iRewSTREAMdFF;iNoRewSTREAMdFF};
        outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
            iRewSTREAMz;iNoRewSTREAMz};
        outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF};
        outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz};
        outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF};
        outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz};
         
        % data.analysis.sessionID = session_identifiers;
        
    elseif isfield(data.streams, 'x405C')
        ISOS = 'x405C';
        SIGNAL = 'x465C';

        cue = data.epocs.St2_.onset;
        cRew = data.epocs.cRewC.onset;
        cNoRew = data.epocs.cNoRewC.onset;
        iRew = data.epocs.iRewC.onset;
        iNoRew = data.epocs.iNoRewC.onset;

        cueSTREAMraw = zeros(height(cue),epocArrayLen);
        cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);

        cueSTREAMdFF = zeros(height(cue),epocArrayLen);
        cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
        
        cueSTREAMz = zeros(height(cue),epocArrayLen);
        cRewSTREAMz = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMz = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);

        cueAMPdFF = zeros(height(cue),1);
        cRewAMPdFF = zeros(height(cRew),1);
        cNoRewAMPdFF = zeros(height(cNoRew),1);
        iRewAMPdFF = zeros(height(iRew),1);
        iNoRewAMPdFF = zeros(height(iNoRew),1);
        
        cueAMPz = zeros(height(cue),1);
        cRewAMPz = zeros(height(cRew),1);
        cNoRewAMPz = zeros(height(cNoRew),1);
        iRewAMPz = zeros(height(iRew),1);
        iNoRewAMPz = zeros(height(iNoRew),1);

        cueAUCdFF = zeros(height(cue),1);
        cRewAUCdFF = zeros(height(cRew),1);
        cNoRewAUCdFF = zeros(height(cNoRew),1);
        iRewAUCdFF = zeros(height(iRew),1);
        iNoRewAUCdFF = zeros(height(iNoRew),1);
        
        cueAUCz = zeros(height(cue),1);
        cRewAUCz = zeros(height(cRew),1);
        cNoRewAUCz = zeros(height(cNoRew),1);
        iRewAUCz = zeros(height(iRew),1);
        iNoRewAUCz = zeros(height(iNoRew),1);

        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);
        epocList = {cue;cRew;cNoRew;iRew;iNoRew};
        outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
            iRewSTREAMraw;iNoRewSTREAMraw};
        outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
            iRewSTREAMdFF;iNoRewSTREAMdFF};
        outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
            iRewSTREAMz;iNoRewSTREAMz};
        outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF};
        outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz};
        outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF};
        outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz};
        
        % data.analysis.sessionID = session_identifiers;
    end
    %time array used for all streams%
    session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    ind = find(session_time>t,1);% find first index of when time crosses threshold
    session_time = session_time(ind:end); % reformat vector to only include allowed time
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
    data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
   
    %downsample streams and time array by N times%
    data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
    data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
    minStreamLength = min(length(data.streams.(ISOS).data),length(data.streams.(SIGNAL).data));
    data.streams.(ISOS).data = data.streams.(ISOS).data(1:minStreamLength);
    data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(1:minStreamLength);

    ISOS_raw = data.streams.(ISOS).data;
    SIGNAL_raw = data.streams.(SIGNAL).data;
    
    % % Perform ratiometric correction
    % ratio = SIGNAL_raw ./ ISOS_raw; % calculate ratio of raw signal to isosbestic signal
    % ratio_smoothed = smoothdata(ratio); % smooth the ratio signal to remove noise
    % SIGNAL_filtered = ratio_smoothed .* ISOS_raw; % multiply smoothed ratio signal by isosbestic signal
    
    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;

    
    %detrend & dFF%
    bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,1);
    Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
    Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    detrend_dFF = detrend(dFF);
    
    cueArray = session_identifiers(1:2:end,:);
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(1:end-1,:);
        cueArray = cueArray(1:end-1,1);
    end
    leverArray = session_identifiers(2:2:end,:);
    levers_z = zeros(height(leverArray),epocArrayLen);
    lever_raw = zeros(height(leverArray),epocArrayLen);
    for m = 1:height(leverArray)
        
        cueBase1 = cueArray(m,1) - 3;
        cueBase2 = cueArray(m,1) - 1;
        [~,cueSt] = min(abs(session_time - cueBase1));
        [~,cueEn] = min(abs(session_time - cueBase2));
        cueBaseMean(m,1) = mean(SIGNAL_raw(1,cueSt:cueEn));
        cueBaseStd(m,1) = std(SIGNAL_raw(1,cueSt:cueEn));

        leverStart = leverArray(m,1) - baseWindow;
        leverEnd = leverStart + timeWindow + baseWindow;
        [~,levSt] = min(abs(session_time - leverStart));
        [~,levEn] = min(abs(session_time - leverEnd));
        leverSigRaw = SIGNAL_raw(1,levSt:levEn);
        if length(leverSigRaw) < epocArrayLen
            mn = mean(leverSigRaw(1,end-10:end));
            leverSigRaw(1,end:epocArrayLen) = mn;
        elseif length(leverSigRaw) > epocArrayLen
            op = length(leverSigRaw);
            arrayDif = op - epocArrayLen;
            leverSigRaw = leverSigRaw(1,1:end-arrayDif);
        end
            
        lever_raw(m,:) = leverSigRaw;
                
    end
    sessionSTREAMSraw{i,1} = lever_raw;
    [~,baseSt] = min(abs(ts1 - baseline(1)));
    [~,baseEn] = min(abs(ts1 - baseline(2)));
    epocON = find(ts1>=0,1);
    epocOFF = find(ts1<=timeWindow,1,'last');
    for n = 1:height(lever_raw)
        % dF/F
        meanCue = cueBaseMean(n,1);
        stdCue = cueBaseStd(n,1);
        
        levers_dFF(n,1:epocArrayLen) = lever_raw(n, 1:epocArrayLen) - meanCue;
        levers_dFF(n,1:epocArrayLen) = 100*(levers_dFF(n,1:epocArrayLen) / meanCue);
        % z-Score
        meanCueBase_dFF = mean(levers_dFF(n,baseSt:baseEn));
        stdCueBase_dFF = std(levers_dFF(n,baseSt:baseEn));
        levers_z(n,1:epocArrayLen) = (levers_dFF(n,1:epocArrayLen) - meanCueBase_dFF) / stdCueBase_dFF;
    end
    
    sessionSTREAMSz{i,1} = mean(levers_z);
    % trialNames = cellstr(num2str(session_identifiers(2:2:end,2)));
    trialNames = array2table(session_identifiers(2:2:end,2),'VariableNames',{'Trial_Type'});
    % epocTime = cellstr(num2str(ts1'));
    % trialNames = cell2table(trialNames,'VariableNames',{'Trial Type'});
    levers_z = array2table(levers_z);
    levers_z = horzcat(trialNames,levers_z);
    cRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 1, 2:end));
    cNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 2, 2:end));
    iRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 3, 2:end));
    iNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 4, 2:end));
    
    leverLatencies(i,1) = mean(leverArray(:,1) - cueArray(:,1));
    % data.analysis.lever_z = levers_z;
    for k = 1:numel(epocList)
            epoc = double(epocList{k});
            streams_raw = double(outputSTREAMSraw{k,1});
            streams_dFF = double(outputSTREAMSdFF{k,1});
            streams_z = double(outputSTREAMSz{k,1});
            ampdFF = double(outputAMPdFF{k,1});
            ampZ = double(outputAMPz{k,1});
            aucdFF = double(outputAUCdFF{k,1});
            aucZ = double(outputAUCz{k,1});
        for ii = 1:height(epoc)
            if epoc(ii) == 0
                streams_dFF(ii,1:epocArrayLen) = NaN;
                ampdFF(ii) = zeros(1);
                streams_z(ii,1:epocArrayLen) = NaN;
                ampZ(ii) = zeros(1);
                aucdFF(ii) = zeros(1);
                aucZ(ii) = zeros(1);
                continue
            end
            
            windowStart = epoc(ii)-baseWindow;
            windowEnd = windowStart+timeWindow+baseWindow;
            [~,windSt] = min(abs(session_time - windowStart));
            [~,windEn] = min(abs(session_time - windowEnd));
            epocSigRaw = SIGNAL_raw(1,windSt:windEn);


            if length(epocSigRaw) < epocArrayLen
                mn = mean(epocSigRaw(1,end-10:end));
                epocSigRaw(1,end:epocArrayLen) = mn;
            elseif length(epocSigRaw) > epocArrayLen
                op = length(epocSigRaw);
                arrayDif = op - epocArrayLen;
                epocSigRaw = epocSigRaw(1,1:end-arrayDif);
            end
            streams_raw(ii,1:epocArrayLen) = epocSigRaw;
        end

        [~,baseSt] = min(abs(ts1 - baseline(1)));
        [~,baseEn] = min(abs(ts1 - baseline(2)));
        epocON = find(ts1>=0,1);
        epocOFF = find(ts1<=timeWindow,1,'last');
        % streams_dFF = zeros(size(streams_raw));
        % streams_z = zeros(size(streams_raw));
        
        for j = 1:height(streams_raw)
            % dF/F
            meanBase = mean(streams_raw(j,baseSt:baseEn));
            stdBase = std(streams_raw(j,baseSt:baseEn));
            streams_dFF(j,1:epocArrayLen) = streams_raw(j,1:epocArrayLen) - meanBase;
            streams_dFF(j,1:epocArrayLen) = 100*(streams_dFF(j,1:epocArrayLen) / meanBase);
            % z-Score
            meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
            stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
            streams_z(j,1:epocArrayLen) = (streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
            % amplitude
            ampdFF(j) = max(streams_dFF(j,1:epocArrayLen));
            ampZ(j) = max(streams_z(j,1:epocArrayLen));
            % AUC
            aucdFF(j) = trapz(ts1(1,epocON:epocOFF),streams_dFF(j,epocON:epocOFF));
            aucZ(j) = trapz(ts1(1,epocON:epocOFF),streams_z(j,epocON:epocOFF));
        end
        
            outputSTREAMSraw{k,1} = streams_raw(:,1:epocArrayLen);
            outputSTREAMSdFF{k,1} = streams_dFF(:,1:epocArrayLen);
            outputSTREAMSz{k,1} = streams_z(:,1:epocArrayLen);
            outputAMPdFF{k,1} = ampdFF;
            outputAMPz{k,1} = ampZ;
            outputAUCdFF{k,1} = aucdFF;
            outputAUCz{k,1} = aucZ;
        
    end
    



    master_cue_STREAMraw(i,:) = mean(outputSTREAMSraw{1,1},1);
    master_cRew_STREAMraw(i,:) = mean(outputSTREAMSraw{2,1},1);
    master_cNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{3,1},1);
    master_iRew_STREAMraw(i,:) = mean(outputSTREAMSraw{4,1},1);
    master_iNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{5,1},1);

    master_cue_STREAMdFF(i,:) = mean(outputSTREAMSdFF{1,1},1);
    master_cRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{2,1},1);
    master_cNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{3,1},1);
    master_iRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{4,1},1);
    master_iNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{5,1},1);

    master_cue_STREAMz(i,:) = mean(outputSTREAMSz{1,1},1);
    master_cRew_STREAMz(i,:) = mean(outputSTREAMSz{2,1},1);
    master_cNoRew_STREAMz(i,:) = mean(outputSTREAMSz{3,1},1);
    master_iRew_STREAMz(i,:) = mean(outputSTREAMSz{4,1},1);
    master_iNoRew_STREAMz(i,:) = mean(outputSTREAMSz{5,1},1);

    master_cRew_cueBase_STREAMz(i,:) = mean(cRew_cueBase,1);
    master_cNoRew_cueBase_STREAMz(i,:) = mean(cNoRew_cueBase,1);
    master_iRew_cueBase_STREAMz(i,:) = mean(iRew_cueBase,1);
    master_iNoRew_cueBase_STREAMz(i,:) = mean(iNoRew_cueBase,1);
   
    
    AMPdFF_analysis(i,1:5) = {mean(outputAMPdFF{1},1) mean(outputAMPdFF{2},1)...
        mean(outputAMPdFF{3},1) mean(outputAMPdFF{4},1) mean(outputAMPdFF{5},1)};
    AMPz_analysis(i,1:5) = {mean(outputAMPz{1},1) mean(outputAMPz{2},1)...
        mean(outputAMPz{3},1) mean(outputAMPz{4},1) mean(outputAMPz{5},1)};

    AUCdFF_analysis(i,1:5) = {mean(outputAUCdFF{1},1) mean(outputAUCdFF{2},1)...
        mean(outputAUCdFF{3},1) mean(outputAUCdFF{4},1) mean(outputAUCdFF{5},1)};
    AUCz_analysis(i,1:5) = {mean(outputAUCz{1},1) mean(outputAUCz{2},1)...
        mean(outputAUCz{3},1) mean(outputAUCz{4},1) mean(outputAUCz{5},1)};

    % save(filename,'data');
        
end
IDlist = cell2table(IDs,'VariableNames',{'ID'});
treatList = cell2table(treatList,'VariableNames',{'Treatment'});

leverLatencies = array2table(leverLatencies,'VariableNames',{'Lever Press Latency'});
leverLatencies = horzcat(IDlist,treatList,leverLatencies);

%% Amplitude Table %%
% AMPdFF_analysis(cellfun(@(x) x==0,AMPdFF_analysis,'UniformOutput',false)) = {NaN};
AMPdFF_analysis_table = cell2table(AMPdFF_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AMPdFF_analysis_table = horzcat(IDlist,treatList,AMPdFF_analysis_table);
% AMPz_analysis(cellfun(@(x) x==0,AMPz_analysis,'UniformOutput',false)) = {NaN};
AMPz_analysis_table = cell2table(AMPz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AMPz_analysis_table = horzcat(IDlist,treatList,AMPz_analysis_table);
%% Area Under Curve Table %%
% AUCdFF_analysis(cellfun(@(x) x==0,AUCdFF_analysis,'UniformOutput',false)) = {NaN};
AUCdFF_analysis_table = cell2table(AUCdFF_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AUCdFF_analysis_table = horzcat(IDlist,treatList,AUCdFF_analysis_table);
% AUCz_analysis(cellfun(@(x) x==0,AUCz_analysis,'UniformOutput',false)) = {NaN};
AUCz_analysis_table = cell2table(AUCz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
AUCz_analysis_table = horzcat(IDlist,treatList,AUCz_analysis_table);

%% Plots for each epoc %%
f1 = figure;
subplot(5,1,1)
plot(ts1,mean(master_cue_STREAMz,1,'omitnan'),'b')
hold on
cue_std = std(master_cue_STREAMz,'omitnan')/sqrt(height(master_cue_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cue_STREAMz,1,'omitnan') + cue_std, fliplr(mean(master_cue_STREAMz,1,'omitnan') - cue_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,2)
plot(ts1,mean(master_cRew_STREAMz,1,'omitnan'),'b')
hold on
cRew_std = std(master_cRew_STREAMz,'omitnan')/sqrt(height(master_cRew_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cRew_STREAMz,1,'omitnan') + cRew_std, fliplr(mean(master_cRew_STREAMz,1,'omitnan') - cRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Correct Reward')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,3)
plot(ts1,mean(master_cNoRew_STREAMz,1,'omitnan'),'b')
hold on
cNoRew_std = std(master_cNoRew_STREAMz,'omitnan')/sqrt(height(master_cNoRew_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cNoRew_STREAMz,1,'omitnan') + cNoRew_std, fliplr(mean(master_cNoRew_STREAMz,1,'omitnan') - cNoRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Correct No Reward')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,4)
plot(ts1,mean(master_iRew_STREAMz,1,'omitnan'),'b')
hold on
iRew_std = std(master_iRew_STREAMz,'omitnan')/sqrt(height(master_iRew_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_iRew_STREAMz,1,'omitnan') + iRew_std, fliplr(mean(master_iRew_STREAMz,1,'omitnan') - iRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Incorrect Reward')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,5)
plot(ts1,mean(master_iNoRew_STREAMz,1,'omitnan'),'b')
hold on
iNoRew_std = std(master_iNoRew_STREAMz,'omitnan')/sqrt(height(master_iNoRew_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_iNoRew_STREAMz,1,'omitnan') + iNoRew_std, fliplr(mean(master_iNoRew_STREAMz,1,'omitnan') - iNoRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Incorrect No Reward')
ylabel('GrabDA dF/F (z-Score)')
xlabel('Time (s)')

f2 = figure;
subplot(5,1,1)
plot(ts1,mean(master_cue_STREAMz,1,'omitnan'),'b')
hold on
cue_std = std(master_cue_STREAMz,'omitnan')/sqrt(height(master_cue_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cue_STREAMz,1,'omitnan') + cue_std, fliplr(mean(master_cue_STREAMz,1,'omitnan') - cue_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Cue')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,2)
plot(ts1,mean(master_cRew_cueBase_STREAMz,1,'omitnan'),'b')
hold on
cRew_std = std(master_cRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_cRew_cueBase_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cRew_cueBase_STREAMz,1,'omitnan') + cRew_std, fliplr(mean(master_cRew_cueBase_STREAMz,1,'omitnan') - cRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Correct Reward (F0 = Cue)')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,3)
plot(ts1,mean(master_cNoRew_cueBase_STREAMz,1,'omitnan'),'b')
hold on
cNoRew_std = std(master_cNoRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_cNoRew_cueBase_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_cNoRew_cueBase_STREAMz,1,'omitnan') + cNoRew_std, fliplr(mean(master_cNoRew_cueBase_STREAMz,1,'omitnan') - cNoRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Correct No Reward (F0 = Cue)')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,4)
plot(ts1,mean(master_iRew_cueBase_STREAMz,1,'omitnan'),'b')
hold on
iRew_std = std(master_iRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_iRew_cueBase_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_iRew_cueBase_STREAMz,1,'omitnan') + iRew_std, fliplr(mean(master_iRew_cueBase_STREAMz,1,'omitnan') - iRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Incorrect Reward (F0 = Cue)')
ylabel('GrabDA dF/F (z-Score)')

subplot(5,1,5)
plot(ts1,mean(master_iNoRew_cueBase_STREAMz,1,'omitnan'),'b')
hold on
iNoRew_std = std(master_iNoRew_cueBase_STREAMz,'omitnan')/sqrt(height(master_iNoRew_cueBase_STREAMz));
XX = [ts1, fliplr(ts1)];
YY = [mean(master_iNoRew_cueBase_STREAMz,1,'omitnan') + iNoRew_std, fliplr(mean(master_iNoRew_cueBase_STREAMz,1,'omitnan') - iNoRew_std)];
% Plot filled standard error bands.
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')
% Plot line at x=0
xline(0,'LineWidth',2,'Color','black')
title('Incorrect No Reward (F0 = Cue)')
ylabel('GrabDA dF/F (z-Score)')
xlabel('Time (s)')
toc

disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);