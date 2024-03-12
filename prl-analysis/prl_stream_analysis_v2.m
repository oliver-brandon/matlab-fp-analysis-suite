clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-5 -2]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0 5]; % time window to grab amplitude from
auc_window = [0 5];
tau_window = [0 2];
t = 5; % seconds to clip from start of streams
N = 10; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
baseAdjust = -1; % adjust baseline of signals to 'baseAdjust' seconds
toPlot = 0; % 1 = plot figures, 0 = don't plot
dualFiber = 0; % 1 = dual fiber, 0 = single fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir(...
    '/Users/brandon/My Drive (bloliv95@gmail.com)/prl/GrabDA/dls_jzl/mats (processed)/','Choose the .mat files you want to analyze.'...
    ); %gets directory%
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
phaseList = {};
treatList = {};

AMPdFF_analysis = cell(numFiles,7);
AMPz_analysis = cell(numFiles,7);
AUCdFF_analysis = cell(numFiles,7);
AUCz_analysis = cell(numFiles,7);
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
    fprintf('Analyzing %s (%d of %d)\n', name, i, numFiles)
    brokenID = strsplit(name,'_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    phaseList{i} = cellstr(strtrim(brokenID{2}));
    treatList{i} = cellstr(strtrim(brokenID{3}));

    load(filename)
    if dualFiber == 1
        if isfield(data.streams,'x405A')
            ISOS = 'x405A';
            SIGNAL = 'x465A';
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
        elseif isfield(data.streams,'x405C')
            ISOS = 'x405C';
            SIGNAL = 'x465C';
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
        else
            disp('No streams available.')
            break
        end
        if ~isfield(data.epocs,'CL1_')
            Correct = 0;
        else
            Correct = data.epocs.CL1_.onset;
        end
        if ~isfield(data.epocs,'IL1_')
            Incorrect = 0;
        else
            Incorrect = data.epocs.IL1_.onset;
        end

        [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);

        cueSTREAMraw = zeros(height(cue),epocArrayLen);
        cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);
        CorrectSTREAMraw = zeros(height(Correct),epocArrayLen);
        IncorrectSTREAMraw = zeros(height(Incorrect),epocArrayLen);

        cueSTREAMdFF = zeros(height(cue),epocArrayLen);
        cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
        CorrectSTREAMdFF = zeros(height(Correct),epocArrayLen);
        IncorrectSTREAMdFF = zeros(height(Incorrect),epocArrayLen);
      
        
        cueSTREAMz = zeros(height(cue),epocArrayLen);
        cRewSTREAMz = zeros(height(cRew),epocArrayLen);
        cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
        iRewSTREAMz = zeros(height(iRew),epocArrayLen);
        iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);
        CorrectSTREAMz = zeros(height(Correct),epocArrayLen);
        IncorrectSTREAMz = zeros(height(Incorrect),epocArrayLen);

        cueAMPdFF = zeros(height(cue),1);
        cRewAMPdFF = zeros(height(cRew),1);
        cNoRewAMPdFF = zeros(height(cNoRew),1);
        iRewAMPdFF = zeros(height(iRew),1);
        iNoRewAMPdFF = zeros(height(iNoRew),1);
        CorrectAMPdFF = zeros(height(Correct),1);
        IncorrectAMPdFF = zeros(height(Incorrect),1);

        cueAMPz = zeros(height(cue),1);
        cRewAMPz = zeros(height(cRew),1);
        cNoRewAMPz = zeros(height(cNoRew),1);
        iRewAMPz = zeros(height(iRew),1);
        iNoRewAMPz = zeros(height(iNoRew),1);
        CorrectAMPz = zeros(height(Correct),1);
        IncorrectAMPz = zeros(height(Incorrect),1);

        cueAUCdFF = zeros(height(cue),1);
        cRewAUCdFF = zeros(height(cRew),1);
        cNoRewAUCdFF = zeros(height(cNoRew),1);
        iRewAUCdFF = zeros(height(iRew),1);
        iNoRewAUCdFF = zeros(height(iNoRew),1);
        CorrectAUCdFF = zeros(height(Correct),1);
        IncorrectAUCdFF = zeros(height(Incorrect),1);
        
        cueAUCz = zeros(height(cue),1);
        cRewAUCz = zeros(height(cRew),1);
        cNoRewAUCz = zeros(height(cNoRew),1);
        iRewAUCz = zeros(height(iRew),1);
        iNoRewAUCz = zeros(height(iNoRew),1);
        CorrectAUCz = zeros(height(Correct),1);
        IncorrectAUCz = zeros(height(Incorrect),1);

        epocList = {cue;cRew;cNoRew;iRew;iNoRew;Correct;Incorrect};
        outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
            iRewSTREAMraw;iNoRewSTREAMraw;CorrectSTREAMraw;IncorrectSTREAMraw};
        outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
            iRewSTREAMdFF;iNoRewSTREAMdFF;CorrectSTREAMdFF;IncorrectSTREAMdFF};
        outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
            iRewSTREAMz;iNoRewSTREAMz;CorrectSTREAMz;IncorrectSTREAMz};
        outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF;CorrectAMPdFF;IncorrectAMPdFF};
        outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz;CorrectAMPz;IncorrectAMPz};
        outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF;CorrectAUCdFF;IncorrectAUCdFF};
        outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz;CorrectAUCz;IncorrectAUCz};
    else
        if isfield(data.streams, 'x405A')
            ISOS = 'x405A';
            SIGNAL = 'x465A';
    
            cue = data.epocs.St1_.onset;
            cRew = data.epocs.cRewA.onset;
            cNoRew = data.epocs.cNoRewA.onset;
            iRew = data.epocs.iRewA.onset;
            iNoRew = data.epocs.iNoRewA.onset;
            if ~isfield(data.epocs,'CL1_')
                Correct = 0;
            else
                Correct = data.epocs.CL1_.onset;
            end
            if ~isfield(data.epocs,'IL1_')
                Incorrect = 0;
            else
                Incorrect = data.epocs.IL1_.onset;
            end
    
            cueSTREAMraw = zeros(height(cue),epocArrayLen);
            cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMraw = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMraw = zeros(height(Incorrect),epocArrayLen);
    
            cueSTREAMdFF = zeros(height(cue),epocArrayLen);
            cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMdFF = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMdFF = zeros(height(Incorrect),epocArrayLen);
          
            
            cueSTREAMz = zeros(height(cue),epocArrayLen);
            cRewSTREAMz = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMz = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMz = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMz = zeros(height(Incorrect),epocArrayLen);
    
            cueAMPdFF = zeros(height(cue),1);
            cRewAMPdFF = zeros(height(cRew),1);
            cNoRewAMPdFF = zeros(height(cNoRew),1);
            iRewAMPdFF = zeros(height(iRew),1);
            iNoRewAMPdFF = zeros(height(iNoRew),1);
            CorrectAMPdFF = zeros(height(Correct),1);
            IncorrectAMPdFF = zeros(height(Incorrect),1);
    
            cueAMPz = zeros(height(cue),1);
            cRewAMPz = zeros(height(cRew),1);
            cNoRewAMPz = zeros(height(cNoRew),1);
            iRewAMPz = zeros(height(iRew),1);
            iNoRewAMPz = zeros(height(iNoRew),1);
            CorrectAMPz = zeros(height(Correct),1);
            IncorrectAMPz = zeros(height(Incorrect),1);
    
            cueAUCdFF = zeros(height(cue),1);
            cRewAUCdFF = zeros(height(cRew),1);
            cNoRewAUCdFF = zeros(height(cNoRew),1);
            iRewAUCdFF = zeros(height(iRew),1);
            iNoRewAUCdFF = zeros(height(iNoRew),1);
            CorrectAUCdFF = zeros(height(Correct),1);
            IncorrectAUCdFF = zeros(height(Incorrect),1);
            
            cueAUCz = zeros(height(cue),1);
            cRewAUCz = zeros(height(cRew),1);
            cNoRewAUCz = zeros(height(cNoRew),1);
            iRewAUCz = zeros(height(iRew),1);
            iNoRewAUCz = zeros(height(iNoRew),1);
            CorrectAUCz = zeros(height(Correct),1);
            IncorrectAUCz = zeros(height(Incorrect),1);
    
    
    
            [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
                cNoRew,iRew,iNoRew);
            epocList = {cue;cRew;cNoRew;iRew;iNoRew;Correct;Incorrect};
            outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
                iRewSTREAMraw;iNoRewSTREAMraw;CorrectSTREAMraw;IncorrectSTREAMraw};
            outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
                iRewSTREAMdFF;iNoRewSTREAMdFF;CorrectSTREAMdFF;IncorrectSTREAMdFF};
            outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
                iRewSTREAMz;iNoRewSTREAMz;CorrectSTREAMz;IncorrectSTREAMz};
            outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF;CorrectAMPdFF;IncorrectAMPdFF};
            outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz;CorrectAMPz;IncorrectAMPz};
            outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF;CorrectAUCdFF;IncorrectAUCdFF};
            outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz;CorrectAUCz;IncorrectAUCz};
             
            data.analysis.sessionID = session_identifiers;
            
        elseif isfield(data.streams, 'x405C')
            ISOS = 'x405C';
            SIGNAL = 'x465C';
    
            cue = data.epocs.St2_.onset;
            cRew = data.epocs.cRewC.onset;
            cNoRew = data.epocs.cNoRewC.onset;
            iRew = data.epocs.iRewC.onset;
            iNoRew = data.epocs.iNoRewC.onset;
            if ~isfield(data.epocs,'CL2_')
                Correct = 0;
            else
                Correct = data.epocs.CL2_.onset;
            end
            if ~isfield(data.epocs,'IL2_')
                Incorrect = 0;
            else
                Incorrect = data.epocs.IL2_.onset;
            end
    
            cueSTREAMraw = zeros(height(cue),epocArrayLen);
            cRewSTREAMraw = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMraw = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMraw = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMraw = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMraw = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMraw = zeros(height(Incorrect),epocArrayLen);
    
            cueSTREAMdFF = zeros(height(cue),epocArrayLen);
            cRewSTREAMdFF = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMdFF = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMdFF = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMdFF = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMdFF = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMdFF = zeros(height(Incorrect),epocArrayLen);
          
            
            cueSTREAMz = zeros(height(cue),epocArrayLen);
            cRewSTREAMz = zeros(height(cRew),epocArrayLen);
            cNoRewSTREAMz = zeros(height(cNoRew),epocArrayLen);
            iRewSTREAMz = zeros(height(iRew),epocArrayLen);
            iNoRewSTREAMz = zeros(height(iNoRew),epocArrayLen);
            CorrectSTREAMz = zeros(height(Correct),epocArrayLen);
            IncorrectSTREAMz = zeros(height(Incorrect),epocArrayLen);
    
            cueAMPdFF = zeros(height(cue),1);
            cRewAMPdFF = zeros(height(cRew),1);
            cNoRewAMPdFF = zeros(height(cNoRew),1);
            iRewAMPdFF = zeros(height(iRew),1);
            iNoRewAMPdFF = zeros(height(iNoRew),1);
            CorrectAMPdFF = zeros(height(Correct),1);
            IncorrectAMPdFF = zeros(height(Incorrect),1);
    
            cueAMPz = zeros(height(cue),1);
            cRewAMPz = zeros(height(cRew),1);
            cNoRewAMPz = zeros(height(cNoRew),1);
            iRewAMPz = zeros(height(iRew),1);
            iNoRewAMPz = zeros(height(iNoRew),1);
            CorrectAMPz = zeros(height(Correct),1);
            IncorrectAMPz = zeros(height(Incorrect),1);
    
            cueAUCdFF = zeros(height(cue),1);
            cRewAUCdFF = zeros(height(cRew),1);
            cNoRewAUCdFF = zeros(height(cNoRew),1);
            iRewAUCdFF = zeros(height(iRew),1);
            iNoRewAUCdFF = zeros(height(iNoRew),1);
            CorrectAUCdFF = zeros(height(Correct),1);
            IncorrectAUCdFF = zeros(height(Incorrect),1);
            
            cueAUCz = zeros(height(cue),1);
            cRewAUCz = zeros(height(cRew),1);
            cNoRewAUCz = zeros(height(cNoRew),1);
            iRewAUCz = zeros(height(iRew),1);
            iNoRewAUCz = zeros(height(iNoRew),1);
            CorrectAUCz = zeros(height(Correct),1);
            IncorrectAUCz = zeros(height(Incorrect),1);
    
            [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
                cNoRew,iRew,iNoRew);
            epocList = {cue;cRew;cNoRew;iRew;iNoRew;Correct;Incorrect};
            outputSTREAMSraw = {cueSTREAMraw;cRewSTREAMraw;cNoRewSTREAMraw;...
                iRewSTREAMraw;iNoRewSTREAMraw;CorrectSTREAMraw;IncorrectSTREAMraw};
            outputSTREAMSdFF = {cueSTREAMdFF;cRewSTREAMdFF;cNoRewSTREAMdFF;...
                iRewSTREAMdFF;iNoRewSTREAMdFF;CorrectSTREAMdFF;IncorrectSTREAMdFF};
            outputSTREAMSz = {cueSTREAMz;cRewSTREAMz;cNoRewSTREAMz;...
                iRewSTREAMz;iNoRewSTREAMz;CorrectSTREAMz;IncorrectSTREAMz};
            outputAMPdFF = {cueAMPdFF;cRewAMPdFF;cNoRewAMPdFF;iRewAMPdFF;iNoRewAMPdFF;CorrectAMPdFF;IncorrectAMPdFF};
            outputAMPz = {cueAMPz;cRewAMPz;cNoRewAMPz;iRewAMPz;iNoRewAMPz;CorrectAMPz;IncorrectAMPz};
            outputAUCdFF = {cueAUCdFF;cRewAUCdFF;cNoRewAUCdFF;iRewAUCdFF;iNoRewAUCdFF;CorrectAUCdFF;IncorrectAUCdFF};
            outputAUCz = {cueAUCz;cRewAUCz;cNoRewAUCz;iRewAUCz;iNoRewAUCz;CorrectAUCz;IncorrectAUCz};
            
            data.analysis.sessionID = session_identifiers;
        end
    end
    %time array used for all streams%
    session_time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
    ind = find(session_time>t,1);% find first index of when time crosses threshold
    session_time = session_time(ind:end); % reformat vector to only include allowed time
    SIGNAL_raw = data.streams.(SIGNAL).data(ind:end);
    ISOS_raw = data.streams.(ISOS).data(ind:end);
   
    %downsample streams and time array by N times%
    ISOS_raw = downsample(ISOS_raw, N);
    SIGNAL_raw = downsample(SIGNAL_raw, N);
    minStreamLength = min(length(ISOS_raw),length(SIGNAL_raw));
    ISOS_raw = ISOS_raw(1:minStreamLength);
    SIGNAL_raw = SIGNAL_raw(1:minStreamLength);
    
    session_time = downsample(session_time, N);
    ts1 = -baseWindow + (1:epocArrayLen) / data.streams.(SIGNAL).fs*N;
    % establish baseline windows
    [~,baseSt] = min(abs(ts1 - (baseline(1))));
    [~,baseEn] = min(abs(ts1 - (baseline(2))));
    [~,ampSt] = min(abs(ts1 - (amp_window(1))));
    [~,ampEn] = min(abs(ts1 - (amp_window(2))));
    [~,aucSt] = min(abs(ts1 - (auc_window(1))));
    [~,aucEn] = min(abs(ts1 - (auc_window(2))));
    [~,tauSt] = min(abs(ts1 - (tau_window(1))));
    [~,tauEn] = min(abs(ts1 - (tau_window(2))));
    
    %% Streams baselined to cue preceding it %%
    idx = find(ts1>baseAdjust,1);
    cueArray = session_identifiers(1:2:end,:);
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(1:end-1,:);
        cueArray = cueArray(1:end-1,1);
    end
    leverArray = session_identifiers(2:2:end,:);
    levers_z = zeros(height(leverArray),epocArrayLen);
    levers_raw = zeros(height(leverArray),epocArrayLen);
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
        levers_raw(m,:) = leverSigRaw;
                
    end
    sessionSTREAMSraw{i,1} = levers_raw;
    amp_cueBase = [];
    auc_cueBase = [];
    for n = 1:height(levers_raw)
        % dF/F
        meanCue = cueBaseMean(n,1);
        stdCue = cueBaseStd(n,1);
        
        levers_dFF(n,1:epocArrayLen) = levers_raw(n, 1:epocArrayLen) - meanCue;
        levers_dFF(n,1:epocArrayLen) = 100*(levers_dFF(n,1:epocArrayLen) / meanCue);
        % z-Score
        meanCueBase_dFF = mean(levers_dFF(n,baseSt:baseEn));
        stdCueBase_dFF = std(levers_dFF(n,baseSt:baseEn));
        levers_z(n,1:epocArrayLen) = (levers_dFF(n,1:epocArrayLen) - meanCueBase_dFF) / stdCueBase_dFF;
        % adjusts streams to baseline of zero at -0.5s %
        % if levers_z(n,idx) < 0
        %     val = levers_z(n,idx);
        %     diff = 0 - val;
        %     levers_z(n,1:epocArrayLen) = levers_z(n,1:epocArrayLen) + abs(diff);
        % elseif levers_z(n,idx) > 0
        %     val = levers_z(n,idx);
        %     diff = 0 - val;
        %     levers_z(n,1:epocArrayLen) = levers_z(n,1:epocArrayLen) - abs(diff);
        % end
        % amp_cueBase(n,1) = max(levers_z(n,ampSt:ampEn));
        % % Calculate AUC above x=0 %
        % positive_indices = levers_z(n,:) > 0;
        % positive_indices = positive_indices(aucSt:aucEn);
        % numPos = sum(positive_indices);
        % if numPos < 3
        %     auc_cueBase(n,1) = 0;
        %     continue
        % end
        % 
        % ts2 = linspace(auc_window(1),auc_window(2),length(positive_indices));
        % y_pos = levers_z(n,positive_indices);
        % x_pos = ts2(1,positive_indices);
        % auc_cueBase(n,1) = trapz(x_pos,y_pos);
        % auc_cueBase(n,1) = trapz(levers_z(n,aucSt:aucEn),ts1(1,aucSt:aucEn));
    end

    sessionSTREAMSz{i,1} = mean(levers_z);
    trialNames = array2table(session_identifiers(2:2:end,2),'VariableNames',{'Trial_Type'});
    levers_z = array2table(levers_z);
    levers_z = horzcat(trialNames,levers_z);
    cRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 1, 2:end));
    cNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 2, 2:end));
    iRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 3, 2:end));
    iNoRew_cueBase = table2array(levers_z(levers_z.Trial_Type == 4, 2:end));

    amp_cueBase = array2table(amp_cueBase);
    amp_cueBase = horzcat(trialNames,amp_cueBase);
    auc_cueBase = array2table(auc_cueBase);
    auc_cueBase = horzcat(trialNames,auc_cueBase);
        
    if isempty(cRew_cueBase)
        amp_cRew_cueBase(i,:) = nan;
        auc_cRew_cueBase(i,:) = nan;
        tau_cRew_cueBase(i,:) = nan;
    else
        % amp_cRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        % auc_cRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
        amp_cRew_cueBase(i,1) = calculateAMP(mean(cRew_cueBase(:,ampSt:ampEn),1,'omitnan'));
        auc_cRew_cueBase(i,1) = calculateAUC(mean(cRew_cueBase(:,aucSt:aucEn),1,'omitnan'),ts1(:,aucSt:aucEn));
        tau_cRew_cueBase(i,1) = calculateTAU(mean(cRew_cueBase(:,tauSt:tauEn),1,'omitnan'),ts1(:,tauSt:tauEn));
    end

    if isempty(cNoRew_cueBase)
        amp_cNoRew_cueBase = nan;
        auc_cNoRew_cueBase = nan;
        tau_cNoRew_cueBase = nan;
    else
        % amp_cNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        % auc_cNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
        amp_cNoRew_cueBase(i,1) = calculateAMP(mean(cNoRew_cueBase(:,ampSt:ampEn),1,'omitnan'));
        auc_cNoRew_cueBase(i,1) = calculateAUC(mean(cNoRew_cueBase(:,aucSt:aucEn),1,'omitnan'),ts1(:,aucSt:aucEn));
        tau_cNoRew_cueBase(i,1) = calculateTAU(mean(cNoRew_cueBase(:,tauSt:tauEn),1,'omitnan'),ts1(:,tauSt:tauEn));
    end

    if isempty(iRew_cueBase)
        amp_iRew_cueBase = nan;
        auc_iRew_cueBase = nan;
        tau_iRew_cueBase = nan;
    else
        % amp_iRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        % auc_iRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
        amp_iRew_cueBase(i,1) = calculateAMP(mean(iRew_cueBase(:,ampSt:ampEn),1,'omitnan'));
        auc_iRew_cueBase(i,1) = calculateAUC(mean(iRew_cueBase(:,aucSt:aucEn),1,'omitnan'),ts1(:,aucSt:aucEn));
        tau_iRew_cueBase(i,1) = calculateTAU(mean(iRew_cueBase(:,tauSt:tauEn),1,'omitnan'),ts1(:,tauSt:tauEn));
    end

    if isempty(iNoRew_cueBase)
        amp_iNoRew_cueBase = nan;
        auc_iNoRew_cueBase = nan;
        tau_iNoRew_cueBase = nan;
    else
        % amp_iNoRew_cueBase = table2array(amp_cueBase(amp_cueBase.Trial_Type == 1,2));
        % auc_iNoRew_cueBase = table2array(auc_cueBase(auc_cueBase.Trial_Type == 1,2));
        amp_iNoRew_cueBase(i,1) = calculateAMP(mean(iNoRew_cueBase(:,ampSt:ampEn),1,'omitnan'));
        auc_iNoRew_cueBase(i,1) = calculateAUC(mean(iNoRew_cueBase(:,aucSt:aucEn),1,'omitnan'),ts1(:,aucSt:aucEn));
        tau_iNoRew_cueBase(i,1) = calculateTAU(mean(iNoRew_cueBase(:,tauSt:tauEn),1,'omitnan'),ts1(:,tauSt:tauEn));
    end

    %% Streams baslined to before each epoc %%
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
                ampdFF(ii) = NaN;
                streams_z(ii,1:epocArrayLen) = NaN;
                ampZ(ii) = NaN;
                aucdFF(ii) = NaN;
                aucZ(ii) = NaN;
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

        [~,baseSt] = min(abs(ts1 - (baseline(1))));
        [~,baseEn] = min(abs(ts1 - (baseline(2))));

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
            % adjusts streams to baseline of zero at -0.5s %
            if streams_z(j,idx) < 0
                val = streams_z(j,idx);
                diff = 0 - val;
                streams_z(j,1:epocArrayLen) = streams_z(j,1:epocArrayLen) + abs(diff);
            elseif streams_z(j,idx) > 0
                val = streams_z(j,idx);
                diff = 0 - val;
                streams_z(j,1:epocArrayLen) = streams_z(j,1:epocArrayLen) - abs(diff);
            end

            % amplitude
            ampdFF(j) = max(streams_dFF(j,ampSt:ampEn));
            ampZ(j) = max(streams_z(j,ampSt:ampEn));

            % % Calculate AUC above x=0 (Z-score) %
            % positive_indices = streams_dFF(j,:) > 0;
            % y_pos = streams_dFF(j,positive_indices);
            % x_pos = ts1(1,positive_indices);
            % aucZ(j) = trapz(x_pos,y_pos);
            aucZ(j) = calculateAUC(streams_z(j,aucSt:aucEn),ts1(1,aucSt:aucEn));
            % % Calculate AUC above x=0 (dFF) %
            % positive_indices = streams_z(j,:) > 0;
            % y_pos = streams_z(j,positive_indices);
            % x_pos = ts1(1,positive_indices);
            % aucdFF(j) = trapz(x_pos,y_pos);
            % aucdFF(j) = trapz(streams_dFF(j,aucSt:aucEn),ts1(1,aucSt:aucEn));
            
        
            outputSTREAMSraw{k,1} = streams_raw(:,1:epocArrayLen);
            outputSTREAMSdFF{k,1} = streams_dFF(:,1:epocArrayLen);
            outputSTREAMSz{k,1} = streams_z(:,1:epocArrayLen);
            outputAMPdFF{k,1} = ampdFF;
            outputAMPz{k,1} = ampZ;
            outputAUCdFF{k,1} = aucdFF;
            outputAUCz{k,1} = aucZ;
        end
        cue_cueBase = outputSTREAMSz{1,1};
        amp_cue_cueBase(i,1) = calculateAMP(mean(cue_cueBase(:,ampSt:ampEn),1,'omitnan'));
        auc_cue_cueBase(i,1) = calculateAUC(mean(cue_cueBase(:,aucSt:aucEn),1,'omitnan'),ts1(:,aucSt:aucEn));
        tau_cue_cueBase(i,1) = calculateTAU(mean(cue_cueBase(:,tauSt:tauEn),1,'omitnan'),ts1(:,tauSt:tauEn));
    end

    



    master_cue_STREAMraw(i,:) = mean(outputSTREAMSraw{1,1},1);
    master_cRew_STREAMraw(i,:) = mean(outputSTREAMSraw{2,1},1);
    master_cNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{3,1},1);
    master_iRew_STREAMraw(i,:) = mean(outputSTREAMSraw{4,1},1);
    master_iNoRew_STREAMraw(i,:) = mean(outputSTREAMSraw{5,1},1);
    master_correct_STREAMraw(i,:) = mean(outputSTREAMSraw{6,1},1);
    master_incorrect_STREAMraw(i,:) = mean(outputSTREAMSraw{7,1},1);

    master_cue_STREAMdFF(i,:) = mean(outputSTREAMSdFF{1,1},1);
    master_cRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{2,1},1);
    master_cNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{3,1},1);
    master_iRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{4,1},1);
    master_iNoRew_STREAMdFF(i,:) = mean(outputSTREAMSdFF{5,1},1);
    master_correct_STREAMdFF(i,:) = mean(outputSTREAMSdFF{6,1},1);
    master_incorrect_STREAMdFF(i,:) = mean(outputSTREAMSdFF{7,1},1);

    master_cue_STREAMz(i,:) = mean(outputSTREAMSz{1,1},1);
    master_cRew_STREAMz(i,:) = mean(outputSTREAMSz{2,1},1);
    master_cNoRew_STREAMz(i,:) = mean(outputSTREAMSz{3,1},1);
    master_iRew_STREAMz(i,:) = mean(outputSTREAMSz{4,1},1);
    master_iNoRew_STREAMz(i,:) = mean(outputSTREAMSz{5,1},1);
    master_correct_STREAMz(i,:) = mean(outputSTREAMSz{6,1},1);
    master_incorrect_STREAMz(i,:) = mean(outputSTREAMSz{7,1},1);

    master_cRew_cueBase_STREAMz(i,:) = mean(cRew_cueBase,1);
    master_cNoRew_cueBase_STREAMz(i,:) = mean(cNoRew_cueBase,1);
    master_iRew_cueBase_STREAMz(i,:) = mean(iRew_cueBase,1);
    master_iNoRew_cueBase_STREAMz(i,:) = mean(iNoRew_cueBase,1);

    % Uncomment if grabbing amp, auc, tau from each trial
    % amp_cueBase_analysis(i,1:5) = {mean(amp_cue_cueBase,'omitnan') mean(amp_cRew_cueBase,'omitnan')...
    %     mean(amp_cNoRew_cueBase,'omitnan') mean(amp_iRew_cueBase,'omitnan') mean(amp_iNoRew_cueBase,'omitnan')};
    % auc_cueBase_analysis(i,1:5) = {mean(auc_cue_cueBase,'omitnan') mean(auc_cRew_cueBase,'omitnan')...
    %     mean(auc_cNoRew_cueBase,'omitnan') mean(auc_iRew_cueBase,'omitnan') mean(auc_iNoRew_cueBase,'omitnan')};
    % tau_cueBase_analysis(i,1:5) = {mean(tau_cue_cueBase,'omitnan') mean(tau_cRew_cueBase,'omitnan')...
    %     mean(tau_cNoRew_cueBase,'omitnan') mean(tau_iRew_cueBase,'omitnan') mean(tau_iNoRew_cueBase,'omitnan')};
    
    amp_cueBase_analysis(i,1:5) = {mean(amp_cue_cueBase,1,'omitnan') mean(amp_cRew_cueBase,1,'omitnan')...
        mean(amp_cNoRew_cueBase,1,'omitnan') mean(amp_iRew_cueBase,1,'omitnan') mean(amp_iNoRew_cueBase,1,'omitnan')};
    auc_cueBase_analysis(i,1:5) = {mean(auc_cue_cueBase,1,'omitnan') mean(auc_cRew_cueBase,1,'omitnan')...
        mean(auc_cNoRew_cueBase,1,'omitnan') mean(auc_iRew_cueBase,1,'omitnan') mean(auc_iNoRew_cueBase,1,'omitnan')};
    tau_cueBase_analysis(i,1:5) = {mean(tau_cue_cueBase,1,'omitnan') mean(tau_cRew_cueBase,1,'omitnan')...
        mean(tau_cNoRew_cueBase,1,'omitnan') mean(tau_iRew_cueBase,1,'omitnan') mean(tau_iNoRew_cueBase,1,'omitnan')};
   
    
    AMPdFF_analysis(i,1:7) = {mean(outputAMPdFF{1},'omitnan') mean(outputAMPdFF{2},'omitnan')...
        mean(outputAMPdFF{3},'omitnan') mean(outputAMPdFF{4},'omitnan') mean(outputAMPdFF{5},'omitnan') mean(outputAMPdFF{6},'omitnan')... 
        mean(outputAMPdFF{7},'omitnan')};
    AMPz_analysis(i,1:7) = {mean(outputAMPz{1},'omitnan') mean(outputAMPz{2},'omitnan')...
        mean(outputAMPz{3},'omitnan') mean(outputAMPz{4},'omitnan') mean(outputAMPz{5},'omitnan') mean(outputAMPz{6},'omitnan')... 
        mean(outputAMPz{7},'omitnan')};

    AUCdFF_analysis(i,1:7) = {mean(outputAUCdFF{1},'omitnan') mean(outputAUCdFF{2},'omitnan')...
        mean(outputAUCdFF{3},'omitnan') mean(outputAUCdFF{4},'omitnan') mean(outputAUCdFF{5},'omitnan') mean(outputAUCdFF{6},'omitnan')... 
        mean(outputAUCdFF{7},'omitnan')};
    AUCz_analysis(i,1:7) = {mean(outputAUCz{1},'omitnan') mean(outputAUCz{2},'omitnan')...
        mean(outputAUCz{3},'omitnan') mean(outputAUCz{4},'omitnan') mean(outputAUCz{5},'omitnan') mean(outputAUCz{6},'omitnan')... 
        mean(outputAUCz{7},'omitnan')};

  
    % save(filename,'data');
        
end
idList = cell2table(IDs','VariableNames',{'ID'});
phaseList = cell2table(phaseList','VariableNames',{'Phase'});
treatList = cell2table(treatList','VariableNames',{'Treatment'});


%% Amplitude Table %%
AMPz_analysis_table = cell2table(AMPz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew','Correct Lever','Incorrect Lever'});
AMPz_analysis_table = horzcat(idList,treatList,phaseList,AMPz_analysis_table);
AMPz_analysis_table = sortrows(AMPz_analysis_table,{'Phase','Treatment'},{'ascend','descend'});
amp_cueBase_table = cell2table(amp_cueBase_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
amp_cueBase_table = horzcat(idList,treatList,phaseList,amp_cueBase_table);
%% Area Under Curve Table %%
AUCz_analysis_table = cell2table(AUCz_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew','Correct Lever','Incorrect Lever'});
AUCz_analysis_table = horzcat(idList,treatList,phaseList,AUCz_analysis_table);
AUCz_analysis_table = sortrows(AUCz_analysis_table,{'Phase','Treatment'},{'ascend','descend'});
auc_cueBase_table = cell2table(auc_cueBase_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
auc_cueBase_table = horzcat(idList,treatList,phaseList,auc_cueBase_table);
%% Tau Table %%
tau_cueBase_table = cell2table(tau_cueBase_analysis,'VariableNames',{'Cue','cRew','cNoRew','iRew','iNoRew'});
tau_cueBase_table = horzcat(idList,treatList,phaseList,tau_cueBase_table);
%% Correct/Incorrect Stream Table %%
master_correct_STREAMz = array2table(master_correct_STREAMz);
master_correct_STREAMz = horzcat(idList,treatList,phaseList,master_correct_STREAMz);
master_correct_STREAMz = sortrows(master_correct_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
master_incorrect_STREAMz = array2table(master_incorrect_STREAMz);
master_incorrect_STREAMz = horzcat(idList,treatList,phaseList,master_incorrect_STREAMz);
master_incorrect_STREAMz = sortrows(master_incorrect_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
%% Trial Type Stream Tables %%
master_cRew_STREAMz = array2table(master_cRew_STREAMz);
master_cRew_STREAMz = horzcat(idList,treatList,phaseList,master_cRew_STREAMz);
master_cRew_STREAMz = sortrows(master_cRew_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
master_cNoRew_STREAMz = array2table(master_cNoRew_STREAMz);
master_cNoRew_STREAMz = horzcat(idList,treatList,phaseList,master_cNoRew_STREAMz);
master_cNoRew_STREAMz = sortrows(master_cNoRew_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
master_iRew_STREAMz = array2table(master_iRew_STREAMz);
master_iRew_STREAMz = horzcat(idList,treatList,phaseList,master_iRew_STREAMz);
master_iRew_STREAMz = sortrows(master_iRew_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
master_iNoRew_STREAMz = array2table(master_iNoRew_STREAMz);
master_iNoRew_STREAMz = horzcat(idList,treatList,phaseList,master_iNoRew_STREAMz);
master_iNoRew_STREAMz = sortrows(master_iNoRew_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
%% Epoc Stream Tables Baselined to Cue %%
a_cue_STREAMz = array2table(master_cue_STREAMz);
a_cue_STREAMz = horzcat(idList,phaseList,treatList,a_cue_STREAMz);
a_cRew_cueBase_STREAMz = array2table(master_cRew_cueBase_STREAMz);
a_cRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_cRew_cueBase_STREAMz);
a_cNoRew_cueBase_STREAMz = array2table(master_cNoRew_cueBase_STREAMz);
a_cNoRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_cNoRew_cueBase_STREAMz);
a_iRew_cueBase_STREAMz = array2table(master_iRew_cueBase_STREAMz);
a_iRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_iRew_cueBase_STREAMz);
a_iNoRew_cueBase_STREAMz = array2table(master_iNoRew_cueBase_STREAMz);
a_iNoRew_cueBase_STREAMz = horzcat(idList,phaseList,treatList,a_iNoRew_cueBase_STREAMz);

a_cue_STREAMz = sortrows(a_cue_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_cRew_cueBase_STREAMz = sortrows(a_cRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_cNoRew_cueBase_STREAMz = sortrows(a_cNoRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_iRew_cueBase_STREAMz = sortrows(a_iRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});
a_iNoRew_cueBase_STREAMz = sortrows(a_iNoRew_cueBase_STREAMz,{'Phase','Treatment'},{'ascend','descend'});


if toPlot == 1
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
    % sgtitle(treatment)
    
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
    % sgtitle(treatment)
    
end
amp_cueBase_table = sortrows(amp_cueBase_table,{'Phase','Treatment'},{'ascend','descend'});
auc_cueBase_table = sortrows(auc_cueBase_table,{'Phase','Treatment'},{'ascend','descend'});
tau_cueBase_table = sortrows(tau_cueBase_table,{'Phase','Treatment'},{'ascend','descend'});

prl_stream_analysis.acq1.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Acq1'),:);
prl_stream_analysis.acq2.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Acq2'),:);
prl_stream_analysis.rev1.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev1'),:);
prl_stream_analysis.rev2.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev2'),:);
prl_stream_analysis.rev3.amplitude = amp_cueBase_table(strcmp(amp_cueBase_table.Phase,'Rev3'),:);
prl_stream_analysis.acq1.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Acq1'),:);
prl_stream_analysis.acq2.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Acq2'),:);
prl_stream_analysis.rev1.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev1'),:);
prl_stream_analysis.rev2.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev2'),:);
prl_stream_analysis.rev3.auc = auc_cueBase_table(strcmp(auc_cueBase_table.Phase,'Rev3'),:);
prl_stream_analysis.acq1.tau = tau_cueBase_table(strcmp(tau_cueBase_table.Phase,'Acq1'),:);
prl_stream_analysis.acq2.tau = tau_cueBase_table(strcmp(tau_cueBase_table.Phase,'Acq2'),:);
prl_stream_analysis.rev1.tau = tau_cueBase_table(strcmp(tau_cueBase_table.Phase,'Rev1'),:);
prl_stream_analysis.rev2.tau = tau_cueBase_table(strcmp(tau_cueBase_table.Phase,'Rev2'),:);
prl_stream_analysis.rev3.tau = tau_cueBase_table(strcmp(tau_cueBase_table.Phase,'Rev3'),:);

prl_stream_analysis.acq1.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.correct = master_correct_STREAMz(strcmp(master_correct_STREAMz.Phase,'Acq1'),:);
prl_stream_analysis.acq1.incorrect = master_incorrect_STREAMz(strcmp(master_incorrect_STREAMz.Phase,'Acq1'),:);

prl_stream_analysis.acq2.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.correct = master_correct_STREAMz(strcmp(master_correct_STREAMz.Phase,'Acq2'),:);
prl_stream_analysis.acq2.incorrect = master_incorrect_STREAMz(strcmp(master_incorrect_STREAMz.Phase,'Acq2'),:);

prl_stream_analysis.rev1.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.correct = master_correct_STREAMz(strcmp(master_correct_STREAMz.Phase,'Rev1'),:);
prl_stream_analysis.rev1.incorrect = master_incorrect_STREAMz(strcmp(master_incorrect_STREAMz.Phase,'Rev1'),:);

prl_stream_analysis.rev2.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.correct = master_correct_STREAMz(strcmp(master_correct_STREAMz.Phase,'Rev2'),:);
prl_stream_analysis.rev2.incorrect = master_incorrect_STREAMz(strcmp(master_incorrect_STREAMz.Phase,'Rev2'),:);

prl_stream_analysis.rev3.cue = a_cue_STREAMz(strcmp(a_cue_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.cRew = a_cRew_cueBase_STREAMz(strcmp(a_cRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.cNoRew = a_cNoRew_cueBase_STREAMz(strcmp(a_cNoRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.iRew = a_iRew_cueBase_STREAMz(strcmp(a_iRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.iNoRew = a_iNoRew_cueBase_STREAMz(strcmp(a_iNoRew_cueBase_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.correct = master_correct_STREAMz(strcmp(master_correct_STREAMz.Phase,'Rev3'),:);
prl_stream_analysis.rev3.incorrect = master_incorrect_STREAMz(strcmp(master_incorrect_STREAMz.Phase,'Rev3'),:);

prl_stream_analysis.metadata.time = ts1;
prl_stream_analysis.metadata.timeWindow = timeWindow;
prl_stream_analysis.metadata.baseWindow = baseWindow;
prl_stream_analysis.metadata.baseline = baseline;
prl_stream_analysis.metadata.amp_window = amp_window;
prl_stream_analysis.metadata.auc_window = auc_window;
prl_stream_analysis.metadata.tau_window = tau_window;
prl_stream_analysis.metadata.trimStart = t;
prl_stream_analysis.metadata.downsample = N;
prl_stream_analysis.metadata.Hz = sigHz;
prl_stream_analysis.metadata.baseAdjust = baseAdjust;

% sfn2023.AMP = AMPz_analysis_table;
% sfn2023.AMP = sortrows(sfn2023.AMP,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.AUC = AUCz_analysis_table;
% sfn2023.AUC = sortrows(sfn2023.AUC,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.correctLever = master_correct_STREAMz;
% sfn2023.correctLever = sortrows(sfn2023.correctLever,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.incorrectLever = master_incorrect_STREAMz;
% sfn2023.incorrectLever = sortrows(sfn2023.incorrectLever,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.cue = a_cue_STREAMz;
% sfn2023.cue = sortrows(sfn2023.cue,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.cRew = master_cRew_STREAMz;
% sfn2023.cRew = sortrows(sfn2023.cRew,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.cNoRew = master_cNoRew_STREAMz;
% sfn2023.cNoRew = sortrows(sfn2023.cNoRew,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.iRew = master_iRew_STREAMz;
% sfn2023.iRew = sortrows(sfn2023.iRew,{'Phase','Treatment'},{'ascend','descend'});
% sfn2023.iNoRew = master_iNoRew_STREAMz;
% sfn2023.iNoRew = sortrows(sfn2023.iNoRew,{'Phase','Treatment'},{'ascend','descend'});

% save('../data-files/prl_stream_analysis.mat','prl_stream_analysis')

toc

disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);
clearvars -except prl_stream_analysis