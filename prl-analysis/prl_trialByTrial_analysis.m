clear; close all;
warning off
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseWindow = 5; % baseline signal to include before TTL 
baseline = [-3 -1]; % baseline signal for dFF/zscore (seconds before onset, positive integer)
amp_window = [0 1]; % time window to grab amplitude from
auc_window = [0 timeWindow];
baseAdjust = -0.5;
t = 5; % seconds to clip from start of streams
N = 100; %Downsample N times
sigHz = 1017/N;
epocArrayLen = round(sigHz * (timeWindow + baseWindow));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myDir = uigetdir(...
    '/Users/brandon/My Drive (bloliv95@gmail.com)/prl/PRL_GRABDA/test','Choose the .mat files you want to analyze.'...
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
idTable_cue = [];
idTable_lev = [];
phaseTable_cue = [];
phaseTable_lev = [];
treatmentTable_cue = [];
treatmentTable_lev = [];
epocAMP_cue = [];
epocAUC_cue = [];
epocAMP_lev = [];
epocAUC_lev = [];
trialByTrial_cue = [];
trialByTrial_lev = [];
trialNumCue = [];
trialNumLev = [];
trialTypeCue = [];
trialTypeLev = [];
for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    brokenID = strsplit(name,'_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    phaseList{i} = cellstr(strtrim(brokenID{2}));
    treatList{i} = cellstr(strtrim(brokenID{3}));

    load(filename)


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
    elseif isfield(data.streams,'x405C')
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
    end
    [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(cue,cRew,...
            cNoRew,iRew,iNoRew);

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

    idx = find(ts1>baseAdjust,1);
    cueArray = session_identifiers(1:2:end,:);
    if session_identifiers(end,2) == 0
        session_identifiers = session_identifiers(1:end-1,:);
        cueArray = cueArray(1:end-1,:);
    end
    leverArray = session_identifiers(2:2:end,:);
    epocSig_z = zeros(height(leverArray),epocArrayLen);
    epocSig_raw = zeros(height(leverArray),epocArrayLen);
    for j = 1:height(session_identifiers(:,1))
        

        epocStart = session_identifiers(j,1) - baseWindow;
        leverEnd = epocStart + timeWindow + baseWindow;
        [~,epoSt] = min(abs(session_time - epocStart));
        [~,epoEn] = min(abs(session_time - leverEnd));
        epocSnip = SIGNAL_raw(1,epoSt:epoEn);
        if length(epocSnip) < epocArrayLen
            mn = mean(epocSnip(1,end-10:end));
            epocSnip(1,end:epocArrayLen) = mn;
        elseif length(epocSnip) > epocArrayLen
            op = length(epocSnip);
            arrayDif = op - epocArrayLen;
            epocSnip = epocSnip(1,1:end-arrayDif);
        end
        epocSig_raw(j,:) = epocSnip;
                
    end
    amp = [];
    auc = [];
    for k = 1:height(epocSig_raw)
        % dF/F
        meanBase = mean(epocSig_raw(k,baseSt:baseEn));
        epoc_dFF(k,1:epocArrayLen) = epocSig_raw(k, 1:epocArrayLen) - meanBase;
        epoc_dFF(k,1:epocArrayLen) = 100*(epoc_dFF(k,1:epocArrayLen) / meanBase);
        % z-Score
        meanBase_dFF = mean(epoc_dFF(k,baseSt:baseEn));
        stdBase_dFF = std(epoc_dFF(k,baseSt:baseEn));
        epocSig_z(k,1:epocArrayLen) = (epoc_dFF(k,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
        % adjusts streams to baseline of zero at -0.5s %
        if epocSig_z(k,idx) < 0
            val = epocSig_z(k,idx);
            diff = 0 - val;
            epocSig_z(k,1:epocArrayLen) = epocSig_z(k,1:epocArrayLen) + abs(diff);
        elseif epocSig_z(k,idx) > 0
            val = epocSig_z(k,idx);
            diff = 0 - val;
            epocSig_z(k,1:epocArrayLen) = epocSig_z(k,1:epocArrayLen) - abs(diff);
        end
    end

    %% Amplitude/AUC %%
    amp = [];
    auc = [];
    for l = 1:height(epocSig_z)
        amp(l,1) = calculateAMP(epocSig_z(l,ampSt:ampEn));
        auc(l,1) = calculateAUC(epocSig_z(l,aucSt:aucEn),ts1);

    end
    cueIdx = find(session_identifiers(:,2) == 0);
    levIdx = find(session_identifiers(:,2) > 0);

    tmp_epocAMP_cue = amp(cueIdx,:);
    tmp_epocAUC_cue = auc(cueIdx,:);
    tmp_epocAMP_lev = amp(levIdx,:);
    tmp_epocAUC_lev = auc(levIdx,:);
    
    epocAMP_cue = [epocAMP_cue;tmp_epocAMP_cue];
    epocAUC_cue = [epocAUC_cue;tmp_epocAUC_cue];
    epocAMP_lev = [epocAMP_lev;tmp_epocAMP_lev];
    epocAUC_lev = [epocAUC_lev;tmp_epocAUC_lev];

    tBYt_cue = epocSig_z(cueIdx,:);
    tBYt_lev = epocSig_z(levIdx,:);
    


    trialByTrial_cue = [trialByTrial_cue;tBYt_cue];
    trialByTrial_lev = [trialByTrial_lev;tBYt_lev];

    id = IDs{i};
    phase = phaseList{i};
    treatment = treatList{i};
    numRows_cue = height(cueIdx);
    numRows_lev = height(levIdx);
    id_cue = repmat({id},numRows_cue,1);
    id_lev = repmat({id},numRows_lev,1);
    id_cue = table(id_cue,'VariableNames',{'ID'});
    id_lev = table(id_lev,'VariableNames',{'ID'});
    idTable_cue = [idTable_cue;id_cue];
    idTable_lev = [idTable_lev;id_lev];
    
    
    tmpphase_cue = repmat({phase},numRows_cue,1);
    tmpphase_cue = table(tmpphase_cue,'VariableNames',{'Phase'});
    phaseTable_cue = [phaseTable_cue;tmpphase_cue];
    tmpphase_lev = repmat({phase},numRows_lev,1);
    tmpphase_lev = table(tmpphase_lev,'VariableNames',{'Phase'});
    phaseTable_lev = [phaseTable_lev;tmpphase_lev];
    
    
    tmptreatment_cue = repmat({treatment},numRows_cue,1);
    tmptreatment_cue = table(tmptreatment_cue,'VariableNames',{'Treatment'});
    treatmentTable_cue = [treatmentTable_cue;tmptreatment_cue];
    tmptreatment_lev = repmat({treatment},numRows_lev,1);
    tmptreatment_lev = table(tmptreatment_lev,'VariableNames',{'Treatment'});
    treatmentTable_lev = [treatmentTable_lev;tmptreatment_lev];

    tempTrialCue = linspace(1,height(cueIdx),height(cueIdx))';
    trialNumCue = [trialNumCue;tempTrialCue];
    tempTrialLev = linspace(1,height(levIdx),height(levIdx))';
    trialNumLev = [trialNumLev;tempTrialLev];

    tempTrialTypeCue = cueArray(:,2);
    trialTypeCue = [trialTypeCue;tempTrialTypeCue];
    tempTrialTypeLev = leverArray(:,2);
    trialTypeLev = [trialTypeLev;tempTrialTypeLev];

end

trialNumCue = array2table(trialNumCue,'VariableNames',{'Trial Num'});
trialNumLev = array2table(trialNumLev,'VariableNames',{'Trial Num'});
trialTypeCue = array2table(trialTypeCue,'VariableNames',{'Trial Type'});
trialTypeLev = array2table(trialTypeLev,'VariableNames',{'Trial Type'});
epocAMP_cue = array2table(epocAMP_cue,'VariableNames',{'Amp'});
epocAMP_lev = array2table(epocAMP_lev,'VariableNames',{'Amp'});
epocAUC_cue = array2table(epocAUC_cue,'VariableNames',{'AUC'});
epocAUC_lev = array2table(epocAUC_lev,'VariableNames',{'AUC'});
totalCols = epocArrayLen;
streamNumbers = 1:totalCols;
streamNames = arrayfun(@(x) ['Stream', num2str(x)], streamNumbers, 'UniformOutput',false);
columnNames = {'Trial Type','Amp','AUC'};
columnNames = [columnNames,streamNames];

trialByTrial_cue = array2table(trialByTrial_cue,"VariableNames",streamNames);
trialByTrial_lev = array2table(trialByTrial_lev,"VariableNames",streamNames);

trialByTrial_cue = horzcat(idTable_cue,phaseTable_cue,treatmentTable_cue,trialNumCue,trialTypeCue,epocAMP_cue,epocAUC_cue,trialByTrial_cue);
trialByTrial_lev = horzcat(idTable_lev,phaseTable_lev,treatmentTable_lev,trialNumLev,trialTypeLev,epocAMP_lev,epocAUC_lev,trialByTrial_lev);

tBYt_cue_acq2 = trialByTrial_cue(strcmpi(trialByTrial_cue.Phase,'Acq2'),:);
tBYt_cue_rev1 = trialByTrial_cue(strcmpi(trialByTrial_cue.Phase,'Rev1'),:);

tBYt_lev_acq2 = trialByTrial_lev(strcmpi(trialByTrial_lev.Phase,'Acq2'),:);
tBYt_lev_rev1 = trialByTrial_lev(strcmpi(trialByTrial_lev.Phase,'Rev1'),:);