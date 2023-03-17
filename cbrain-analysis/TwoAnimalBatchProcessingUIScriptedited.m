%% Fiber Photometry Epoch Averaging Example
%
% <html>
% This example goes through fiber photometry analysis using techniques <br>
% such as data smoothing, bleach detrending, and z-score analysis. <br>
% The epoch averaging was done using TDTfilter. <br><br>
% Author Contributions: <br>
% TDT, David Root, and the Morales Lab contributed to the writing and/or conceptualization of the code. <br>
% The signal processing pipeline was inspired by the workflow developed by <a href="https://doi.org/10.1016/j.celrep.2017.10.066">David Barker et al. (2017)</a> for the Morales Lab. <br>
% The data used in the example were provided by David Root. <br><br>
% Author Information: <br>
% David H. Root <br>
% Assistant Professor <br>
% Department of Psychology & Neuroscience <br>
% University of Colorado, Boulder <br>
% Lab Website: <a href="https://www.root-lab.org">https://www.root-lab.org</a> <br>
% david.root@colorado.edu <br><br>clea
% About the authors: <br>
% The Root lab and Morales lab investigate the neurobiology of reward, aversion, addiction, and depression. <br>
% <br> TDT edits all user submissions in coordination with the contributing
% author(s) prior to publishing.
% </html>

% % Edited / adapted for multi-bin analysis by AK - 2/2020, automated batch processing by SR - 3/2021

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab path.
close all;
clear all; 
clc;

%% User Input Section

% IMPORTANT!!!! - Folder containing the individual recording files must be named
% with the following convention: 

%                BoxA_AnimalNameA_BoxB_AnimalNameB_Treatment_Date

% Example -      Box5_JM3_Box6_JF7_Vehicle_12-25

% Underscores MUST be in appropriate places, variables are made based on
% placement of underscores

% Note - All recordings in a single batch must be from the same rig. The
% current script cannot distinguish between TTLs of different rigs.

% Location of parent directory containing TDTbin2mat folder
%addpath('/Users/brandon/Desktop/DA_PRL_DATA/FP_SCRIPTS');
addpath('C:\Users\brand\iCloudDrive\Desktop\DA_PRL_DATA\FP_SCRIPTS');

% Default File Save Location
%SavePATH = fullfile('/Users/brandon/Desktop');
SavePATH = fullfile('C:\Users\brand\iCloudDrive\Desktop')

% If using new TTLs, see section near line 100 of this script to input new TTL values 

%% UI to generate checkbox for operant box selection
cbxfig = uifigure('Name','Check which operant boxes to analyze','Position',[600 500 400 100]);
cbx1 = uicheckbox(cbxfig,'Text','Box A','Position',[30 50 91 15]);
cbx2 = uicheckbox(cbxfig,'Text','Box B','Position',[130 50 91 15]);
subbtn = uibutton(cbxfig,'Position',[250 50 91 20],'Text','OK',...
              'ButtonPushedFcn','uiresume(gcbf)');
disp('Please select box(es) to analyze');
uiwait(cbxfig);

if cbx1.Value == 1
    LoopValue1 = 1;
elseif cbx1.Value == 0
    LoopValue1 = 2;
end

if cbx2.Value == 1
    LoopValue2 = 2;
elseif cbx2.Value == 0
    LoopValue2 = 1;
end

if LoopValue1 == 2 && LoopValue2 == 1
    disp('Error: Please check a box to analyze');
end
close(cbxfig);

%% UI to pick which loops to run: Old Rig (3) or New Rig (4) and associated TTLs

prompt = {'Enter Rig Number'};
dlgtitle = 'User Input';
dims = [1 35];
definput = {'4'};
UserInput = inputdlg(prompt,dlgtitle,dims,definput);
RigNumber = str2num(UserInput{1});

% Name of TTLs
TTL1 = 'Reward_Dispensed';
TTL2 = 'Lever_Press';
TTL3 = 'Left_Cue_Light';
TTL4 = 'Interval_End';

% If using a different rig, please add correct TTLs in same convention as
% below
if RigNumber == 2
    %                    {Box A, Box B} 
    T.Reward_Dispensed = {'PC0','PC1'};
    T.House_Light = {'PC2','PC3'};
    T.Right_Nosepoke = {'PC6','PC7'};
    T.Left_Nosepoke = {'PC4','PC5'};
elseif RigNumber == 3
    %                    {Box A, Box B} 
    T.Reward_Dispensed = {'RwA','RwB'};
    T.Lever_Press = {'RLA','RLB'};
    T.Left_Cue_Light = {'ItA','ItB'};
    T.Interval_End = {'leA','leB'};
elseif RigNumber == 4
    %                    {Box A, Box B} 
    T.Reward_Dispensed = {'aRw','bRw'};
    T.Lever_Press = {'RLA','RLB'};
    T.Left_Cue_Light = {'ItA','ItB'};
    T.Interval_End = {'IeA','IeB'};
elseif RigNumber == 10 || 11
%     Stationary rigs in mouse ephys room
    T.Reward_Dispensed = {'RwA','RwB'};
    T.Lever_Press = {'RLA','RLB'};
    T.Left_CueLight = {'ItA','ItB'};
    T.Interval_End = {'leA','leB'};
else 
    disp('Rig Number Invalid');
end

% UI to generate drop-down list for TTL selection
TTLList = fieldnames(T);
[indx,tf] = listdlg('PromptString',{'Select TTLs to Plot, Ctrl+Click for multiple'},'ListString',TTLList);
TTLNum = numel(TTLList(indx));
RunNum = 1;

%% Select folder copied from Tanks on the rig
% Parent folder where recordings are kept (Copied from Tanks on the rig)
DATAPATH = uigetdir(SavePATH,'Select parent folder containing recordings from same rig on your computer');
root = DATAPATH;
FolderList = dir(root);

%% Beginning of loop to run different files in parent directory

for k = 3 :  length(FolderList)
animID = FolderList(k).name; % Name of folder containing your specific recording

% Folder name used to create naming variables
BrokenAnimID = strsplit(animID,'_');
animName_A = char(BrokenAnimID(1,2));
animName_B = char(BrokenAnimID(1,6));
% Treatment used in recording
TreatmentTypeA = char(BrokenAnimID(1,3));
TreatmentTypeB = char(BrokenAnimID(1,7));
% Date of recording
RecDateA = char(BrokenAnimID(1,4));
RecDateB = char(BrokenAnimID(1,8));

%% Creating directories based on user input

BLOCKPATH = fullfile(DATAPATH,animID);
animDir_A = fullfile(SavePATH,animName_A,RecDateA);
SavePATH_A = fullfile(animDir_A);
animDir_B = fullfile(SavePATH,animName_B,RecDateB);
SavePATH_B = fullfile(animDir_B);



% Beginning of loop to run both boxes
for LoopNumber = LoopValue1:LoopValue2

if LoopNumber == 1
    SigSave = SavePATH_A;
    animName = animName_A;
    TreatmentType = TreatmentTypeA;
    RecDate = RecDateA;
    mkdir(animDir_A);
elseif LoopNumber == 2
    SigSave = SavePATH_B;
    animName = animName_B;
    TreatmentType = TreatmentTypeB;
    RecDate = RecDateB;
    mkdir(animDir_B);
end   

for TTLLoop = 1:TTLNum
    
    if indx(TTLLoop) == 1
        TTL_A = T.(TTL1){1};
        TTL_B = T.(TTL1){2};
        TTLTitle = 'Reward Dispensed';
    elseif indx(TTLLoop) == 2
        TTL_A = T.(TTL2){1};
        TTL_B = T.(TTL2){2};
        TTLTitle = 'Lever Press';
    elseif indx(TTLLoop) == 3
        TTL_A = T.(TTL3){1};
        TTL_B = T.(TTL3){2};
        TTLTitle = 'Left Cue Light';
    elseif indx(TTLLoop) == 4
        TTL_A = T.(TTL4){1};
        TTL_B = T.(TTL4){2};
        TTLTitle = 'Interval End';
    end
    
    %% Setup the variables for the data you want to extract
% We will extract two different stream stores surrounding the TTL epoch
% event. We are interested in a specific event code for reaching PR.    
if RigNumber == 2
    
    if LoopNumber == 1
        REF_EPOC = append(TTL_A,'/'); % TTLs - BOX A
        STREAM_STORE1 = 'x465A'; % name of the 405 store - BOX A
        STREAM_STORE2 = 'x405A'; % name of the 470 store - BOX A
        TTLName = TTL_A;
    elseif LoopNumber == 2
        REF_EPOC = append(TTL_B,'/'); % TTLs - BOX B
        STREAM_STORE1 = 'x465C'; % name of the 405 store - BOX B
        STREAM_STORE2 = 'x405C'; % name of the 470 store - BOX B
        TTLName = TTL_B;
    else
        disp('Error in Box-specific TTL determination');
    end
    
elseif RigNumber == 3
    
    if LoopNumber == 1
        REF_EPOC = append(TTL_A,'/'); % TTLs - BOX A
        STREAM_STORE1 = 'x405A'; % name of the 405 store - BOX A
        STREAM_STORE2 = 'x470A'; % name of the 470 store - BOX A
        TTLName = TTL_A;
    elseif LoopNumber == 2
        REF_EPOC = append(TTL_B,'/'); % TTLs - BOX B
        STREAM_STORE1 = 'V2_B'; % name of the 405 store - BOX B
        STREAM_STORE2 = 'B2_B'; % name of the 470 store - BOX B
        TTLName = TTL_B;
    else
        disp('Error in Box-specific TTL determination');
    end

elseif RigNumber == 4
    
    % MUST INCLUDE 'x' before stream name!!!!!!
    if LoopNumber == 1
        REF_EPOC = append(TTL_A,'/'); % REWARD TTLs - BOX A
        STREAM_STORE1 = 'x405A'; % name of the 405 store - BOX A
        STREAM_STORE2 = 'x465A'; % name of the 470 store - BOX A
        TTLName = TTL_A;
    elseif LoopNumber == 2
        REF_EPOC = append(TTL_B,'/'); % REWARD TTLs - BOX B
        STREAM_STORE1 = 'x405C'; % name of the 405 store - BOX B
        STREAM_STORE2 = 'x465C'; % name of the 470 store - BOX B
        TTLName = TTL_B;
    else 
        disp('Error in Box-specific TTL determination');
    end
    
elseif RigNumber == 10 || 11
%     Stationary rigs in mouse ephys room
    if LoopNumber == 1
        REF_EPOC = append(TTL_A,'/'); % TTLs - BOX A
        STREAM_STORE1 = 'x405A'; % name of the 405 store - BOX A
        STREAM_STORE2 = 'x465A'; % name of the 470 store - BOX A
        TTLName = TTL_A;
    elseif LoopNumber == 2
        REF_EPOC = append(TTL_B,'/'); % TTLs - BOX B
        STREAM_STORE1 = 'x405C'; % name of the 405 store - BOX B
        STREAM_STORE2 = 'x465C'; % name of the 470 store - BOX B
        TTLName = TTL_B;
    else
        disp('Error in Box-specific TTL determination');
    end
else 
    disp('Error in Rig-specific TTL Determination');
end

disp(append(append(animName,' ',TreatmentType,' ',TTLName,' ',RecDate),' Start')); % Display plotting status in command window

TRANGE = [ -2 37]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-2 -1.5]; % baseline period within our window
ARTIFACT = Inf; % optionally set an artifact rejection level

% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'});

%% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event. Use the 'VALUES' parameter to specify allowed values of
% the REF_EPOC to extract.  For stream events, the chunks of data are
% stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);

%%
% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);

%%
% Applying a time filter to a uniformly sampled signal means that the
% length of each segment could vary by one sample.  Let's find the minimum
% length so we can trim the excess off before calculating the mean.
minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);

allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);

% downsample 10x and average 465 signal
allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F470 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F470(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength2 = size(F470,2);

% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal2 = mean(F470);
stdSignal2 = std(double(F470))/sqrt(size(F470,1));
dcSignal2 = mean(meanSignal2);

% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;


% bls = polyfit(F405(1:end), F470(1:end), 1);
bls = polyfit(F470(1:end), F405(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
dFF = F470 - Y_fit_all;

% z-scores for dF/F matrix
zall = zeros(size(dFF));
for i = 1:size(dFF,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(dFF(i,ind)); % baseline period mean (-10sec to -6sec)
    zsd = std(dFF(i,ind)); % baseline period stdev
    zall(i,:)=(dFF(i,:) - zb)/zsd; % Z score per bin
end

zerror = std(zall)/sqrt(size(zall,1));

% Fill band values for second subplot. Doing here to scale onset bar
% correctly
XX =[ts2, fliplr(ts2)];
YY =[mean(zall)-zerror, fliplr(mean(zall)+zerror)];

% PLOT All trials    
figure;
% subplot(4,1,3)
plot(ts2, mean(zall), 'color',[.3 .75, .93], 'LineWidth', 3); hold on;
%Blue color for WT
plot(ts2, mean(zall), 'color', [0 0.4470 0.7410], 'LineWidth', 3); hold on;

%Red color for Het
% plot(ts2, mean(zall), 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;

% line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)
line([0 0], [floor(min(YY)), ceil(max(YY))], 'Color', [0 0 0], 'LineWidth', 1.3, 'LineStyle', '--')

% Blue shading for WT
h = fill(XX, YY, 'b');
set(h, 'facealpha',.25,'edgecolor','none')

% %Red shading for Het
% h = fill(XX, YY, 'r');
% set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
ylim([floor(min(YY)), ceil(max(YY))]);
xlabel('Time [s]','FontSize',12)
ylabel('Z-score', 'FontSize', 12)
%title(sprintf('Reward Response, All Trials', size(zall,1)));
title(append(animName,' Cue - ',TreatmentType,', ',TTLName,' - ',RecDate));
% title(append(animName,' DA Response - ',TTLTitle));
%c2 = colorbar;

fig1Name = append('CuelightResponse',animName,'_',TreatmentType,'_',TTLName,'_',RecDate);
saveas(gcf,fullfile(SigSave,fig1Name),'bmp');
saveas(gcf,fullfile(SigSave,fig1Name),'png');

        %zall = zeros(size(Y_dF_all));

figure;
imagesc(ts2, 1, zall);
colormap('jet'); 
TrialNum = sprintf('%d',size(zall,1));
title(append(animName,' Z-Score Heat Map ',TrialNum,' Trials - ',TreatmentType,'Cue',', ',TTLName,' - ',RecDate));
% title(sprintf('Z-Score Heat Map, %d Trials', size(zall,1) ));
colorbar;
% title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 12);
xlabel('Time [s]','FontSize',12)
set(gcf,'renderer','painters')

fig2Name = append('HeatMap',animName,'_',TreatmentType,'Cuelight','_',TTLName,'_',RecDate);
 saveas(gcf,fullfile(SigSave,fig2Name),'bmp');
 saveas(gcf,fullfile(SigSave,fig2Name),'png');

%% TRACES & ZScores - in n BINS
% where n <= total events / 2
% aka the limit to smallest bin size is 2 events
% Split into third bins

bin = 1;
trialBin_lim = floor(size(F470,1)/bin);
trialBin_last = trialBin_lim + mod(size(F470,1), bin);

tsZero = find(abs(ts1-0.001) < 0.01);
tsUpper = tsZero(1)+205;
tsLower = tsZero(1)-100;
TRANGE2 = [tsZero(1) : tsUpper]; 

tsZero = find(abs(ts1-0.001) < 0.01);
tsUpper1 = tsZero(1)+905;
tsLower1 = tsZero(1)+300;
TRANGE3 = [tsLower1(1) : tsUpper1]; 


rowIndex = 1;
totBKG = [];
totSIG = [];
tot_zall = [];



for b = 1:bin

        
    % for last bin - accounting for Total Trial # (ODD)
    if  b == bin
        trialBin_lim = trialBin_last;
    end
    
    c = rowIndex + trialBin_lim - 1;
    % z-scores for dF/F matrix, baseline as defined by epoch @ start [-5 0.5]
    temp_zall = zall(rowIndex:c, :);
    
    temp_zerror = std(temp_zall)/sqrt(size(temp_zall,1));
    SIGz(b,:) = {temp_zall, temp_zerror};

%     
%     % Fitting 405 channel onto 465 channel to detrend signal bleaching
%     % Scale and fit data
%     % Algorithm sourced from Tom Davidson's Github:
%     % https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
%     
% %     temp_bls = polyfit(temp_BKG(1:end), temp_SIG(1:end), 1);
% temp_bls = polyfit(temp_SIG(1:end),temp_BKG(1:end), 1);
%     temp_Y_fit_all = temp_bls(1) .* temp_BKG + temp_bls(2);
%     temp_Y_dF_all = temp_SIG - temp_Y_fit_all;
%     
%     % dF/F
%     temp_dFF = 100*(temp_Y_dF_all)./temp_Y_fit_all;
    
%    temp_dFF = temp_Y_dF_all;


% Changes bin within for loop
    temp_zall = zall(rowIndex: c , :);
    temp_zerror = std(temp_zall)/sqrt(size(temp_zall,1));    
    SIGz(b,:) = {temp_zall, temp_zerror};    
    rowIndex = (c + 1);
    
%% SAVE SIG matrix containing z-scores to excel files in created directory
%  Comment out which sheets you don't need

cd(SigSave)
filename = append('Z-Scores_',animName,'_',TreatmentType,'_','Cue',RecDate,'.xlsx'); 
 filename2 = fullfile(SavePATH,append(TreatmentType,'_','PeakDAData.xlsx'));

% Matrices containing all z-scores, binned z-scores, and mean z-scores 
zallmatrix = [ts1;zall];
writematrix(zallmatrix,filename,'Sheet','All Z-scores','Range','A2');
% 
 sigz = cell2mat(SIGz(b,1)); 
 sigmatrix = [ts1;sigz];

 writematrix(sigmatrix,fullfile(SavePATH,filename),'Sheet',append('Bin ',string(b),' Z-scores'),'Range','A2');

 %if b == 3
    
cd(SavePATH)
filename3 = append('DA_Data_B1',animName,TreatmentType,RecDate,'Cue','.xlsx'); 
meanzall = mean(zall);
mbzallmatrix = [meanzall];
PlusCell = ['A',num2str(RunNum + 1)];
writematrix(ts1,fullfile(SavePATH,filename3),'Sheet',append('Bin ',string(b),' Mean z-scores'),'Range','A1');
writematrix(mbzallmatrix,fullfile(SavePATH,filename3),'Sheet',append('Bin ',string(b),' Mean z-scores'),'Range',PlusCell);

% Maximum Signal Peaks within 2 seconds of selected TTL onset
cd(SavePATH)
filename4 = 'MaxPeaksMasterRig1.5.5.xlsx';
MaxPeakCue = max(meanzall(TRANGE2));
ExportMaxPeaks = {animName,RecDate,TreatmentType,MaxPeakCue,'Cue'};
EMPcell = ['A',num2str(RunNum + 1)];
EMPtitles = {'Subject','Date','Treatment','Max DA Peak'};
writecell(EMPtitles,fullfile(SavePATH,filename4),'Sheet','Max Peak Cue','Range','A1');
writecell(ExportMaxPeaks,fullfile(SavePATH,filename4),'Sheet','Max Peak Cue','Range',EMPcell); 

cd(SavePATH)
MaxPeakReward = max(meanzall(TRANGE3));
ExportMaxPeaks = {animName,RecDate,TreatmentType,MaxPeakReward,'Reward'};
EMPcell = ['A',num2str(RunNum + 1)];
EMPtitles = {'Subject','Date','Treatment','Max DA Peak'};
writecell(EMPtitles,fullfile(SavePATH,filename4),'Sheet','Max Peak Reward','Range','A1');
writecell(ExportMaxPeaks,fullfile(SavePATH,filename4),'Sheet','Max Peak Reward','Range',EMPcell); 

 %end

end % loop to split into bins
disp(append(append(animName,' ',TreatmentType,' ',TTLName,' ',RecDate),' Done')); % Display plotting status in command window
end % loop to run all TTLs
RunNum = RunNum + 1;
end % loop to run both boxes
end % loop to run all recordings in parent folder



