clear
myDir = uigetdir('','Choose the mat file(s) you want to infuse.'); %gets directory%
disp('Choose a folder containing one or more mat files that you wish to infuse with wheel epochs.')
if myDir == 0
    disp("Select a directory of mat files to start")
    return
end
savDir = uigetdir('','Choose where you want to save the infused mat file(s).'); %gets directory%
disp('Choose a folder to save the infused mat files to.')
if savDir == 0
    disp("Select a valid save directory")
    return
end
tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name},'.mat'));
numFiles = length(myFiles);
missing_notes = {};
for i = 1:numFiles
    FILEPATH = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(FILEPATH);
    load(FILEPATH);
    if isfield(data.epocs, 'onWheel') || isfield(data.epocs, 'offWheel')...
            || isfield(data.epocs, 'runStart') || isfield(data.epocs, 'runStop')
        fprintf('%s already infused...skipping (%d of %d)\n',name,i,numFiles)
        continue
    else
        fprintf('Infusing %s...(%d of %d)\n',name,i,numFiles)
    end
    %% makes new wheel data epocs %%
    %% epocs created in TDT OpenScope save in the notes of Cam1 %% 
    
    %combine index with timestamp data from Cam1 notes%
    
    if ~isfield(data.epocs.Cam1,'notes')
        missing_notes = [missing_notes; {name}];
        fprintf('%s missing notes...skipping (%d of %d)\n',name,i,numFiles)
        continue
    else
        ind = double(data.epocs.Cam1.notes.index);
        ts = data.epocs.Cam1.notes.ts;
        var1 = [ind ts];
    end
    
    %separate by index to make separate epocs%
    onWheel = var1(ismember(var1(:,1),[1]),:);
    offWheel = var1(ismember(var1(:,1),[2]),:);
    runStart = var1(ismember(var1(:,1),[3]),:);
    runStop = var1(ismember(var1(:,1),[4]),:);
    
    %extract time stamps%
    onWheelTs = onWheel(:,2);
    offWheelTs = offWheel(:,2);
    runStartTs = runStart(:,2);
    runStopTs = runStop(:,2);
    %extract indicies%
    onWheelInd = onWheel(:,1);
    offWheelInd = offWheel(:,1);
    runStartInd = runStart(:,1);
    runStopInd = runStop(:,1);
    %make a new epoc structure based on Cam1 notes extracted data%
    %onWheel%
    data.epocs.onWheel.onset = onWheelTs;
    data.epocs.onWheel.offset = onWheelTs + 0.01;
    data.epocs.onWheel.name = 'onWheel';
    data.epocs.onWheel.data = onWheelInd;
    %offWheel%
    data.epocs.offWheel.onset = offWheelTs;
    data.epocs.offWheel.offset = offWheelTs + 0.01;
    data.epocs.offWheel.name = 'offWheel';
    data.epocs.offWheel.data = offWheelInd;
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

    save(FILEPATH,"data")
end
toc
disp("Successfully infused epocs into mat files")
fprintf("Files infused: %d\n",numFiles)
fprintf("Save location: %s\n",savDir)

NERD_STATS(toc,numFiles);