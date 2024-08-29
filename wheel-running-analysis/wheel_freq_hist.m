clear
myDir = uigetdir('/Users/brandon/personal-drive/hot-wheels/wheel-running-mats/unlocked'); %gets directory%
if myDir == 0
    disp("Select a folder of files to start")
    return
end

tic
myFiles = dir(myDir); %gets all tanks in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._'}));
myFiles = myFiles(endsWith({myFiles.name}, {'.mat'}));
numFiles = length(myFiles);
IDs = cell(size(numFiles,1));
dayList = cell(size(numFiles,1));
counts = [];

for i = 1:numFiles
    durationsTemp = [];
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);

    brokenID = strsplit(name,'_');
    IDs{i} = cellstr(strtrim(brokenID{1}));
    dayList{i} = cellstr(strtrim(brokenID{3}));

    load(filename)
    fprintf('Analyzing %s (%d of %d)\n', name, i, numFiles)
    onWheel = data.epocs.onWheel.onset;
    runStart = data.epocs.runStart.onset;

    for j = 1:length(runStart)
        rS = runStart(j);
        idx = find(onWheel<rS,1,'last');
        oW = onWheel(idx);
        dur = rS - oW;
        durationsTemp = [durationsTemp;dur];
    end
    var1 = sum(durationsTemp > 0 & durationsTemp < 3);
    var2 = sum(durationsTemp > 3 & durationsTemp < 6);
    var3 = sum(durationsTemp > 6 & durationsTemp < 9);
    var4 = sum(durationsTemp > 9 & durationsTemp < 12);
    var5 = sum(durationsTemp > 12 & durationsTemp < 15);
    var6 = sum(durationsTemp > 15 & durationsTemp < 20);
    var7 = sum(durationsTemp > 20);

    counts(i,1:7) = [var1,var2,var3,var4,var5,var6,var7];

end
IDs = cell2table(IDs','VariableNames',{'ID'});
dayList = cell2table(dayList','VariableNames',{'Day'});
counts = array2table(counts,'VariableNames',{'0-3s','3-6s','6-9s','9-12s','12-15s','15-20s','20s+'});

dataTable = horzcat(IDs,dayList,counts);
dataTable = sortrows(dataTable,{'ID','Day'},{'ascend','ascend'});
dataSum(:,1:7) = [sum(dataTable.("0-3s")),sum(dataTable.("3-6s")),sum(dataTable.("6-9s")),...
    sum(dataTable.("9-12s")),sum(dataTable.("12-15s")),sum(dataTable.("15-20s")),sum(dataTable.("20s+"))];
toc
disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)

NERD_STATS(toc,numFiles);