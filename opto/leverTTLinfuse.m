clear

myDir = uigetdir(...
    '/Users/brandon/personal-drive/optomouse-prime/opto-reversal','Choose the .mat files you want to analyze.'...
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

for i = 1:numFiles
    filename = fullfile(myDir,myFiles(i).name);
    [~,name,~] = fileparts(filename);
    fprintf('Infusing %s (%d of %d)\n', name, i, numFiles)
    load(filename)
    if ~isfield(data.epocs, 'CL1_')
        levers = data.epocs.IL1_.onset;
    elseif ~isfield(data.epocs, 'IL1_')
        levers = data.epocs.CL1_.onset;
    elseif isfield(data.epocs, 'CL1_') && isfield(data.epocs, 'IL1_')
        correct = data.epocs.CL1_.onset;
        incorrect = data.epocs.IL1_.onset;
        levers = sort([correct;incorrect]);
    else
        disp('missing epocs')
        break
    end
    data.epocs.levers.onset = levers;
    data.epocs.levers.offset = levers + 1;
    data.epocs.levers.name = 'levers';
    data.epocs.levers.data = ones(height(levers),1);
    save(filename, 'data');
end
disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)
toc
NERD_STATS(toc,numFiles);