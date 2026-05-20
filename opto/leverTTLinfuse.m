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
    filename = fullfile(myDir, myFiles(i).name);
    [~, name, ~] = fileparts(filename);
    fprintf('Infusing %s (%d of %d)\n', name, i, numFiles)

    S = load(filename);          % safer than bare load()
    data = S.data;

    % Skip if already infused
    if isfield(data.epocs, 'levers') && isfield(data.epocs.levers, 'onset')
        disp("already infused")
        continue
    end

    hasCL = isfield(data.epocs, 'CL1_');
    hasIL = isfield(data.epocs, 'IL1_');

    % Determine lever onsets
    if hasCL && hasIL
        levers = sort([data.epocs.CL1_.onset; data.epocs.IL1_.onset]);
    elseif hasCL
        levers = data.epocs.CL1_.onset;
    elseif hasIL
        levers = data.epocs.IL1_.onset;
    else
        disp('missing epocs')
        continue
    end

    % Write infused epoc
    data.epocs.levers.onset  = levers;
    data.epocs.levers.offset = levers + 1;
    data.epocs.levers.name   = 'levers';
    data.epocs.levers.data   = ones(numel(levers), 1);

    save(filename, 'data');
end
disp("Successfully analyzed .mat files")
fprintf("Files analyzed: %d\n", numFiles)
toc
NERD_STATS(toc,numFiles);