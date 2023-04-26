warning off
myDir = uigetdir(pwd); %gets directory%
if myDir == 0
    disp("Select a folder to start")
    return
end
tic
myFiles = dir(myDir); %gets all files in directory%
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','._','~'}));
myFiles = myFiles(endsWith({myFiles.name},'.xlsx'));
numFiles = length(myFiles);
for i = 1:numFiles
    % Load data from Excel file
    filename = fullfile(myDir,myFiles(i).name);
    data = readtable(filename);
    
    % Get the number of columns in the data
    numCols = size(data, 2);
    
    % Create a new matrix to hold the summed and divided values
    newData = zeros(1, numCols);
    
    % Loop through each column
    for i = 2:numCols
        % Get the column data
        colData = data{:, i};
        
        % Find all the numbers that end with ".9"
        idx = round(mod(colData, 1), 1) == 0.9;
        
        % Sum the numbers in front of ".9" and divide by 100
        if any(idx)
            newData(i) = sum(floor(colData(idx))) / sum(idx) / 100;
        end
    end
    
    % Write the new data to a new sheet in the excel file
    writetable(array2table(newData), filename, 'Sheet', 'AvgValues');
end
toc
NERD_STATS(toc,numFiles);