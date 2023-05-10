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
    fprintf('Loading file %d of %d.\n',i,numFiles)
    
    % Load data from Excel file
    filename = fullfile(myDir,myFiles(i).name);
    data = readtable(filename);
    
    % Get the number of columns in the data
    numCols = size(data, 2);
    
    % Create a new matrix to hold the summed and divided values
    avgHE = zeros(1, numCols);
    heTS = zeros(numCols,1);

    disp('Calculating...')
    
    % Loop through each column
    for j = 2:numCols
        f = 1;
        % Get the column data
        colData = data{:, j};
        
        % Find all the numbers that end with ".9"
        idx = round(mod(colData, 1), 1) == 0.9;
        
        % Sum the numbers in front of ".9" and divide by 100
        if any(idx)
            avgHE(j) = sum(floor(colData(idx))) / sum(idx) / 100;
            
        end
        for k = 1:height(idx)
            
            if idx(k) == 1
                heTS(f,j) = floor(colData(k)) / 100;
                f = f + 1;
            end
        end
    end
    
    disp('Writing...')
    % Write the new data to a new sheet in the excel file
    writetable(array2table(avgHE), filename, 'Sheet', 'Avg HE Per Session');
    writetable(array2table(heTS), filename, 'Sheet', 'HE Timestamps');
end
toc
NERD_STATS(toc,numFiles);