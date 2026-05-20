function [data] = wd_epocs(filename, tsDir, data)

tsLeftPath = fullfile(tsDir, strcat(filename,'_','leftCup','.csv'));
tsRightPath = fullfile(tsDir, strcat(filename,'_','rightCup','.csv'));
if ~isfile(tsLeftPath) || ~isfile(tsRightPath)
    error('wd_epocs:MissingCsv', 'Missing WD timestamp CSV for %s.', filename);
end
leftCup = readtable(tsLeftPath);
rightCup = readtable(tsRightPath);

requiredLeft = {'LeftStart','LeftStop'};
requiredRight = {'RightStart','RightStop'};
if ~all(ismember(requiredLeft, leftCup.Properties.VariableNames)) || ...
        ~all(ismember(requiredRight, rightCup.Properties.VariableNames))
    error('wd_epocs:MissingColumns', 'WD timestamp CSV files are missing required columns.');
end

data.epocs.lApp.onset = table2array(leftCup(:,"LeftStart"));
data.epocs.lApp.offset = table2array(leftCup(:,"LeftStop"));
data.epocs.rApp.onset = table2array(rightCup(:,"RightStart"));
data.epocs.rApp.offset = table2array(rightCup(:,"RightStop"));

