

cRew = data.epocs.cRewA.onset;
cNoRew = data.epocs.cNoRewA.onset;
iRew = data.epocs.iRewA.onset;
iNoRew = data.epocs.iNoRewA.onset;

cRew(1:end,2) = 1;
cRew(1:end,3) = 1;
cNoRew(1:end,2) = 1;
cNoRew(1:end,3) = 0;
iRew(1:end,2) = 0;
iRew(1:end,3) = 1;
iNoRew(1:end,2) = 0;
iNoRew(1:end,3) = 0;

% Combine the onset data into a single matrix for further analysis
allTrials = [cRew; cNoRew; iRew; iNoRew];
allTrials = sortrows(allTrials,1,"ascend");

% Eliminate rows that have 0 in column 1
allTrials(allTrials(:,1) == 0, :) = [];