function [data] = wd_epocs(filename, tsDir, data)

tsLeftPath = strcat(tsDir,'/',filename,'_','leftCup','.csv');
tsRightPath = strcat(tsDir,'/',filename,'_','rightCup','.csv');
leftCup = readtable(tsLeftPath);
rightCup = readtable(tsRightPath);


data.epocs.lApp.onset = table2array(leftCup(:,"LeftStart"));
data.epocs.lApp.offset = table2array(leftCup(:,"LeftStop"));
data.epocs.rApp.onset = table2array(rightCup(:,"RightStart"));
data.epocs.rApp.offset = table2array(rightCup(:,"RightStop"));

