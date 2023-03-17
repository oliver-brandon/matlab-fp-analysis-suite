clear all;
load("Z:\DA_PRL\PRL_Mats_Split\Rev1\124F_Rev1_JZL8.mat")

A1 = data.epocs.cRewA.onset;
A2 = data.epocs.cNoRewA.onset;
A3 = data.epocs.iRewA.onset;
A4 = data.epocs.iNoRewA.onset;

X1 = ones(height(A1),1);
X2 = ones(height(A2),1)*2;
X3 = ones(height(A3),1)*3;
X4 = ones(height(A4),1)*4;
% Append first column of each array to new array
combined_array1 = [A1(:,1); A2(:,1); A3(:,1); A4(:,1)];
combined_array2 = [X1(:,1); X2(:,1); X3(:,1); X4(:,1)];

final_array = [combined_array1 combined_array2];
% Sort the first column of newArray in ascending order
final_array = sortrows(final_array, 1);
% Find the rows in A where the first column is not zero
idx = final_array(:,1) ~= 0;
% Use logical indexing to remove rows with zero values
final_array = final_array(idx,:);

session_ts = final_array(:,1);
trial_type = final_array(:,2);

cue_lever_ts(:,1) = data.epocs.St1_.onset;
session_ts(2,:) = [];
cue_lever_ts(:,2) = session_ts;
lever_ts = cue_lever_ts(:,2) - cue_lever_ts(:,1);