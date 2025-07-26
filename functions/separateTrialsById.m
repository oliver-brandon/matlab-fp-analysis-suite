function [group1, group2, meansGroup1, meansGroup2] = separateTrialsById(inputTable)
    % Extract unique ids
    uniqueIds = unique(inputTable.id);
    
    % Initialize output cell arrays to hold data for each id
    group1 = cell(length(uniqueIds), 2); % For trialNum 1-14
    group2 = cell(length(uniqueIds), 2); % For trialNum 15-29
    
    % Initialize arrays for storing column means
    meansGroup1 = cell(length(uniqueIds), 2);
    meansGroup2 = cell(length(uniqueIds), 2);
    
    % Iterate over each unique id
    for i = 1:length(uniqueIds)
        % Extract rows corresponding to the current id
        currentId = uniqueIds(i);
        idRows = inputTable(inputTable.id == currentId, :);
        
        % Separate rows by trialOutcome
        rowsOutcome1 = idRows(idRows.trialOutcome == 1, :);
        rowsOutcome0 = idRows(idRows.trialOutcome == 0, :);
        
        % Further separate by trialNum 1-14 and 15-29 for both trialOutcome
        % For trialNum 1-14
        group1_1to14 = rowsOutcome1(rowsOutcome1.trialNum >= 1 & rowsOutcome1.trialNum <= 14, 4:end);
        group1_0to14 = rowsOutcome0(rowsOutcome0.trialNum >= 1 & rowsOutcome0.trialNum <= 14, 4:end);
        
        % For trialNum 15-29
        group2_1to29 = rowsOutcome1(rowsOutcome1.trialNum >= 15 & rowsOutcome1.trialNum <= 29, 4:end);
        group2_0to29 = rowsOutcome0(rowsOutcome0.trialNum >= 15 & rowsOutcome0.trialNum <= 29, 4:end);
        
        % Store in output cell arrays
        group1{i, 1} = group1_1to14; % Outcome 1, trialNum 1-14
        group1{i, 2} = group1_0to14; % Outcome 0, trialNum 1-14
        group2{i, 1} = group2_1to29; % Outcome 1, trialNum 15-29
        group2{i, 2} = group2_0to29; % Outcome 0, trialNum 15-29
        
        % Calculate column means for each group
        if ~isempty(group1_1to14)
            meansGroup1{i, 1} = mean(group1_1to14, 1); % Outcome 1, trialNum 1-14
        else
            meansGroup1{i, 1} = []; % Empty if no data
        end
        
        if ~isempty(group1_0to14)
            meansGroup1{i, 2} = mean(group1_0to14, 1); % Outcome 0, trialNum 1-14
        else
            meansGroup1{i, 2} = []; % Empty if no data
        end
        
        if ~isempty(group2_1to29)
            meansGroup2{i, 1} = mean(group2_1to29, 1); % Outcome 1, trialNum 15-29
        else
            meansGroup2{i, 1} = [];
        end
        
        if ~isempty(group2_0to29)
            meansGroup2{i, 2} = mean(group2_0to29, 1); % Outcome 0, trialNum 15-29
        else
            meansGroup2{i, 2} = [];
        end
    end
end