function errorProbTS = errorProbExtract( ...
    trial_type, ...
    session_ts, ...
    errorType, ...
    lever...
    )
% errorType: 1 = winStay, 2 = winShift, 3 = loseStay, 4 = loseShift
% lever: errors based on correct lever = 1, incorrect lever = 2
% Initialize index array
index_array = [];
if lever == 1
    reward = 1;
    noreward = 2;
    shift1 = 3;
    shift2 = 4;
elseif lever == 2
    reward = 3;
    noreward = 4;
    shift1 = 1;
    shift2 = 2;
end
if errorType == 1
    % Loop through rows of trial_type
    for i = 2:size(session_ts,1)
        if (trial_type(i,:) == reward || trial_type(i,:) == noreward) && trial_type(i-1,:) == reward
            
            index_array(end+1) = i; % add current index to index array
        end
    end
elseif errorType == 2
    % Loop through rows of trial_type
    for i = 2:size(session_ts,1)
        if (trial_type(i,:) == shift1 || trial_type(i,:) == shift2) && trial_type(i-1,:) == reward
            
            index_array(end+1) = i; % add current index to index array
        end
    end
elseif errorType == 3
    % Loop through rows of trial_type
    for i = 2:size(session_ts,1)
        if (trial_type(i,:) == reward || trial_type(i,:) == noreward) && trial_type(i-1,:) == noreward
            
            index_array(end+1) = i; % add current index to index array
        end
    end
elseif errorType == 4
% Loop through rows of trial_type
for i = 2:size(session_ts,1)
    if (trial_type(i,:) == shift1 || trial_type(i,:) == shift2) && trial_type(i-1,:) == noreward
        
        index_array(end+1) = i; % add current index to index array
    end
end
else
    disp('Select a valid error probability to grab')
end
if isempty(index_array)
    errorProbTS = [];
else
    % Extract values from input2 using index_array
    errorProbTS = session_ts(index_array,:);
end