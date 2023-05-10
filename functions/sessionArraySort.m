function [session_identifiers,lever_session_ts,trial_number,trial_name,lever_latency] = sessionArraySort(...
    CUETTL,TTL1,TTL2,TTL3,TTL4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function is designed to extract and sort all TTL events in order
    % of occurance during a session. Currently set up to extract and sort
    % the following TTLs: 
    %
    % Cue (1)
    % Correct Rewarded Lever Press (2)
    % Correct Unrewarded Lever Press (3)
    % Incorrect Rewarded Lever Press (4)
    % Incorrect Unrewarded Lever Press (5)
    % 
    % The TTL variable being passed must be a column vector with timestamps
    % corresponding to the TTL event.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X0 = zeros(height(CUETTL),1);
    X1 = ones(height(TTL1),1);
    X2 = ones(height(TTL2),1)*2;
    X3 = ones(height(TTL3),1)*3;
    X4 = ones(height(TTL4),1)*4;
    % Append first column of each array to new array
    combined_array1 = [CUETTL(:,1); TTL1(:,1); TTL2(:,1); TTL3(:,1); TTL4(:,1)];
    combined_array2 = [X0(:,1); X1(:,1); X2(:,1); X3(:,1); X4(:,1)];
    
    session_identifiers = [combined_array1 combined_array2];
    % Sort the first column of newArray in ascending order
    session_identifiers = sortrows(session_identifiers, 1);
    % removes lever ts preceding missing second cue ts
    session_identifiers = [session_identifiers(1:2,:);session_identifiers(4:end,:)];
    % Find the rows in A where the first column is not zero
    idx = session_identifiers(:,1) ~= 0;
    % Use logical indexing to remove rows with zero values
    session_identifiers = session_identifiers(idx,:);
    
    lever_session_ts = session_identifiers(:,1);
    trial_number = session_identifiers(:,2);
    
    cue_lever_ts(:,1) = CUETTL;
   
    lever_session_ts(2,:) = [];
    
    if size(cue_lever_ts,1) > size(lever_session_ts,1)
        cue_lever_ts = cue_lever_ts(1:size(lever_session_ts,1), :);
    else
        cue_lever_ts(end+1:size(lever_session_ts,1), :) = 0;
    end

    cue_lever_ts(:,2) = lever_session_ts;
    lever_latency = cue_lever_ts(:,2) - cue_lever_ts(:,1);
    % Initialize the new cell array with the same size as the original column
    trial_name = cell(size(trial_number));
    
    % Loop through the original column and replace each number with the corresponding string
    for i = 1:length(trial_number)
        switch trial_number(i)
            case 1
                trial_name{i} = ['C+ Trial ' num2str(i)];
            case 2
                trial_name{i} = ['C- Trial ' num2str(i)];
            case 3
                trial_name{i} = ['I+ Trial ' num2str(i)];
            case 4
                trial_name{i} = ['I- Trial ' num2str(i)];
        end
    end

end