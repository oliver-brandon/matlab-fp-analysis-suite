function [session_ts,trial_type,trial_name,lever_ts] = sessionArraySort(cueTS,A1,A2,A3,A4)
    
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
    
    cue_lever_ts(:,1) = cueTS;
   
    session_ts(2,:) = [];
    cue_lever_ts(:,2) = session_ts;
    lever_ts = cue_lever_ts(:,2) - cue_lever_ts(:,1);
    % Initialize the new cell array with the same size as the original column
    trial_name = cell(size(trial_type));
    
    % Loop through the original column and replace each number with the corresponding string
    for i = 1:length(trial_type)
        switch trial_type(i)
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