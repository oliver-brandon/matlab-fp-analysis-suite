clear
warning off MATLAB:table:ModifiedAndSavedVarnames
% Define the path to the folder containing the Excel files
folder_path = '/Users/brandon/DA_PREF/AM251_Collab/Preference_AM251_RawData/'; % include back/forward slash at end of path
save_path = '/Users/brandon/DA_PREF/AM251_Collab/Preference_AM251_ExtractedTS/';
TSfilter = 3; % removes start/stop times that are less than 'TSfilter' seconds long

% Get a list of all the Excel files in the folder
file_list = dir(fullfile(folder_path, '*.xlsx'));
file_list = file_list(~startsWith({file_list.name},{'.','..','._','~'}));
numFiles = length(file_list);
% Loop through each Excel file in the folder
tic
for i = 1:numFiles
    % Get the filename and full file path of the current Excel file
    file_name = file_list(i).name;
    file_path = fullfile(folder_path, file_name);
    fprintf('Loading file %d of %d\n',i,numFiles)
    % Read the data from the Excel file
    % Read the cell data as a matrix
    trialInfo = readtable(file_path, 'Range', 'A34:B36');
    % Extract the string from the matrix
    animalID = char(trialInfo{1,2}); % Convert to a character array
    data = readtable(file_path, 'Range', 'B37:R73000');
    
    % Initialize variables to store the start and stop times for activity
    L_start_times = [];
    L_stop_times = [];
    disp('Finding timestamps for left cup.')
    % Loop through each row for Left Cup
    for j = 1:size(data, 1)
        if data{j, 'InZone_LeftCup_Center_point_'} == 0 && data{j+1, 'InZone_LeftCup_Center_point_'} == 1
           L_start_times = [L_start_times, data{j+1, 'RecordingTime'}];
        end   
        if data{j, 'InZone_LeftCup_Center_point_'} == 1 && data{j+1, 'InZone_LeftCup_Center_point_'} == 0
           L_stop_times = [L_stop_times, data{j, 'RecordingTime'}];
        end
        if j+1 == size(data, 1)
                break
        end     
    end
    R_start_times = [];
    R_stop_times = [];
    disp('Finding timestamps for right cup.')
    for j = 1:size(data, 1)
        if data{j, 'InZone_RightCup_Center_point_'} == 0 && data{j+1, 'InZone_RightCup_Center_point_'} == 1
            % record the start time
            R_start_times = [R_start_times, data{j+1, 'RecordingTime'}];
        end   
        if data{j, 'InZone_RightCup_Center_point_'} == 1 && data{j+1, 'InZone_RightCup_Center_point_'} == 0
            R_stop_times = [R_stop_times, data{j, 'RecordingTime'}];
        end
        if j+1 == size(data, 1)
            break
        end     
    end
    L_start_times = L_start_times';
    L_stop_times = L_stop_times';
    R_start_times = R_start_times';
    R_stop_times = R_stop_times';
    
    fprintf('Filtering out start/stop timestamps less than %d seconds\n',TSfilter)
    % Initialize variables to store filtered start and stop times for left
    % cup
    L_filtered_start_times = [];
    L_filtered_stop_times = [];
    % Loop through each start time for left cup
    for k = 1:length(L_start_times)
        % Check if there is a corresponding stop time
        if k < length(L_start_times) && L_stop_times(k) >= L_start_times(k)
            % Check if the stop time is at least 10 seconds after the start time
            if L_stop_times(k) - L_start_times(k) >= TSfilter
                % Add the start and stop times to the filtered variables
                L_filtered_start_times = [L_filtered_start_times; L_start_times(k)];
                L_filtered_stop_times = [L_filtered_stop_times; L_stop_times(k)];
            end
        % If this is the last start time, add the corresponding stop time to the filtered variables
        elseif k == length(L_start_times)
            break
        end
    end
    % Initialize variables to store filtered start and stop times for right
    % cup
    R_filtered_start_times = [];
    R_filtered_stop_times = [];
    % Loop through each start time for right cup
    for k = 1:length(R_start_times)
        % Check if there is a corresponding stop time
        if k < length(R_start_times) && R_stop_times(k) >= R_start_times(k)
            % Check if the stop time is at least 10 seconds after the start time
            if R_stop_times(k) - R_start_times(k) >= TSfilter
                % Add the start and stop times to the filtered variables
                R_filtered_start_times = [R_filtered_start_times; R_start_times(k)];
                R_filtered_stop_times = [R_filtered_stop_times; R_stop_times(k)];
            end
        % If this is the last start time, add the corresponding stop time to the filtered variables
        elseif k == length(R_start_times)
            break
        end
    end
    disp('Writing filtered timestamps to new csv files.')

    % Create a table with the timestamp data
    L = table(L_filtered_start_times, L_filtered_stop_times, ...
        'VariableNames', {'Left Start', 'Left Stop'});
    R = table(R_filtered_start_times, R_filtered_stop_times, ...
        'VariableNames', {'Right Start', 'Right Stop'});


   
    
    L_new_filename = strcat(save_path,animalID,'_','leftCup','.csv'); % Add file extension
    % Write the table to a new sheet in the Excel file
    writetable(L, L_new_filename);
    R_new_filename = strcat(save_path,animalID,'_','rightCup','.csv'); % Add file extension
    % Write the table to a new sheet in the Excel file
    writetable(R, R_new_filename);
end
toc
fprintf("Done! Files analyzed: %d\n", numFiles)
NERD_STATS(toc,numFiles);
