function [data, errorProbLeverTS, errorProbCueTS] = errorProbExtract(...
    data,...
    session_identifiers, ...
    errorType, ...
    lever ...
    )
% This function will extract lever and cue timestamps for win-stay, win-shift, lose-stay,
% and lose-shift error probabilities. It will output two arrays:
% errorProbLeverTS and errorProbCueTS. Column 1 is the timestamps
% associated with the win/lose trial and column 2 is the timestamps
% associated with the stay/shift trial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Args:
%% session_identifiers: 
% An Nx2 array with cue and lever time stamps in order
% of occurance in column 1 and a number associated with each timestamp in
% column 2 that represents the TTL (0 = cue, 1 = cRew, 2 = cNoRew, 3 =
% iRew, 4 = iNoRew). Can be extracted using the function sessionArraySort.m
%% errorType: 
% 1 = winStay, 2 = winShift, 3 = loseStay, 4 = loseShift
%% lever: 
% errors based on correct lever = 1, incorrect lever = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize index array
win_index_array = [];
lose_index_array = [];
stay_index_array = [];
shift_index_array = [];
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
    for i = 4:2:size(session_identifiers,1)
        if (session_identifiers(i,2) == reward || session_identifiers(i,2) == noreward) && session_identifiers(i-2,2) == reward
            
            win_index_array(end+1,1) = i - 2;
            stay_index_array(end+1,1) = i; % add current index to index array
        end
    end
    % Extract win-stay timestamps
    errorProbLeverTS(:,1) = session_identifiers(win_index_array,1);
    errorProbLeverTS(:,2) = session_identifiers(stay_index_array,1);
    errorProbCueTS(:,1) = session_identifiers(win_index_array-1,1);
    errorProbCueTS(:,2) = session_identifiers(stay_index_array-1,1);
elseif errorType == 2
    % Loop through rows of trial_type
    for i = 4:2:size(session_identifiers,1)
        if (session_identifiers(i,2) == shift1 || session_identifiers(i,2) == shift2) && session_identifiers(i-2,2) == reward
            
            win_index_array(end+1,1) = i - 2;
            shift_index_array(end+1,1) = i; % add current index to index array
        end
    end
    % Extract win-shift timestamps
    errorProbLeverTS(:,1) = session_identifiers(win_index_array,1);
    errorProbLeverTS(:,2) = session_identifiers(shift_index_array,1);
    errorProbCueTS(:,1) = session_identifiers(win_index_array-1,1);
    errorProbCueTS(:,2) = session_identifiers(shift_index_array-1,1);
elseif errorType == 3
    % Loop through rows of trial_type
    for i = 4:2:size(session_identifiers,1)
        if (session_identifiers(i,2) == reward || session_identifiers(i,2) == noreward) && session_identifiers(i-2,2) == noreward
            
            lose_index_array(end+1,1) = i - 2;
            stay_index_array(end+1,1) = i; % add current index to index array
        end
    end
    % Extract lose-stay timestamps
    errorProbLeverTS(:,1) = session_identifiers(lose_index_array,1);
    errorProbLeverTS(:,2) = session_identifiers(stay_index_array,1);
    errorProbCueTS(:,1) = session_identifiers(lose_index_array-1,1);
    errorProbCueTS(:,2) = session_identifiers(stay_index_array-1,1);
elseif errorType == 4
    % Loop through rows of trial_type
    for i = 4:2:size(session_identifiers,1)
        if (session_identifiers(i,2) == shift1 || session_identifiers(i,2) == shift2) && session_identifiers(i-2,2) == noreward
            
            lose_index_array(end+1,1) = i - 2;
            shift_index_array(end+1,1) = i; % add current index to index array
        end
    end
    % Extract lose-shift timestamps
    errorProbLeverTS(:,1) = session_identifiers(lose_index_array,1);
    errorProbLeverTS(:,2) = session_identifiers(shift_index_array,1);
    errorProbCueTS(:,1) = session_identifiers(lose_index_array-1,1);
    errorProbCueTS(:,2) = session_identifiers(shift_index_array-1,1);
else
    disp('Select a valid error probability to grab')
end

% Creates new TTLs
if errorType == 1
    % lever TTLs
    data.epocs.WIN_stay_LEV.name = 'WIN_stay_LEV';
    data.epocs.WIN_stay_LEV.onset = errorProbLeverTS(:,1);
    data.epocs.WIN_stay_LEV.offset = errorProbLeverTS(:,1) + 1;
    data.epocs.WIN_stay_LEV.data = ones(height(errorProbLeverTS(:,1)))*5;
    data.epocs.win_STAY_LEV.name = 'win_STAY_LEV';
    data.epocs.win_STAY_LEV.onset = errorProbLeverTS(:,2);
    data.epocs.win_STAY_LEV.offset = errorProbLeverTS(:,2) + 1;
    data.epocs.win_STAY_LEV.data = ones(height(errorProbLeverTS(:,2)))*6;
    % cue TTLs
    data.epocs.WIN_stay_CUE.name = 'WIN_stay_CUE';
    data.epocs.WIN_stay_CUE.onset = errorProbCueTS(:,1);
    data.epocs.WIN_stay_CUE.offset = errorProbCueTS(:,1) + 1;
    data.epocs.WIN_stay_CUE.data = ones(height(errorProbCueTS(:,1)))*7;
    data.epocs.win_STAY_CUE.name = 'win_STAY_CUE';
    data.epocs.win_STAY_CUE.onset = errorProbCueTS(:,2);
    data.epocs.win_STAY_CUE.offset = errorProbCueTS(:,2) + 1;
    data.epocs.win_STAY_CUE.data = ones(height(errorProbCueTS(:,2)))*8;
elseif errorType == 2
    % lever TTLs
    data.epocs.WIN_shift_LEV.name = 'WIN_shift_LEV';
    data.epocs.WIN_shift_LEV.onset = errorProbLeverTS(:,1);
    data.epocs.WIN_shift_LEV.offset = errorProbLeverTS(:,1) + 1;
    data.epocs.WIN_shift_LEV.data = ones(height(errorProbLeverTS(:,1))) * 9;
    data.epocs.win_SHIFT_LEV.name = 'win_SHIFT_LEV';
    data.epocs.win_SHIFT_LEV.onset = errorProbLeverTS(:,2);
    data.epocs.win_SHIFT_LEV.offset = errorProbLeverTS(:,2) + 1;
    data.epocs.win_SHIFT_LEV.data = ones(height(errorProbLeverTS(:,2))) * 10;
    % cue TTLs
    data.epocs.WIN_shift_CUE.name = 'WIN_shift_CUE';
    data.epocs.WIN_shift_CUE.onset = errorProbCueTS(:,1);
    data.epocs.WIN_shift_CUE.offset = errorProbCueTS(:,1) + 1;
    data.epocs.WIN_shift_CUE.data = ones(height(errorProbCueTS(:,1))) * 11;
    data.epocs.win_SHIFT_CUE.name = 'win_SHIFT_CUE';
    data.epocs.win_SHIFT_CUE.onset = errorProbCueTS(:,2);
    data.epocs.win_SHIFT_CUE.offset = errorProbCueTS(:,2) + 1;
    data.epocs.win_SHIFT_CUE.data = ones(height(errorProbCueTS(:,2))) * 12;
elseif errorType == 3
    % lever TTLs
    data.epocs.LOS_stay_LEV.name = 'LOS_stay_LEV';
    data.epocs.LOS_stay_LEV.onset = errorProbLeverTS(:,1);
    data.epocs.LOS_stay_LEV.offset = errorProbLeverTS(:,1) + 1;
    data.epocs.LOS_stay_LEV.data = ones(height(errorProbLeverTS(:,1))) * 13;
    data.epocs.los_STAY_LEV.name = 'los_STAY_LEV';
    data.epocs.los_STAY_LEV.onset = errorProbLeverTS(:,2);
    data.epocs.los_STAY_LEV.offset = errorProbLeverTS(:,2) + 1;
    data.epocs.los_STAY_LEV.data = ones(height(errorProbLeverTS(:,2))) * 14;
    % cue TTLs
    data.epocs.LOS_stay_CUE.name = 'LOS_stay_CUE';
    data.epocs.LOS_stay_CUE.onset = errorProbCueTS(:,1);
    data.epocs.LOS_stay_CUE.offset = errorProbCueTS(:,1) + 1;
    data.epocs.LOS_stay_CUE.data = ones(height(errorProbCueTS(:,1))) * 15;
    data.epocs.los_STAY_CUE.name = 'los_STAY_CUE';
    data.epocs.los_STAY_CUE.onset = errorProbCueTS(:,2);
    data.epocs.los_STAY_CUE.offset = errorProbCueTS(:,2) + 1;
    data.epocs.los_STAY_CUE.data = ones(height(errorProbCueTS(:,2))) * 16;
elseif errorType == 4
    % lever TTLs
    data.epocs.LOS_shift_LEV.name = 'LOS_shift_LEV';
    data.epocs.LOS_shift_LEV.onset = errorProbLeverTS(:,1);
    data.epocs.LOS_shift_LEV.offset = errorProbLeverTS(:,1) + 1;
    data.epocs.LOS_shift_LEV.data = ones(height(errorProbLeverTS(:,1))) * 17;
    data.epocs.los_SHIFT_LEV.name = 'los_SHIFT_LEV';
    data.epocs.los_SHIFT_LEV.onset = errorProbLeverTS(:,2);
    data.epocs.los_SHIFT_LEV.onset = errorProbLeverTS(:,2) + 1;
    data.epocs.los_SHIFT_LEV.data = ones(height(errorProbLeverTS(:,2))) * 18;
    % cue TTLs
    data.epocs.LOS_shift_CUE.name = 'LOS_shift_CUE';
    data.epocs.LOS_shift_CUE.onset = errorProbCueTS(:,1);
    data.epocs.LOS_shift_CUE.offset = errorProbCueTS(:,1) + 1;
    data.epocs.LOS_shift_CUE.data = ones(height(errorProbCueTS(:,1))) * 19;
    data.epocs.los_SHIFT_CUE.name = 'los_SHIFT_CUE';
    data.epocs.los_SHIFT_CUE.onset = errorProbCueTS(:,2);
    data.epocs.los_SHIFT_CUE.offset = errorProbCueTS(:,2) + 1;
    data.epocs.los_SHIFT_CUE.data = ones(height(errorProbCueTS(:,2))) * 20;
end