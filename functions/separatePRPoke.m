function [rewardedNosepokes, nonrewardedNosepokes, timeoutNosepokes] = separatePRPoke(activePort, rewards)
    % Parameters
    rewardTimeout = 2; % seconds after reward delivery
    ITI = 10; % inter-trial interval in seconds
    PRlist = [1,2,3,5,8,12,18,27,40,60,90,135,200,300,450,675,1000];

    % Initialize outputs
    nonrewardedNosepokes = [];
    timeoutNosepokes = [];
    
    % Find closest values from rewards to activePort
    rewardedNosepokes = zeros(size(rewards));
    for i = 1:length(rewards)
        [~, idx] = min(abs(activePort - rewards(i)));
        rewardedNosepokes(i) = activePort(idx);
    end
    
    
    % Find values in activePort that are less than or equal to ITI seconds after a value in rewards
    for i = 1:length(rewards)
        % Calculate the time window for each reward
        timeWindowStart = rewards(i);
        timeWindowEnd = timeWindowStart + ITI + rewardTimeout;
        
        % Find activePort values within the time window
        validNosepokes = activePort(activePort >= timeWindowStart & activePort <= timeWindowEnd);
        
        % Append to nonrewardedNosepokes
        timeoutNosepokes = [timeoutNosepokes; validNosepokes];
    end

        
        
    % Filter out values in activePort that exist in rewardedNosepokes and timeoutNosepokes
    nonrewardedNosepokes = activePort(~ismember(activePort, rewardedNosepokes) & ~ismember(activePort, timeoutNosepokes));
end
