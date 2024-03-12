function [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(epocOnset, timeout)
    % separateTimestamps separates lever press epocOnset into two categories.
    % Inputs:
    %   epocOnset - a column vector of epocOnset.
    % Outputs:
    %   rewardTimestamps - epocOnset that resulted in a reward.
    %   timeoutTimestamps - epocOnset during the timeout period.

    % Initialize the output vectors
    rewardTimestamps = [];
    rewardTimeout = [];
    timeoutTimestamps = [];

    % Check if the epocOnset vector is empty
    if isempty(epocOnset)
        return;
    end

    % The first timestamp always results in a reward
    lastRewardTimestamp = epocOnset(1);
    rewardTimestamps = [rewardTimestamps; lastRewardTimestamp];

    % Iterate through the rest of the epocOnset
    for i = 2:length(epocOnset)
        currentTime = epocOnset(i);
        
        if (currentTime - lastRewardTimestamp) <= 3
            rewardTimeout = [rewardTimeout; currentTime];
        % Check if the current timestamp is during the timeout period
        elseif (currentTime - lastRewardTimestamp) > (timeout + 3)
            % More than timout seconds have passed, count as a reward
            rewardTimestamps = [rewardTimestamps; currentTime];
            lastRewardTimestamp = currentTime;
        else
            % Within the timeout-second timeout, no reward
            timeoutTimestamps = [timeoutTimestamps; currentTime];
        end
    end
end

