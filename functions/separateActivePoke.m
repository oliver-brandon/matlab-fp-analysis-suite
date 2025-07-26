function [rewardedNosepokes, nonrewardedNosepokes, timeoutNosepokes] = separateActivePoke(timestamps,FR)
    % Parameters
    rewardTimeout = 0; % seconds after reward delivery
    ITI = 10; % inter-trial interval in seconds

    % Initialize outputs
    rewardedNosepokes = [];
    nonrewardedNosepokes = [];
    timeoutNosepokes = [];

    idx = 2; % Start at the first nosepoke
    n = length(timestamps);

    while idx <= n
        % Start of the current trial is the first nosepoke not in timeout/ITI
        trialPokes = [];
        count = 0;

        % Collect FR number of active pokes (outside of timeout/ITI)
        while idx <= n && count < FR
            t = timestamps(idx);
            trialPokes = [trialPokes; t];
            count = count + 1;
            idx = idx + 1;
        end

        if count == FR
            % The last poke in this set is rewarded
            rewardedNosepokes = [rewardedNosepokes; trialPokes(end)];
            % The preceding pokes are non-rewarded
            if FR > 1
                nonrewardedNosepokes = [nonrewardedNosepokes; trialPokes(1:end-1)];
            end

            % Mark timeout and ITI window
            rewardTime = trialPokes(end);
            timeoutEnd = rewardTime + rewardTimeout;
            ITIEnd = timeoutEnd + ITI;

            % Collect pokes during timeout/ITI
            while idx <= n && timestamps(idx) < ITIEnd
                timeoutNosepokes = [timeoutNosepokes; timestamps(idx)];
                idx = idx + 1;
            end
        else
            % If we exit early (shouldn't happen in typical FR sessions), append all to nonrewarded
            nonrewardedNosepokes = [nonrewardedNosepokes; trialPokes];
        end
    end
end
