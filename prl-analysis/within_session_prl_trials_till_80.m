clearvars -except input
window = 10;
threshold = 0.80;
minCorrect = ceil(window * threshold);
% Copy/paste into command window: input = [];
% Paste choice data into input
%   - Assumes input is correctness (1 = correct, 0 = incorrect)
choices = input(:,1);
nTrials = length(choices);

nWindows = nTrials - window + 1;
windowAcc = nan(nWindows,1);

trial80_exact = NaN;   % this will hold the exact trial

for k = 1:nWindows
    
    currentWindow = choices(k:k+window-1);
    windowAcc(k) = mean(currentWindow);
    
    if round(windowAcc(k),2) >= threshold
        
        % Find the trial inside this window where
        % cumulative correct reaches the minimum required
        cumulativeCorrect = cumsum(currentWindow);
        idxWithinWindow = find(cumulativeCorrect >= minCorrect, 1, 'first');
        
        trial80_exact = k + idxWithinWindow - 1;
        break
    end
end

if isnan(trial80_exact)
    fprintf('80%% criterion was NOT achieved.\n');
else
    fprintf('80%% criterion first achieved on trial %d\n', trial80_exact);
end
choices = [];