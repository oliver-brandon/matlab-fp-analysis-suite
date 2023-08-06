function [signal_dFF] = deltaFF(signal)
%% Computes delta F/F0 for a given signal %%
%%     (signal - baseline) / baseline     %%

baseline = mean(signal); % calculate the mean baseline
deltaF = signal - baseline; % calculate the change in the signal
deltaFF = deltaF ./ baseline; % calculate delta F/F0
signal_dFF = deltaFF;

end