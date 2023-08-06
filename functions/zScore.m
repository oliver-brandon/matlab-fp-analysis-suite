function [zScoreSignal] = zScore(signal)
%% Computes z-score for a given signal %%
%% (signal - mean signal) / standard deviation of signal %%

meanSignal = mean(signal);
stdSignal = std(signal);
zScoreSignal = (signal - meanSignal) / stdSignal;

end