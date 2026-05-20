function idx = timeToIndex(targetTime, timeStart, fs, nSamples)
% timeToIndex Convert seconds to bounded sample indices.

idx = round((targetTime - timeStart) .* fs) + 1;
idx = max(1, min(nSamples, idx));

end
