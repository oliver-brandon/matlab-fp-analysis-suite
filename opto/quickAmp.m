clearvars -except sigs
numSigs = size(sigs,2);
amplitudes = [];
for i = 1:numSigs
    amp = max(sigs(204:end,i));
    % Store the amplitude in an array for later use
    amplitudes(i,1) = amp;
end

sigs = [];