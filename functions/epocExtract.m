function [epocSTREAM] = epocExtract( ...
    sessionSignal, ...
    sessionTime,...
    TTLarray, ...
    preEpochWindow, ...
    epocWindow, ...
    zBaseWindow, ...
    ampWindow, ...
    minArrayLen, ...
    ts1...
    )

streams_z = [];
streams_dFF = [];
epoc = TTLarray(:);
streams_raw = NaN(height(epoc), minArrayLen);
if isempty(epoc)
    epocSTREAM = streams_z;
    return
end

fs = 1 / median(diff(sessionTime));
sessionStart = sessionTime(1);
nSamples = numel(sessionTime);
for ii = 1:height(epoc)
    if epoc(ii) == 0
        streams_raw(ii,1:minArrayLen) = NaN;
        continue
    end
    windowStart = epoc(ii)-preEpochWindow;
    windowEnd = windowStart+preEpochWindow+epocWindow;
    windSt = timeToIndex(windowStart, sessionStart, fs, nSamples);
    windEn = timeToIndex(windowEnd, sessionStart, fs, nSamples);
    if windEn < windSt
        streams_raw(ii,1:minArrayLen) = NaN;
        continue
    end
    epocSigRaw = sessionSignal(1,windSt:windEn);

    if length(epocSigRaw) < minArrayLen
        if isempty(epocSigRaw)
            mn = NaN;
        else
            tailStart = max(1, length(epocSigRaw) - 10);
            mn = mean(epocSigRaw(1,tailStart:end), 'omitnan');
        end
        epocSigRaw(1,end:minArrayLen) = mn;
    elseif length(epocSigRaw) > minArrayLen
        op = length(epocSigRaw);
        arrayDif = op - minArrayLen;
        epocSigRaw = epocSigRaw(1,1:end-arrayDif);
    end
    streams_raw(ii,1:minArrayLen) = epocSigRaw;
end
[~,baseSt] = min(abs(ts1 - (zBaseWindow(1))));
[~,baseEn] = min(abs(ts1 - (zBaseWindow(2))));

for j = 1:height(streams_raw)
    if all(isnan(streams_raw(j,:)))
        streams_dFF(j,1:minArrayLen) = NaN;
        streams_z(j,1:minArrayLen) = NaN;
        continue
    end
    % dF/F
    meanBase = mean(streams_raw(j,baseSt:baseEn), 'omitnan');
    if ~isfinite(meanBase) || meanBase == 0
        streams_dFF(j,1:minArrayLen) = NaN;
        streams_z(j,1:minArrayLen) = NaN;
        continue
    end
    streams_dFF(j,1:minArrayLen) = streams_raw(j,1:minArrayLen) - meanBase;
    streams_dFF(j,1:minArrayLen) = 100*(streams_dFF(j,1:minArrayLen) / meanBase);
    % z-Score
    meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn), 'omitnan');
    stdBase_dFF = std(streams_dFF(j,baseSt:baseEn), 0, 'omitnan');
    if ~isfinite(meanBase_dFF) || ~isfinite(stdBase_dFF) || stdBase_dFF == 0
        streams_z(j,1:minArrayLen) = NaN;
        continue
    end
    streams_z(j,1:minArrayLen) = (streams_dFF(j,1:minArrayLen) - meanBase_dFF) / stdBase_dFF;

end
epocSTREAM = streams_z;
