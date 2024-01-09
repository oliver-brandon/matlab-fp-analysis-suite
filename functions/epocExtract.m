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
epoc = TTLarray;
idx = find(ts1>-0.5,1);
for ii = 1:height(epoc)
    if epoc(ii) == 0
        streams_raw(ii,1:minArrayLen) = NaN;
        break
    end
    windowStart = epoc(ii)-preEpochWindow;
    windowEnd = windowStart+preEpochWindow+epocWindow;
    [~,windSt] = min(abs(sessionTime - windowStart));
    [~,windEn] = min(abs(sessionTime - windowEnd));
    epocSigRaw = sessionSignal(1,windSt:windEn);

    if length(epocSigRaw) < minArrayLen
        mn = mean(epocSigRaw(1,end-10:end));
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
    if isnan(streams_raw)
        streams_raw(1,1:minArrayLen) = NaN;
    end
    % dF/F
    meanBase = mean(streams_raw(j,baseSt:baseEn));
    stdBase = std(streams_raw(j,baseSt:baseEn));
    streams_dFF(j,1:minArrayLen) = streams_raw(j,1:minArrayLen) - meanBase;
    streams_dFF(j,1:minArrayLen) = 100*(streams_dFF(j,1:minArrayLen) / meanBase);
    % z-Score
    meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
    stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
    streams_z(j,1:minArrayLen) = (streams_dFF(j,1:minArrayLen) - meanBase_dFF) / stdBase_dFF;
    % adjusts streams to baseline of zero at -0.5s %
    if streams_z(j,idx) < 0
        val = streams_z(j,idx);
        diff = 0 - val;
        streams_z(j,1:minArrayLen) = streams_z(j,1:minArrayLen) + abs(diff);
    elseif streams_z(j,idx) > 0
        val = streams_z(j,idx);
        diff = 0 - val;
        streams_z(j,1:minArrayLen) = streams_z(j,1:minArrayLen) - abs(diff);
    end

end
epocSTREAM = streams_z;