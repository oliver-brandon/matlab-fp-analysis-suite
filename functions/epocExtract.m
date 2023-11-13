function [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
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

% epocSTREAM = [];
% epocAMP = [];
% epocAUC = [];
epoc = TTLarray;
idx = find(ts1>-0.5,1);
for ii = 1:height(epoc)
    % if epoc(ii) == 0 || isempty(epoc)
    %     epocSTREAM(ii,1:minArrayLen) = NaN;
    %     epocAMP(ii) = NaN;
    %     epocAUC(ii) = NaN;
    %     break
    % end
    
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
[~,ampSt] = min(abs(ts1 - (ampWindow(1))));
[~,ampEn] = min(abs(ts1 - (ampWindow(2))));
for j = 1:height(streams_raw)
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
    
    % amplitude
    ampdFF(j) = max(streams_dFF(j,ampSt:ampEn));
    ampZ(j) = max(streams_z(j,ampSt:ampEn));

    % Calculate AUC above x=0 (dFF) %
    positive_indices = streams_dFF(j,:) > 0;
    y_pos = streams_dFF(j,positive_indices);
    x_pos = ts1(1,positive_indices);
    aucZ(j) = trapz(x_pos,y_pos);

    % Calculate AUC above x=0 (dFF) %
    positive_indices = streams_z(j,:) > 0;
    y_pos = streams_z(j,positive_indices);
    x_pos = ts1(1,positive_indices);
    aucdFF(j) = trapz(x_pos,y_pos);
end
epocSTREAM = streams_z;
epocAMP = ampZ;
epocAUC = aucZ;