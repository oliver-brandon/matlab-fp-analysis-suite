function [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
    rawSignal, ...
    sessionTime,...
    ts1,...
    TTLarray, ...
    preEpochWindow, ...
    postEpocWindow, ...
    baseWindow, ...
    ampWindow, ...
    epocArrayLen ...
    )

epocSTREAM = [];
epocAMP = [];
epocAUC = [];
for i = 1:height(TTLarray)
    if isempty(TTLarray)
        epocSTREAM(1,1:epocArrayLen) = zeros;
        epocAMP = zeros(1);
        epocAUC = zeros(1);
        break
    elseif TTLarray == 0
        epocSTREAM(1,1:epocArrayLen) = zeros;
        epocAMP = zeros(1);
        epocAUC = zeros(1);
        break
    end
    windowStart = TTLarray(i,1)-preEpochWindow;
    windowEnd = windowStart+postEpocWindow+preEpochWindow;
    [~,windSt] = min(abs(sessionTime-windowStart));
    [~,windEn] = min(abs(sessionTime-windowEnd));
    epocSigRaw = rawSignal(1,windSt:windEn);
    if length(epocSigRaw) < epocArrayLen
        mn = mean(epocSigRaw(1,end-10:end));
        epocSigRaw(1,end:epocArrayLen) = mn;
    elseif length(epocSigRaw) > epocArrayLen
        op = length(epocSigRaw);
        arrayDif = op - epocArrayLen;
        epocSigRaw = epocSigRaw(1,1:end-arrayDif);
    end
    streams_raw(i,1:epocArrayLen) = epocSigRaw;
end
    

    [~,baseSt] = min(abs(ts1 - baseWindow(1)));
    [~,baseEn] = min(abs(ts1 - baseWindow(2)));
    [~,ampSt] = min(abs(ts1 - (ampWindow(1))));
    [~,ampEn] = min(abs(ts1 - (ampWindow(2))));
   
for j = 1:height(streams_raw)
    % dF/F
    meanBase = mean(streams_raw(j,baseSt:baseEn));
    stdBase = std(streams_raw(j,baseSt:baseEn));
    streams_dFF(j,1:epocArrayLen) = streams_raw(j,1:epocArrayLen) - meanBase;
    streams_dFF(j,1:epocArrayLen) = 100*(streams_dFF(j,1:epocArrayLen) / meanBase);
    % z-Score
    meanBase_dFF = mean(streams_dFF(j,baseSt:baseEn));
    stdBase_dFF = std(streams_dFF(j,baseSt:baseEn));
    streams_z(j,1:epocArrayLen) = (streams_dFF(j,1:epocArrayLen) - meanBase_dFF) / stdBase_dFF;
    % amplitude
    ampdFF(j) = max(streams_dFF(j,ampSt:ampEn));
    ampZ(j) = max(streams_z(j,ampSt:ampEn));
    % AUC
    aucdFF(j) = trapz(ts1,streams_dFF);
    aucZ(j) = trapz(ts1,streams_z);
end
    
    % epocSTREAM = mean(streams_z,1);
    % epocAMP = mean(ampZ,1);
    % epocAUC = mean(aucZ,1);

    epocSTREAM = streams_z;
    epocAMP = ampZ;
    epocAUC = aucZ;
