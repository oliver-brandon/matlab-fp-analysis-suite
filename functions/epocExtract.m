function [epocSTREAM,epocAMP,epocAUC] = epocExtract( ...
    sessionSignal, ...
    sessionTime,...
    TTLarray, ...
    preEpochWindow, ...
    epocWindow, ...
    zBaseWindow, ...
    ampWindow, ...
    minArrayLen ...
    )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    epocSTREAM = [];
    epocAMP = [];
    epocAUC = [];
    for e = 1:height(TTLarray)
            e1 = TTLarray(e,1)-preEpochWindow;
            e2 = e1+epocWindow+preEpochWindow;
            if isempty(TTLarray)
                epocSTREAM(1,1:minArrayLen) = zeros;
                epocAMP = zeros(1);
                epocAUC = zeros(1);
                break
            elseif TTLarray == 0
                epocSTREAM(1,1:minArrayLen) = zeros;
                epocAMP = zeros(1);
                epocAUC = zeros(1);
                break
            end
            e3 = TTLarray(e,1)-zBaseWindow(:,2);
            e4 = TTLarray(e,1)-zBaseWindow(:,1);
            [~,ind1] = min(abs(sessionTime-e1));
            [~,ind2] = min(abs(sessionTime-e2));
            [~,ind3] = min(abs(sessionTime-e3));
            [~,ind4] = min(abs(sessionTime-e4));
            pre_signal = sessionSignal(1,ind4:ind3);
            epoc_signal = sessionSignal(1,ind1:ind2);
            epoc_time = sessionTime(1,ind1:ind2);
            [~,x] = find(epoc_time>ampWindow(:,1));
            [~,y] = find(epoc_time<ampWindow(:,2));
            windowAMP = [x,y];

            zb = mean(pre_signal);
            zsd = std(pre_signal);
            zfinal = (epoc_signal - zb)/zsd;
            if length(zfinal) < minArrayLen
                continue
            end
            epocSTREAM(e,:) = zfinal(1:minArrayLen);
            epocAMP(e,:) = max(zfinal(1, windowAMP(:,1):windowAMP(:,2)));
            epocAUC(e,:) = trapz(epoc_time,zfinal);
    end
end