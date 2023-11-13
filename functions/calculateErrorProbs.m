function [errProb] = calculateErrorProbs(...
    errorProbLeverTS,...
    session_identifiers,...
    errorType,...
    lever)

if errorType == 1
    if lever == 1
        val1 = 1;
    else
        val1 = 3;
    end
elseif errorType == 2
    if lever == 1
        val1 = 1;
    else
        val1 = 3;
    end
elseif errorType == 3
    if lever == 1
        val1 = 2;
    else
        val1 = 4;
    end
elseif errorType == 4
    if lever == 1
        val1 = 2;
    else
        val1 = 4;
    end
end
errProb = height(errorProbLeverTS) / sum(session_identifiers(:,2) == val1);


end