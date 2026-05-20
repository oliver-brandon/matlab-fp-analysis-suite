clearvars -except input;

% pull out choice/outcome columns
choice  = input(:,1);
outcome = input(:,2);
nTrials = numel(choice);

% how many full 30-trial blocks?
nBlocks = floor(nTrials/30);

% Preallocate result array:
% [percent_cor, per, reg, winStay, loseShift, winStay_inc, loseShift_inc, winStay_both, loseShift_both,
%  winShift, loseStay, winShift_inc, loseStay_inc, winShift_both, loseStay_both,
%  nextCorrect_afterFirstCorrectReward, percentCorrect_afterFirstCorrectReward, trialsTo3Correct]
results = nan(nBlocks,18);

for b = 1:nBlocks
    % indices for this block
    i0 = (b-1)*30 + 1;
    i1 = b*30;

    c = choice(i0:i1);
    o = outcome(i0:i1);
    T = numel(c);

    % initialize counters
    per = 0;     reg = 0;
    wS  = 0;     win_tot  = 0;
    lS  = 0;     lose_tot = 0;
    wS_inc = 0;  win_tot_inc  = 0;
    lS_inc = 0;  lose_tot_inc = 0;
    wShift = 0;  lStay = 0;
    wShift_inc = 0;  lStay_inc = 0;

    % 1) percent correct
    percent_cor = sum(c)/T;

    % 2) perseveration: count zeros until first run of three 1's
    for i = 1:T-2
        if c(i)==0
            per = per + 1;
        end
        if c(i)==1 && c(i+1)==1 && c(i+2)==1
            break
        end
    end
    if percent_cor == 1
        per = 0;
    end

    % 3) regressive errors: zeros after first run of three 1's
    for i = 1:T-2
        if c(i)==1 && c(i+1)==1 && c(i+2)==1
            reg = sum(c(i+3:end)==0);
            break
        end
    end
    if percent_cor == 1
        reg = 0;
    end


    %% 3b) trials to 3 consecutive correct (within block)
    trialsTo3Correct = NaN;
    for i = 1:(T-2)
        if c(i)==1 && c(i+1)==1 && c(i+2)==1
            trialsTo3Correct = i+2;
            break;
        end
    end

    % 4) win-stay (correct lever)
    for i = 1:T-1
        if c(i)==1 && o(i)==1
            win_tot = win_tot + 1;
            if c(i+1)==1, wS = wS + 1; end
        end
    end
    percent_wS = wS / max(win_tot,1);

    % 5) lose-shift (correct lever)
    for i = 1:T-1
        if c(i)==1 && o(i)==0
            lose_tot = lose_tot + 1;
            if c(i+1)==0, lS = lS + 1; end
        end
    end
    percent_lS = lS / max(lose_tot,1);


    % 5b) win-shift (correct lever): rewarded then switch levers next trial
    for i = 1:T-1
        if c(i)==1 && o(i)==1
            if c(i+1)==0, wShift = wShift + 1; end
        end
    end
    percent_wShift = wShift / max(win_tot,1);

    % 5c) lose-stay (correct lever): not rewarded then stay on same lever next trial
    for i = 1:T-1
        if c(i)==1 && o(i)==0
            if c(i+1)==1, lStay = lStay + 1; end
        end
    end
    percent_lStay = lStay / max(lose_tot,1);

    % 6) win-stay (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==1
            win_tot_inc = win_tot_inc + 1;
            if c(i+1)==0, wS_inc = wS_inc + 1; end
        end
    end
    percent_wS_inc = wS_inc / max(win_tot_inc,1);

    percent_wS_both = (wS + wS_inc) / max((win_tot + win_tot_inc),1);

    % 7) lose-shift (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==0
            lose_tot_inc = lose_tot_inc + 1;
            if c(i+1)==1, lS_inc = lS_inc + 1; end
        end
    end
    percent_lS_inc = lS_inc / max(lose_tot_inc,1);


    % 7b) win-shift (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==1
            if c(i+1)==1, wShift_inc = wShift_inc + 1; end
        end
    end
    percent_wShift_inc = wShift_inc / max(win_tot_inc,1);

    % 7c) lose-stay (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==0
            if c(i+1)==0, lStay_inc = lStay_inc + 1; end
        end
    end
    percent_lStay_inc = lStay_inc / max(lose_tot_inc,1);

    percent_wShift_both = (wShift + wShift_inc) / max((win_tot + win_tot_inc),1);
    percent_lStay_both  = (lStay  + lStay_inc)  / max((lose_tot + lose_tot_inc),1);

    percent_lS_both = (lS + lS_inc) / max((lose_tot + lose_tot_inc),1);

    % 8) NEW: indicator that the trial AFTER the first correct+rewarded trial is correct
    % and percent correct on the remaining trials after that first correct+rewarded trial.
    idxFirstCR = find(c==1 & o==1, 1, 'first');
    if isempty(idxFirstCR)
        nextCorrect_afterFirstCR = 0;
        percentCorrect_afterFirstCR = 0;
    else
        if idxFirstCR < T
            nextCorrect_afterFirstCR = double(c(idxFirstCR+1)==1);
            percentCorrect_afterFirstCR = mean(c(idxFirstCR+1:end));
        else
            % first correct+rewarded trial happened on the last trial of the block
            nextCorrect_afterFirstCR = 0;
            percentCorrect_afterFirstCR = 0;
        end
    end

    % save into results matrix
    results(b,:) = [ percent_cor, per, reg, ...
                     percent_wS, percent_lS, ...
                     percent_wS_inc, percent_lS_inc, percent_wS_both, ...
                     percent_lS_both, ...
                     percent_wShift, percent_lStay, ...
                     percent_wShift_inc, percent_lStay_inc, ...
                     percent_wShift_both, percent_lStay_both, ...
                     nextCorrect_afterFirstCR, percentCorrect_afterFirstCR, trialsTo3Correct ];
end

% convert to a table
outputTable = array2table( results, ...
    'VariableNames', { ...
      'percent_cor', 'perseveration', 'regressive', ...
      'winStay',     'loseShift',     ...
      'winStay_inc', 'loseShift_inc', 'winStay_both', 'loseShift_both', ...
      'winShift',    'loseStay',      ...
      'winShift_inc','loseStay_inc',  ...
      'winShift_both','loseStay_both', ...
      'nextCorrect_afterFirstCR', 'percentCorrect_afterFirstCR', 'trialsTo3Correct'} );

% display first few rows
if exist('head','file') == 2
    disp( head(outputTable) )
else
    disp( outputTable(1:min(8,height(outputTable)),:) )
end
