clearvars -except input;

% pull out choice/outcome columns
choice  = input(:,1);
outcome = input(:,2);
nTrials = numel(choice);

% how many full 30‐trial blocks?
nBlocks = floor(nTrials/30);

% Preallocate result array: [percent_cor, per, reg, percent_wS, percent_lS, percent_wS_inc, percent_lS_inc]
results = nan(nBlocks,9);

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
    
    % 1) percent correct
    percent_cor = sum(c)/T;
    
    % 2) perseveration: count zeros until first run of three 1's
    tillthree = T+1;
    for i = 1:T-2
        if c(i)==0
            per = per + 1;
        end
        if c(i)==1 && c(i+1)==1 && c(i+2)==1
            tillthree = i+2;
            break
        end
    end
    if percent_cor == 1
        per = 0;
    else
        per = per;
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
    else
        reg = reg;
    end
    % 4) win‐stay (correct lever)
    for i = 1:T-1
        if c(i)==1 && o(i)==1
            win_tot = win_tot + 1;
            if c(i+1)==1, wS = wS + 1; end
        end
    end
    percent_wS = wS / max(win_tot,1);
    
    % 5) lose‐shift (correct lever)
    for i = 1:T-1
        if c(i)==1 && o(i)==0
            lose_tot = lose_tot + 1;
            if c(i+1)==0, lS = lS + 1; end
        end
    end
    percent_lS = lS / max(lose_tot,1);
    
    % 6) win‐stay (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==1
            win_tot_inc = win_tot_inc + 1;
            if c(i+1)==0, wS_inc = wS_inc + 1; end
        end
    end
    percent_wS_inc = wS_inc / max(win_tot_inc,1);

    percent_wS_both = (wS + wS_inc) / (win_tot + win_tot_inc);
    
    % 7) lose‐shift (incorrect lever)
    for i = 1:T-1
        if c(i)==0 && o(i)==0
            lose_tot_inc = lose_tot_inc + 1;
            if c(i+1)==1, lS_inc = lS_inc + 1; end
        end
    end
    percent_lS_inc = lS_inc / max(lose_tot_inc,1);

    percent_lS_both = (lS + lS_inc) / (lose_tot + lose_tot_inc);
    
    % save into results matrix
    results(b,:) = [ percent_cor, per, reg, ...
                     percent_wS, percent_lS, ...
                     percent_wS_inc, percent_lS_inc, percent_wS_both, ...
                     percent_lS_both ];
end

% convert to a table
outputTable = array2table( results, ...
    'VariableNames', { ...
      'percent_cor', 'perseveration', 'regressive', ...
      'winStay',     'loseShift',     ...
      'winStay_inc', 'loseShift_inc', 'winStay_both', 'loseShift_both'} );

% display first few rows
disp( head(outputTable) )