clearvars -except input;
choice = input(:,1);
outcome = input(:,2);

tot_rows = size(choice,1);
num_blocks = ceil(tot_rows / 30);

% Preallocate arrays: one value per block
primeGONOGO_full = nan(num_blocks, 1);
prime2AFC_full = nan(num_blocks, 1);

primeGONOGO_first = nan(num_blocks, 1);
prime2AFC_first = nan(num_blocks, 1);

primeGONOGO_second = nan(num_blocks, 1);
prime2AFC_second = nan(num_blocks, 1);

for block = 1:num_blocks
    idx_start = (block - 1) * 30 + 1;
    idx_end = min(block * 30, tot_rows);
    trials = choice(idx_start:idx_end);
    n_trials = idx_end - idx_start + 1;

    % Full block d'
    tot_correct = sum(trials);
    primeGONOGO_full(block) = norminv((tot_correct + 0.5) / (n_trials + 1)) ...
                            - norminv((n_trials - tot_correct + 0.5) / (n_trials + 1));
    prime2AFC_full(block) = sqrt(2) * norminv((tot_correct + 0.5) / (n_trials + 1));

    % First half
    half_point = floor(n_trials / 2);
    trials_first = trials(1:half_point);
    tot_correct_first = sum(trials_first);
    primeGONOGO_first(block) = norminv((tot_correct_first + 0.5) / (half_point + 1)) ...
                             - norminv((half_point - tot_correct_first + 0.5) / (half_point + 1));
    prime2AFC_first(block) = sqrt(2) * norminv((tot_correct_first + 0.5) / (half_point + 1));

    % Second half
    trials_second = trials(half_point+1:end);
    n_second = length(trials_second);
    tot_correct_second = sum(trials_second);
    primeGONOGO_second(block) = norminv((tot_correct_second + 0.5) / (n_second + 1)) ...
                              - norminv((n_second - tot_correct_second + 0.5) / (n_second + 1));
    prime2AFC_second(block) = sqrt(2) * norminv((tot_correct_second + 0.5) / (n_second + 1));
end

primeGONOGO = [primeGONOGO_full,primeGONOGO_first,primeGONOGO_second];
primeGONOGO = array2table(primeGONOGO,"VariableNames",{'GOfull','GOfirst','GOsecond'});

prime2AFC = [prime2AFC_full,prime2AFC_first,prime2AFC_second];
prime2AFC = array2table(prime2AFC,'VariableNames',{'2AFCfull','2AFCfirst','2AFCsecond'});

prime_table = horzcat(primeGONOGO,prime2AFC);
    

