clearvars -except input;
correct = input(:,1);
rewards = input(:,2);
n = size(correct,1);

% Determine number of columns needed
num_blocks = ceil(n/30);

% Preallocate
cummulative_cor = nan(30, num_blocks);
cummulative_rew = nan(30, num_blocks);

for block = 1:num_blocks
    cor_count = 0;
    rew_count = 0;
    
    row_start = (block-1)*30 + 1;
    row_end = min(block*30, n);
    
    for i = row_start:row_end
        local_row = i - row_start + 1;
        
        if correct(i) == 1
            cor_count = cor_count + 1;
        end
        cummulative_cor(local_row, block) = cor_count;
        
        if rewards(i) == 1
            rew_count = rew_count + 1;
        end
        cummulative_rew(local_row, block) = rew_count;
    end
end
