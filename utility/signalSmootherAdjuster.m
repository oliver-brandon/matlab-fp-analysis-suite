% Paste into command window: sigs = [];
% Paste signals into sigs variable (long format, not wide)
% Adjusted/smoothed signals will be in smooth_sigs1 or smooth_sigs2

close all;
clear smooth_sigs1 smooth_sigs2;
adjust_only = 1; % 0 = smooth + adjust, 1 = adjust only
zeroto = 0; % time in seconds to zero signal to
plotSigs = 0; % 1 = yes, 0 = no

if adjust_only == 1
    smooth_sigs1 = sigs;
    smooth_sigs2 = sigs;
else
    for i = 1:size(sigs,2)
        if sigs(:,i) == 0
            smooth_sigs1(:,i) = NaN;
            smooth_sigs2(:,i) = NaN;
            continue
        end
        smooth_sigs1(:,i) = smoothdata(sigs(:,i),'movmean',20);
        smooth_sigs2(:,i) = smoothdata(sigs(:,i),'movmean',60);
    end
end

ts1 = -2 + (1:length(sigs)) / 1017*10;
idx = find(ts1>zeroto,1);
for j = 1:size(smooth_sigs1,2)
    % adjusts streams to baseline of zero at -0.5s %
    if isnan(smooth_sigs1(idx,j))
        continue
    elseif smooth_sigs1(idx,j) < 0
        val = smooth_sigs1(idx,j);
        diff = 0 - val;
        smooth_sigs1(1:end,j) = smooth_sigs1(1:end,j) + abs(diff);
    elseif smooth_sigs1(idx,j) > 0
        val = smooth_sigs1(idx,j);
        diff = 0 - val;
        smooth_sigs1(1:end,j) = smooth_sigs1(1:end,j) - abs(diff);
    end    
end
for j = 1:size(smooth_sigs2,2)
    % adjusts streams to baseline of zero at -0.5s %
    if isnan(smooth_sigs2(idx,j))
        continue
    elseif smooth_sigs2(idx,j) < 0
        val = smooth_sigs2(idx,j);
        diff = 0 - val;
        smooth_sigs2(1:end,j) = smooth_sigs2(1:end,j) + abs(diff);
    elseif smooth_sigs2(idx,j) > 0
        val = smooth_sigs2(idx,j);
        diff = 0 - val;
        smooth_sigs2(1:end,j) = smooth_sigs2(1:end,j) - abs(diff);
    end    
end
if plotSigs == 1
    subplot(3,1,1)
    plot(sigs)
    subplot(3,1,2)
    plot(smooth_sigs1)
    subplot(3,1,3)
    plot(smooth_sigs2)
else
    disp('')
end
sigs = [];