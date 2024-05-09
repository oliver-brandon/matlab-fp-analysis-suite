close all;
clear smooth_sigs1 smooth_sigs2;
adjust_only = 0;
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
        smooth_sigs1(:,i) = smoothdata(sigs(:,i),'movmean',25);
        smooth_sigs2(:,i) = smoothdata(sigs(:,i),'movmean',75);
    end
end


for j = 1:size(smooth_sigs1,2)
    % adjusts streams to baseline of zero at -0.5s %
    if isnan(smooth_sigs1(1,j))
        continue
    elseif smooth_sigs1(1,j) < 0
        val = smooth_sigs1(1,j);
        diff = 0 - val;
        smooth_sigs1(1:end,j) = smooth_sigs1(1:end,j) + abs(diff);
    elseif smooth_sigs1(1,j) > 0
        val = smooth_sigs1(1,j);
        diff = 0 - val;
        smooth_sigs1(1:end,j) = smooth_sigs1(1:end,j) - abs(diff);
    end    
end
for j = 1:size(smooth_sigs2,2)
    % adjusts streams to baseline of zero at -0.5s %
    if isnan(smooth_sigs2(1,j))
        continue
    elseif smooth_sigs2(1,j) < 0
        val = smooth_sigs2(1,j);
        diff = 0 - val;
        smooth_sigs2(1:end,j) = smooth_sigs2(1:end,j) + abs(diff);
    elseif smooth_sigs2(1,j) > 0
        val = smooth_sigs2(1,j);
        diff = 0 - val;
        smooth_sigs2(1:end,j) = smooth_sigs2(1:end,j) - abs(diff);
    end    
end
subplot(3,1,1)
plot(sigs)
subplot(3,1,2)
plot(smooth_sigs1)
subplot(3,1,3)
plot(smooth_sigs2)