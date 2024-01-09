close all;
clear smooth_sigs1 smooth_sigs2;
for i = 1:size(sigs,2)
    if sigs(1,i) == 0
        smooth_sigs1(:,i) = NaN;
        smooth_sigs2(:,i) = NaN;
        continue
    end
    smooth_sigs1(:,i) = smoothdata(sigs(:,i),'movmean',10);
    smooth_sigs2(:,i) = smoothdata(sigs(:,i),'movmean',15);
end
subplot(3,1,1)
plot(sigs)
subplot(3,1,2)
plot(smooth_sigs1)
subplot(3,1,3)
plot(smooth_sigs2)