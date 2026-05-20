function y = blockMeanDownsample(x, factor)
% blockMeanDownsample Downsample rows by averaging non-overlapping blocks.

nBins = floor(size(x,2) / factor);
if nBins == 0
    y = zeros(size(x,1), 0);
    return
end

x = x(:,1:nBins*factor);
x = reshape(x, size(x,1), factor, nBins);
y = squeeze(mean(x, 2, 'omitnan'));
if isvector(y)
    y = reshape(y, size(x,1), nBins);
end

end
