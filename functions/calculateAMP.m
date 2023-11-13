function [amplitude] = calculateAMP(stream)

% [pks,locs] = findpeaks(stream,ts1);
% idx = (locs>0) & (locs<0.75);
% amplitude = max(pks(idx));
% if isempty(amplitude)
%     amplitude = NaN;
% end
amplitude = max(stream);
end