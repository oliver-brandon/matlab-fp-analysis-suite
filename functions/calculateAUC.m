function [auc] = calculateAUC(stream,ts1)


auc = trapz(ts1,stream);


% Calculate AUC above x=0 %
% positive_indices = stream > 0;
% if sum(positive_indices) < 3
%     auc = NaN;
% else
%     y_pos = stream(1,positive_indices);
%     x_pos = ts1(1,positive_indices);
%     auc = trapz(x_pos,y_pos);
% end


end