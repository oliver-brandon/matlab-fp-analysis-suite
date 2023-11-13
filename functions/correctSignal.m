function [F465_corrected] = correctSignal(F405,F465)

bls = polyfit(F465(1:end), F405(1:end),1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_df_all = F465 - Y_fit_all;

F465_corrected = Y_df_all;

end