Instructions for error probability stream extracting using the David Root script

1. To ensure all scripts/functions are up to date, download the entire matlab-fp-analysis-suite from https://github.com/oliver-brandon/matlab-fp-analysis-suite

2. Add matlab-fp-analysis-suite/functions to your MATLAB path

3. In the prl-analysis folder, open EpochPlot_HeatMap_PRL_errors.m

4. Copy/paste the path to a tank in the 'BLOCKPATH' variable

5. Choose which stream store to plot/extract

6. Edit 'errorType' for the desired error probability to extract

7. Edit 'lever' to change which lever the errors are based on

8. Edit 'TTL' to choose which operant box TTLs should be used

9. Edit 'REF_EPOC' to include the desired error probability epoc to extract. The available epocs depend on the errorType being extracted. For example, if you extract winStay from the correct lever, the following epocs will be created: WIN_stay_LEV, WIN_stay_CUE, win_STAY_LEV, and win_STAY_CUE. The words in all caps indicate which part of the error to extract. For example, WIN_stay_CUE will extract the cue preceding a 'WIN' and win_STAY_LEV will extract the 'STAY' lever press.
	List of all possible epochs:
		win-stay: WIN_stay_LEV, WIN_stay_CUE, win_STAY_LEV, win_STAY_CUE
		win-shift: WIN_shift_LEV, WIN_shift_CUE, win_SHIFT_LEV, win_SHIFT_CUE
		lose-stay: LOS_stay_LEV, LOS_stay_CUE, los_STAY_LEV, los_STAY_CUE
		lose-shift: LOS_shift_LEV, LOS_shift_CUE, los_SHIFT_LEV, los_SHIFT_CUE

10. The script will plot the averaged epoc event for that tank included in the variable 'meanSignal1' and will also include each individual trial stream for that epoc in the 'zall' variable

