clear
BLOCKPATH = '/Users/brandon/My Drive/prl/dual_fiber/tanks/74MUL_Acq5_NA_Empty_NA_NA';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
data = prl_df_epocs(data,1);

cue = data.epocs.St1_.onset;
cRew = data.epocs.cRewA.onset;
cNoRew = data.epocs.cNoRewA.onset;
iRew = data.epocs.iRewA.onset;
iNoRew = data.epocs.iNoRewA.onset;

[session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(...
    cue, cRew, cNoRew, iRew, iNoRew);
[data, errorProbLeverTS, errorProbCueTS] = errorProbExtract(data, session_identifiers, 4, 2);