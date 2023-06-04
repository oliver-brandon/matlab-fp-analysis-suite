filename = "/Users/brandon/DA_PRL/PRL_GRABDA/Acq1/Vehicle/01M_Acq1_NA.mat";
load(filename)

if isfield(data.streams, 'x405A')
    cue = data.epocs.St1_.onset;
    cRew = data.epocs.cRewA.onset;
    cNoRew = data.epocs.cNoRewA.onset;
    iRew = data.epocs.iRewA.onset;
    iNoRew = data.epocs.iNoRewA.onset;
elseif isfield(data.sreams, 'x405C')
    cue = data.epocs.St2_.onset;
    cRew = data.epocs.cRewC.onset;
    cNoRew = data.epocs.cNoRewC.onset;
    iRew = data.epocs.iRewC.onset;
    iNoRew = data.epocs.iNoRewC.onset;
end


[session_identifiers,lever_session_ts,trial_number,trial_name,lever_latency] = sessionArraySort(...
    cue,cRew,cNoRew,iRew,iNoRew);

session_identifiers = vertcat(session_identifiers(1:2,:),session_identifiers(4:end,:));
cue_cRew = [];
cue_cNoRew = [];
cue_iRew = [];
cue_iNoRew = [];

for i = 1:height(session_identifiers)
    if session_identifiers(i,2) == 0
        continue
    elseif session_identifiers(i,2) == 1
        cue_cRew = [cue_cRew; session_identifiers(i-1:i,1:2)];
    elseif session_identifiers(i,2) == 2
        cue_cNoRew = [cue_cNoRew; session_identifiers(i-1:i,1:2)];
    elseif session_identifiers(i,2) == 3
        cue_iRew = [cue_iRew; session_identifiers(i-1:i,1:2)];
    elseif session_identifiers(i,2) == 4
        cue_iNoRew = [cue_iNoRew; session_identifiers(i-1:i,1:2)];
    end 
end
if isempty(cue_cRew)
    cue_cRew = 0;
elseif isempty(cue_cNoRew)
    cue_cNoRew = 0;
elseif isempty(cue_iRew)
    cue_iRew = 0;
elseif isempty(cue_iNoRew)
    cue_iNoRew = 0;
end

data.epocs.cue2lever.cRew = cue_cRew;
data.epocs.cue2lever.cNoRew = cue_cNoRew;
data.epocs.cue2lever.iRew = cue_iRew;
data.epocs.cue2lever.iNoRew = cue_iNoRew;
