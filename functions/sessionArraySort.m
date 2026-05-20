function [session_identifiers,lever_session_ts,trial_number,trial_name] = sessionArraySort(...
    CUETTL,TTL1,TTL2,TTL3,TTL4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extracts cue/lever TTL events into paired cue -> lever rows.
    %
    % Output session_identifiers is an Nx2 array:
    %   column 1 = timestamp
    %   column 2 = event type (0 cue, 1 cRew, 2 cNoRew, 3 iRew, 4 iNoRew)
    %
    % Only levers whose immediately preceding event is a cue are retained.
    % This avoids positional deletion rules that can silently shift pairings.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cueTs = normalizeTimestamps(CUETTL);
    ttl1Ts = normalizeTimestamps(TTL1);
    ttl2Ts = normalizeTimestamps(TTL2);
    ttl3Ts = normalizeTimestamps(TTL3);
    ttl4Ts = normalizeTimestamps(TTL4);
    leverTs = [
        ttl1Ts, ones(numel(ttl1Ts),1);
        ttl2Ts, ones(numel(ttl2Ts),1)*2;
        ttl3Ts, ones(numel(ttl3Ts),1)*3;
        ttl4Ts, ones(numel(ttl4Ts),1)*4
    ];

    cueEvents = [cueTs, zeros(numel(cueTs),1)];
    allEvents = sortrows([cueEvents; leverTs], 1);

    session_identifiers = [];
    for i = 2:size(allEvents,1)
        isLever = allEvents(i,2) > 0;
        previousIsCue = allEvents(i-1,2) == 0;
        if isLever && previousIsCue
            session_identifiers = [session_identifiers; allEvents(i-1,:); allEvents(i,:)]; %#ok<AGROW>
        end
    end

    if isempty(session_identifiers)
        lever_session_ts = [];
        trial_number = [];
        trial_name = {};
        return
    end

    lever_session_ts = session_identifiers(:,1);
    trial_number = session_identifiers(:,2);

    trial_outcome = trial_number(2:2:end,1);
    trial_name = cell(size(trial_outcome));
    for i = 1:length(trial_outcome)
        switch trial_outcome(i)
            case 1
                trial_name{i} = ['C+ Trial ' num2str(i)];
            case 2
                trial_name{i} = ['C- Trial ' num2str(i)];
            case 3
                trial_name{i} = ['I+ Trial ' num2str(i)];
            case 4
                trial_name{i} = ['I- Trial ' num2str(i)];
        end
    end

end

function ts = normalizeTimestamps(ts)
if isempty(ts) || isequal(ts, 0)
    ts = [];
else
    ts = ts(:);
    ts = ts(~isnan(ts) & ts ~= 0);
end
end
