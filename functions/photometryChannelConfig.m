function cfg = photometryChannelConfig(channel)
% photometryChannelConfig Return stream and TTL field names for one channel.

if isnumeric(channel)
    if channel == 1
        channel = 'A';
    elseif channel == 2
        channel = 'C';
    else
        error('photometryChannelConfig:UnknownChannel', ...
            'Numeric channel must be 1 (A) or 2 (C).');
    end
end

channel = upper(string(channel));

switch channel
    case "A"
        cfg.channel = 'A';
        cfg.isos = 'x405A';
        cfg.signal = 'x465A';
        cfg.altSignal = 'x560A';
        cfg.cue = 'St1_';
        cfg.correct = 'CL1_';
        cfg.incorrect = 'IL1_';
        cfg.pellet = 'Pe1_';
        cfg.cRew = 'cRewA';
        cfg.cNoRew = 'cNoRewA';
        cfg.iRew = 'iRewA';
        cfg.iNoRew = 'iNoRewA';
        cfg.cueCor = 'cueCorA';
        cfg.cueInc = 'cueIncA';
    case "C"
        cfg.channel = 'C';
        cfg.isos = 'x405C';
        cfg.signal = 'x465C';
        cfg.altSignal = 'x560C';
        cfg.cue = 'St2_';
        cfg.correct = 'CL2_';
        cfg.incorrect = 'IL2_';
        cfg.pellet = 'Pe2_';
        cfg.cRew = 'cRewC';
        cfg.cNoRew = 'cNoRewC';
        cfg.iRew = 'iRewC';
        cfg.iNoRew = 'iNoRewC';
        cfg.cueCor = 'cueCorC';
        cfg.cueInc = 'cueIncC';
    otherwise
        error('photometryChannelConfig:UnknownChannel', ...
            'Channel must be A/C or 1/2.');
end

end
