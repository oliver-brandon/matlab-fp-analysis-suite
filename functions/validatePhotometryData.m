function [ok, problems] = validatePhotometryData(data, requiredStreams, requiredEpocs, minDurationSec)
% validatePhotometryData Validate common stream/epoc requirements.

if nargin < 2 || isempty(requiredStreams)
    requiredStreams = {};
end
if nargin < 3 || isempty(requiredEpocs)
    requiredEpocs = {};
end
if nargin < 4 || isempty(minDurationSec)
    minDurationSec = 0;
end

problems = {};

if ~isstruct(data)
    problems{end+1,1} = 'data is not a struct';
    ok = false;
    return
end

if ~isfield(data, 'streams') || ~isstruct(data.streams)
    problems{end+1,1} = 'data.streams is missing';
else
    for i = 1:numel(requiredStreams)
        streamName = requiredStreams{i};
        if ~isfield(data.streams, streamName)
            problems{end+1,1} = sprintf('missing stream %s', streamName);
            continue
        end
        stream = data.streams.(streamName);
        if ~isfield(stream, 'data') || isempty(stream.data)
            problems{end+1,1} = sprintf('stream %s has no data', streamName);
        end
        if ~isfield(stream, 'fs') || isempty(stream.fs) || stream.fs <= 0
            problems{end+1,1} = sprintf('stream %s has invalid fs', streamName);
        elseif minDurationSec > 0 && isfield(stream, 'data')
            durationSec = numel(stream.data) / stream.fs;
            if durationSec < minDurationSec
                problems{end+1,1} = sprintf( ...
                    'stream %s shorter than %.3g seconds', streamName, minDurationSec);
            end
        end
    end
end

if ~isfield(data, 'epocs') || ~isstruct(data.epocs)
    problems{end+1,1} = 'data.epocs is missing';
else
    for i = 1:numel(requiredEpocs)
        epocName = requiredEpocs{i};
        if ~isfield(data.epocs, epocName)
            problems{end+1,1} = sprintf('missing epoc %s', epocName);
        elseif ~isfield(data.epocs.(epocName), 'onset')
            problems{end+1,1} = sprintf('epoc %s has no onset', epocName);
        end
    end
end

ok = isempty(problems);

end
