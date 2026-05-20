function appendStructField(filename, prl_phase, treatment, fieldData)
% APPENDSTRUCTFIELD Appends a new field to an existing structure .mat file
%   APPENDSTRUCTFIELD(FILENAME, FIELDNAME, FIELDDATA) appends a new field 
%   with name FIELDNAME and data FIELDDATA to the structure stored in the
%   .mat file specified by FILENAME in the current working directory.

% Load the existing structure from the .mat file
S = load(filename, 'prl_ERT');
if ~isfield(S, 'prl_ERT')
    error('appendStructField:MissingStruct', ...
        '%s does not contain a prl_ERT variable.', filename);
end
prl_ERT = S.prl_ERT;

% Append the new field to the structure
prl_ERT.(prl_phase).(treatment) = fieldData;

% Save the updated structure back to the .mat file
save(filename, 'prl_ERT', '-append');
end
