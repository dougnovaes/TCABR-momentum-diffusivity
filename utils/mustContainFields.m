% Salvar como: utils/mustContainFields.m

function mustContainFields(s, fields)
%MUSTCONTAINFIELDS Custom validation function for arguments blocks.
%   Checks if the input struct 's' contains all field names specified in
%   the cell array 'fields'. Throws an error if any field is missing.

    if ~isstruct(s)
        eid = 'Validation:notStruct';
        msg = 'Input must be a struct.';
        throwAsCaller(MException(eid, msg));
    end
    
    for i = 1:numel(fields)
        if ~isfield(s, fields{i})
            eid = 'Validation:missingField';
            msg = sprintf('Input struct is missing required field: ''%s''', fields{i});
            throwAsCaller(MException(eid, msg));
        end
    end
end