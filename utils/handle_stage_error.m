function handle_stage_error(ME, stage_num, stop_on_error)
%HANDLE_STAGE_ERROR Provides formatted error reporting for pipeline stages.
%   This function prints a detailed error message to the standard error
%   stream, including the stage number, error message, and the file/line
%   where the error occurred. It then rethrows the error if 'stop_on_error'
%   is true, halting execution.
%
%   Syntax:
%       handle_stage_error(ME, stage_num, stop_on_error)
%
%   Inputs:
%       ME            - The MException object caught in a try-catch block.
%       stage_num     - The number of the stage that failed.
%       stop_on_error - Logical. If true, rethrows the exception to halt the script.

    arguments
        ME (1,1) MException
        stage_num (1,1) double
        stop_on_error (1,1) logical
    end

    fprintf(2, '[X] Stage %.1f failed.\n', stage_num); % Use %.1f to support stages like 7.5
    fprintf(2, '    Error: %s\n', ME.message);

    % Provide file and line number for easier debugging
    if ~isempty(ME.stack)
        file_info = ME.stack(1);
        fprintf(2, '    In file: %s (line %d)\n', file_info.file, file_info.line);
    end

    if stop_on_error
        rethrow(ME);
    end
end