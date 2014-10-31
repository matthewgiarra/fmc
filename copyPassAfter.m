function JOBFILE = copyPassAfter(JOBFILE, PASS)
% This function updates a jobfile by copying a pass and repeating it

% Original number of passes specified
nPasses = length(JOBFILE.Parameters.Processing);

% If the variable "PASS" was specified as the string "end",
% then copy the new pass at the end of the old jobfile.
if ischar(PASS) && strcmpi(PASS, 'end');
    PASS = nPasses;
end

% Copy all the passes after PASS and shift them down by one
JOBFILE.Parameters.Processing(PASS+1:end+1) = JOBFILE.Parameters.Processing(PASS:end);



end