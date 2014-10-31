function compile_measurePeakHeightRatio;

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Set the dynamic variable sizes
coder.varsize('correlationPlane', [inf, inf], [true, true]);

% Run coder to generate the mex file.
codegen -config cfg measurePeakHeightRatio -args {correlationPlane};


end