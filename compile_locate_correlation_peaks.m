function compile_locate_correlation_peaks();

% Example variables
spatial_correlation_plane = coder.typeof(1.00, [inf, inf]);


% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Run coder to generate the mex file.
codegen -config cfg locate_correlation_peaks -args {spatial_correlation_plane};

end