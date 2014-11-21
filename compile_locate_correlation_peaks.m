function compile_locate_correlation_peaks();

% Run coder to generate the mex file.
codegen -config coder.config('mex') locate_correlation_peaks -args {coder.typeof(1.00, [inf, inf])};

end
