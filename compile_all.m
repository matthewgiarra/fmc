function compile_all()
% This function compiles all the compileable codes 
% in the FMC library to MEX files.

% Get the current path
parent_dir = pwd;

% Navigate to the ba_interpolation directory (can't just call the
% code's path directly I guess)
cd ba_interpolation;

% Inform the user
fprintf('Compiling ba_interp2.m to MEX...\n');

% Compile the interpolation code.
compile_ba_interp2;

% Navigate back up to the parent directory
cd(parent_dir);

% Inform the user
fprintf('\nCompiling locate_correlation_peaks.m to MEX...\n');

% Compile the correlation peak locating function
compile_locate_correlation_peaks;

end