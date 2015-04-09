
% Load the particle image MAT files
% Fill this out
path_to_mat_file = '';

% Load the mat file
load(path_to_mat_file);

% Measure region dimensions
[height, width, nImages] = size(region_mat_01);

% Allocate the correlation plane.
correlation_plane = zeros(height, width);

% Make a Gaussian window
g = gaussianWindowFilter([height, width], [0.5, 0.5], 'fraction');

% Loop over all the regions
for k = 1 : 1000
    
    % Window the regions and convert them to doubles.
    region_01 = g .* double(imageMatrix1(:, :, k));
    region_02 = g .* double(imageMatrix2(:, :, k));
    
    % Calculate the phase-only correlation
    % and add it to the ensemble.
    correlation_plane = correlation_plane + ...
        freq2space(fftshift(splitComplex(fftn(double(region_02), ...
        [height, width]) .* conj(fftn(double(region_01), ...
        [height, width])))));

end

% This is how many indices you want to skip for plotting,
% since the mesh plot will be really dense otherwise
% and might not look nice.
Skip = 4;

% Skip some data points.
c = correlation_plane(1 : Skip : end, 1 : Skip : end);

% Normalize the correlation plane to its max value
% and plot it
mesh(c ./ max(c(:)), 'edgecolor', 'black'); 

% Set an axis limit.
zlim([0, 1.1]);

% Remove axis ticks
set(gca, 'xticklabel', '');
set(gca, 'yticklabel', '');
set(gca, 'zticklabel', '');

% Make the axes square
axis square;

% Add a thick box border.
box on



