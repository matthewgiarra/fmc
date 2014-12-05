function [MUTUAL_INFORMATION, MAX_INTENSITY_MEAN, AUTO_CORR_PEAK_DIAMETER_MEAN, CROSS_CORRELATION_REAL] = calculate_mutual_information(REGION_01, REGION_02, COMPILED)
% [MUTUAL_INFORMATION, MAX_INTENSITY_MEAN, AUTO_CORR_PEAK_DIAMETER_MEAN, CROSS_CORRELATION_REAL] = calculate_mutual_information(REGION_01, REGION_02, COMPILED)
% This function is used to calculates the "mutual information" between two
% PIV interrogation regions (IRs). Mutual information is defined here 
% as the ratio of the maximum value of the inverse Fourier transform of the 
% cross correlation between the pair of IRs (after minimum-subtraction) 
% to the maximum value of the auto-correlation of a "standard particle" 
% image. The "standard particle" image is a Gaussian function whose
% maximum value is the average of the maximum intensities of the two
% interrogation regions, and whose "diameter" (defined as four times the 
% standard deviation of the Gaussian function that best fits the
% correlation peak) is the same as the average diameter of the primary
% peaks of the auto-correlations of the two IRs.
%
% INPUTS
%   REGION_01 = [M x N] Matrix containing the intensity data of the first
%       interrogation region, where M and N are the number of rows and columns
%       in the interrogation region. 
%
%   REGION_02 = [M x N] Matrix containing the intensity data of the second
%       interrogation region, where M and N are the number of rows and columns
%       in the interrogation region. 
% 
%   COMPILED = Boolean flag specifying whether or not to run compiled codes
%       within SUBPIXEL. The default value is zero (i.e., do not run
%       compiled codes).
%
% OUTPUTS:
%   MUTUAL_INFORMATION = Mutual information between the two interrogation
%       regions, where "mutual information" is defined in the header of this
%       file.
%
%   MAX_INTENSITY_MEAN = Average of the maximum intensities of the two
%       interrogation regions (scalar).
%
%   AUTO_CORR_PEAK_DIAMETER_MEAN = Average of the diameters of the
%       peaks of the auto correlations of each interrogation region, where the
%       diameter is defined as four times the standard deviation.
%
%   CROSS_CORRELATION_REAL = [M x N] Matrix containing the inverse 
%       Fourier Transform of the full cross correlation plane, with the
%       zero-frequency coordinate located at the center of the matrix.
%
% SEE ALSO:
%   SUBPIXEL, FREQ2SPACE

% Default to not using compiled codes
if nargin < 3
    COMPILED = 0;
end

% Size of the interrogation regions.
[region_height, region_width] = size(REGION_01);

% Particle intensity as the maximum intensity in the particle image
max_intensity_region_01 = max( double(REGION_01(:)) );
max_intensity_region_02 = max( double(REGION_02(:)) );

% Average of the max intensities of the two interrogation regions.
% This is used to specify the intensity of the standard Gaussian particle.
MAX_INTENSITY_MEAN = mean( [max_intensity_region_01, max_intensity_region_02] );

% Calculate the complex spectral Fourier transforms
% of each interogation region
fourier_transform_01 = fftn( double(REGION_01), [region_height region_width] );
fourier_transform_02 = fftn( double(REGION_02), [region_height region_width] );

% This is the complex cross correlation of the two images, calculated by conjugate
% multiplying their complex Fourier transforms.
cross_correlation_complex = (fourier_transform_02 .* conj(fourier_transform_01));

% This is the inverse Fourier transform of the complex cross correlation.
% This is capitalized because it is an output of the function.
CROSS_CORRELATION_REAL = freq2space(cross_correlation_complex);

% This is the height of the minimum-subtracted correlation plane.
cross_correlation_min_subtracted_max = max( CROSS_CORRELATION_REAL(:) ) - min( CROSS_CORRELATION_REAL(:) );

% Calucluate the complex spectral autocorrelations
% of the two interrogation regions.
auto_corr_spectral_01 = fourier_transform_01 .* conj(fourier_transform_01); 
auto_corr_spectral_02 = fourier_transform_02 .* conj(fourier_transform_02);

% Inverse Fourier transforms of the
% complex spectral autocorrelations, 
% followed by FFTSHIFT coordinate shift
auto_corr_spatial_01 = freq2space(auto_corr_spectral_01);                
auto_corr_spatial_02 = freq2space(auto_corr_spectral_02);   

% minimum subtraction to eliminate background noise effect
auto_corr_min_sub_01 = auto_corr_spatial_01 - min(auto_corr_spatial_01(:));
auto_corr_min_sub_02 = auto_corr_spatial_02 - min(auto_corr_spatial_02(:));

% Calculate the diameters of the primary peak of the autocorrelations.
% The function subpixel is way overkill for this one operation. 
% To do: split subpixel into separate subroutines.
[~, ~, ~, auto_corr_peak_diameter_01, ~] = subpixel(auto_corr_min_sub_01, ones(size(auto_corr_min_sub_01)), 1, 0, COMPILED);
[~, ~, ~, auto_corr_peak_diameter_02, ~] = subpixel(auto_corr_min_sub_02, ones(size(auto_corr_min_sub_02)), 1, 0, COMPILED);

% Calculate the average of the diameters of the two autocorrelation peaks.
AUTO_CORR_PEAK_DIAMETER_MEAN = mean([auto_corr_peak_diameter_01, auto_corr_peak_diameter_02]);

% Grid of coordinates over which to calculate the standard particle.
% The speed of this code could be improved by removing meshgrid
% and passing coordinate matrices as inputs to the function. 
% Meshgrid is fairly slow.
[x, y] = meshgrid(1 : region_width, 1 : region_height);

% Calculate the integer nearest to the center of the coordinate matrix
xc = round(region_width  / 2);
yc = round(region_height / 2);

% This makes a Gaussian distribution cenetered at the geometric centroid of
% the image, with max intensity of MAX_INTENSITY_MEAN and standard
% deviation of AUTO_CORR_PEAK_DIAMETER_MEAN
% This is the "standard particle"
standard_particle = MAX_INTENSITY_MEAN * exp(-16*((x - xc).^2 + (y - yc).^2) / AUTO_CORR_PEAK_DIAMETER_MEAN ^ 2);

% Calculate the max value of the auto correlation of the standard particle.
standard_particle_auto_corr_max = sum(standard_particle(:) .^ 2);

% Mutual Information is the ratio between the contribution of all correlated particle and the contribution of one particle
MUTUAL_INFORMATION = cross_correlation_min_subtracted_max / standard_particle_auto_corr_max; 

end

