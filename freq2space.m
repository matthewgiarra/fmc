function SPATIAL_PLANE = freq2space(SPECTRAL_PLANE)
% freq2space(SPECTRUM, SPECTRUM, SPECTRUM_WIDTH) 
% Shifts the zero-frequency component of a correlation in the spectral domain 
% to the center of spectrum, then transforms the correlation into the
% spatial domain.
%
% INPUTS:
%   SPECTRAL_PLANE = Complex matrix with SPECTRUM_HEIGHT rows and SPECTRUM_WIDTH
%       columns. The coordinates of this matrix represent frequencies or
%       wave numbers.
%   
% OUTPUTS:
%   SPATIAL_PLANE = Real matrix of the same dimensions as SPECTRAL_PLANE
%       containing the inverse Fourier transform of SPECTRAL_PLANE after 
%       shifting the zero-frequency component to the center of the matrix.
%
% SEE ALSO:
%   fftshift

% Determine the number of rows and columns of the input array.
[spectrum_height, spectrum_width] = size(SPECTRAL_PLANE);

% Inverse fft of the correlation plane, then take the absolute value
% of its real part.
SPATIAL_PLANE = fftshift(abs(real(ifftn(SPECTRAL_PLANE, [spectrum_height, spectrum_width], 'symmetric'))));

end