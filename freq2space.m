function [spatialCorr] = freq2space(SPECTRUM, SPECTRUM_HEIGHT, SPECTRUM_WIDTH)
% freq2space(SPECTRUM, SPECTRUM, SPECTRUM_WIDTH) 
% Shifts the zero-frequency component of a correlation in the spectral domain 
% to the center of spectrum, then transforms the correlation into the
% spatial domain.
%
% SEE ALSO:
%   fftshift

% Inverse fft of the correlation plane, then take the absolute value
% of its real part.
spatialCorr = fftshift(abs(real(ifftn(SPECTRUM, [SPECTRUM_HEIGHT, SPECTRUM_WIDTH], 'symmetric'))));

end