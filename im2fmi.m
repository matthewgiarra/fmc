function FMI = im2fmi(IMAGE, SPECTRUM_HEIGHT, SPECTRUM_WIDTH, X_LP, Y_LP, COMPILED)
% Calculate the Fourier-Mellin Invariant Descriptor of an image (or other rectilinear data)
% INPUTS:
%   IMAGE = 2D grayscale image to be converted to FMI
%   SPECTRUM_HEIGHT = height (pixels) of the FFT magnitude of the IMAGE
%   SPECTRUM_WIDTH = width (pixels) of the FFT magnitude of the IMAGE
%   X_LP = Non-integer column coordinates at which the FFT magnitude is resampled
%   X_LP = Non-integer row coordinates at which the FFT magnitude is resampled
%   COMPILED = flag specifying whether or not to use the compiled interp2
%   code "ba_interp2". 
%
% OUTPUTS:
%   FMI = FMI descriptor of the image (2D grayscale)

% Case of using compiled code.
if COMPILED
    FMI = ba_interp2(fftshift(abs(fftn(IMAGE, [SPECTRUM_HEIGHT, SPECTRUM_WIDTH]))), X_LP, Y_LP, 'cubic');
else
    % Case of using matlab's built-in interp2 (which is slower than ba_interp2)
    FMI = interp2(fftshift(abs(fftn(IMAGE, [SPECTRUM_HEIGHT, SPECTRUM_WIDTH]))), X_LP, Y_LP, 'bicubic', 0);
end

end
