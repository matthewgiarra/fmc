function [TY, TX, TH, S, FMC_PEAK_HEIGHT_RATIO, RPC_PEAK_HEIGHT_RATIO, RPC_PEAK_HEIGHT, RPC_PEAK_DIAMETER] = FMC(REGION1, REGION2, SPATIALWINDOW, IMAGESPECTRALFILTER, FMI_WINDOW, FMI_SPECTRAL_FILTER, XREGION, YREGION, SPECTRUM_HEIGHT, SPECTRUM_WIDTH, XLP, YLP, RMIN, RMAX, DIFFERENCEMETHOD, COMPILED)
% This function estimates the horizontal and vertial translations, angle of
% rotation, and isotropic scaling factor that relate a pair of images.
% For PIV, this function should be run on pairs of interrogation regions.
%
% Need to add header...

% Determine the height and width of the regions (in pixels)
[regionHeight, regionWidth] = size(REGION1);

% Convert images to their FMI descriptors
[FMI_01] = im2fmi(SPATIALWINDOW .* REGION1, SPECTRUM_HEIGHT, SPECTRUM_WIDTH, XLP, YLP, COMPILED);
[FMI_02] = im2fmi(SPATIALWINDOW .* REGION2, SPECTRUM_HEIGHT, SPECTRUM_WIDTH, XLP, YLP, COMPILED);

% Height and width of the FMI descriptors.
[fmi_height, fmi_width] = size(FMI_01);

% Perform the cross phase-only correlation of the two FMI descriptors.
spectralRPC = fftshift(splitComplex(fftn(FMI_WINDOW .* FMI_02, [fmi_height, fmi_width]) .* conj(fftn(FMI_WINDOW .* FMI_01, [fmi_height, fmi_width])))) .* FMI_SPECTRAL_FILTER;

% Convert the phase correlation plane of the two input images from the spectral domain to the spatial domain
fmi_correlation_plane = freq2space(spectralRPC);

% Deal with +/- 90 degree rotations where the correlation peak can overlap the 
% vertical edges of the FMC plane.
[~, maxLoc] = max(fmi_correlation_plane(:));

% Determine row position of peak.
% This is eqivalent to but faster than sub2ind
r = rem(maxLoc - 1, fmi_height) + 1;

% Automatically output results for multiple peaks from the FMC plane.
multi_peak_flag = 1;

% Shift the FFT correlation plane so that the peak is in the center of the plane.
% Double if statement to avoid logical OR (which is slow)
if r <= 9
   
    % Shift the plane so the peak is in the center.
    fmi_correlation_plane = fftshift(fmi_correlation_plane, 1);
    % Subpixel peak finding
    [fy, fmi_peak_shift_x, fmc_peak_heights] = subpixel(fmi_correlation_plane, ones([fmi_height, fmi_width]), 1, multi_peak_flag);
    fmi_peak_shift_y = fy + (fmi_height / 2);
elseif r >= (fmi_height - 8)
    
    % Shift the plane so the peak is in the center
    fmi_correlation_plane = fftshift(fmi_correlation_plane, 1);
    % Subpixel peak finding
    [fy, fmi_peak_shift_x, fmc_peak_heights] = subpixel(fmi_correlation_plane, ones([fmi_height, fmi_width]), 1, multi_peak_flag);
    fmi_peak_shift_y = fy + (fmi_height / 2);
else
    
    % Subpixel peak finding
    [fmi_peak_shift_y, fmi_peak_shift_x, fmc_peak_heights] = ...
        subpixel(fmi_correlation_plane, ones([fmi_height, fmi_width]), 1, multi_peak_flag);
end

% Measure the ratio of the largest to second largest correlation peaks.
FMC_PEAK_HEIGHT_RATIO = fmc_peak_heights(1) / fmc_peak_heights(2);

% Angular shift in radians. In image coordinates, a positive rotation angle
% indicates a clockwise rotation. This is because the positive direction of the vertical axis is downward in image coordinates. 
estimated_rotation_initial = -1 * fmi_peak_shift_y * 2 * pi * (1 - 1 / fmi_height) / fmi_height;

% Scaling factor
estimated_scaling_initial = exp(-1 * log(RMAX/RMIN) * fmi_peak_shift_x / (fmi_width - 1));

% Calculate the horizontal and vertical coordinate of the
% geometric centroid of the region.
xc = regionWidth / 2 - 0.5;
yc = regionHeight / 2 - 0.5;

% Reset the multi_peak flag depending on the peak height ratio 
% of the FMC plane. Use multiple peaks for peak height ratios
% lower than 1.5 (arbitrary)
multi_peak_flag = FMC_PEAK_HEIGHT_RATIO < 1.5;

% Resolve the rotation ambiguity
[TH, S, TX, TY, RPC_PEAK_HEIGHT_RATIO, RPC_PEAK_HEIGHT, RPC_PEAK_DIAMETER] ...
    = resolveRotationAmbiguity(REGION1, REGION2, ...
    SPATIALWINDOW, IMAGESPECTRALFILTER, ...
    estimated_rotation_initial, estimated_scaling_initial, (YREGION - 1) - yc, (XREGION - 1) - xc, multi_peak_flag, DIFFERENCEMETHOD, COMPILED);
       
end



