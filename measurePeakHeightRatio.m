function [peakRatio, peakHeight1] = measurePeakHeightRatio(SPATIAL_CORRELATION_PLANE, COMPILED)
% This function measures the ratio between the highest and second highest
% peak in a correlation plane. 
% peakHeight is the height of the highest peak
% peakRatio is the ratio of the height of the highest peak to the height of
% the second highest peak.

% Choose whether or not to use compiled codes.
if COMPILED
    %Locate peaks using imregionalmax
    peak_matrix = locate_correlation_peaks_mex(SPATIAL_CORRELATION_PLANE);
else
    peak_matrix = locate_correlation_peaks(SPATIAL_CORRELATION_PLANE);
end

% Find the value and location of the tallest peak.
[peakHeight1, L1] = max(peak_matrix(:)); 

% Remove tallest peak from the second plane, leaving the second-tallest peak 
peak_matrix(L1) = 0; 

% Find the value and location of the second tallest peak.
peakHeight2 = max(peak_matrix(:));

% Calculate the peak height ratio
peakRatio = peakHeight1 / peakHeight2; 

end