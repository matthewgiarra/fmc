function [peakRatio, peakHeight1] = measurePeakHeightRatio(correlationPlane)
% This function measures the ratio between the highest and second highest
% peak in a correlation plane. 
% peakHeight is the height of the highest peak
% peakRatio is the ratio of the height of the highest peak to the height of
% the second highest peak.

% Locate the peaks
Peaks = correlationPlane .* double(imregionalmax(correlationPlane));

% Find the value and location of the tallest peak.
[peakHeight1, L1] = max(Peaks(:)); 

% Remove tallest peak from the second plane, leaving the second-tallest peak 
Peaks(L1) = 0; 

% Find the value and location of the second tallest peak.
peakHeight2 = max(Peaks(:));

% Calculate the peak height ratio
peakRatio = peakHeight1 / peakHeight2; 

end