%#codegen
function PEAK_MATRIX = locate_correlation_peaks(SPATIAL_CORRELATION)
% This function accepts a spatial correlation plane
% and returns the plane with peaks extracted. 
% This tiny function is separate from the main subpixel code
% so that it can be compiled easily.

 %Locate peaks using imregionalmax
PEAK_MATRIX = imregionalmax(SPATIAL_CORRELATION) .* SPATIAL_CORRELATION;

end