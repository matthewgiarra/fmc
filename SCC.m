function [TRANSLATIONY, TRANSLATIONX, PEAKHEIGHTRATIO, SPATIALSCC] = SCC(IMAGE1, IMAGE2)
% [TRANSLATIONY TRANSLATIONX SPATIALRPC] = RPC(IMAGE1, IMAGE2, SPECTRALFILTER) calculates the robust phase correlation between two images
% 
% INPUTS
%   IMAGE1 = First image
%   IMAGE2 = Second image
%   SPECTRALFILTER = RPC Spectral energy filter to apply to the spectral correlation between IMAGE1 and IMAGE2
%
% OUTPUTS
%   TRANSLATIONX = Most probable horizontal translation (in pixels) relating the two images
%   TRANSLATIONY = Most probable vertical translation (in pixels) relating the two images
%   SPATIALRPC = Spatial correlation plane of the RPC between the two images
%
% SEE ALSO
%   spectralEnergyFilter, robustPhaseCorrelation, freq2space, subpixel

% % % % % % % % 
% Begin Function  %
% % % % % % % % 

% Compute the robust phase correlation(RPC) of the windowed images. Report the result in the spectral domain.
spectralScc = crossCorrelation(IMAGE1, IMAGE2);

% Convert the RPC of the two input images from the spectral domain to the spatial domain
SPATIALSCC = freq2space(spectralScc);

% Measure peak height ratio
PEAKHEIGHTRATIO = measurePeakHeightRatio(SPATIALSCC);

% Prana subplixel implmentation
[TRANSLATIONY, TRANSLATIONX] = subpixel(SPATIALSCC, ones(size(SPATIALSCC)), 1, 0); % Subpixel peak location (poorly commented function taken from Prana) 

end


