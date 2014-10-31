
function scalingErrorHistogram_rotation_scaling(JOBLIST)

% Specify number of bins
nScalingBins = 50;
% nRotPixBins = 50;
nRotBins = 50;

% Plot font size
fSize = 16;

% Number of jobs to plot.
nJobs = length(JOBLIST);

for n = 1 : nJobs
    
JOBFILE = JOBLIST(n);

regionHeight = JOBFILE.Parameters.RegionHeight;
regionWidth = JOBFILE.Parameters.RegionWidth;


% Extract parameters from Jobfile
fmiRpcDiameter = JOBFILE.Parameters.Processing.FmiRPCDiameter;
setType = JOBFILE.SetType;
caseName = JOBFILE.CaseName;
correlationType = JOBFILE.CorrelationType;
imageType = JOBFILE.ImageType;
regionHeight = JOBFILE.Parameters.RegionHeight;
regionWidth =  JOBFILE.Parameters.RegionWidth;
startSet = JOBFILE.Parameters.Sets.Start;
endSet = JOBFILE.Parameters.Sets.End;

% Specify the set number format.
setDigits = 5;
setFormat = ['%0' num2str(setDigits) '.0f'];

% Determine the local path to the project repository
Repository = determineLocalRepositoryPath;

% Determine the path to the directory containing the data file to be loaded
DataDir = fullfile(Repository,'analysis', 'data', imageType, setType, caseName, [num2str(regionHeight) 'x' num2str(regionWidth)], 'stat');

% Determine the name of the data file to load
FileName = ['errorStatistics_' setType '_' correlationType...
    '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_' ...
    num2str(startSet, setFormat) '-' num2str(endSet, setFormat)  '.mat']; 

% Load the data file
Data = load(fullfile(DataDir, FileName));

% True and estimated scaling
trueScaling = Data.trueScaling;
estimatedScaling = Data.estimatedScaling;

% True and est. rotations
trueRotationRad = Data.trueRotation;
trueRotationDegrees = rad2deg(trueRotationRad);
estimatedRotationRad = Data.estimatedRotation;
rotationAbsErrorMagnitude = rad2deg(angleAbsDiff(trueRotationRad, estimatedRotationRad));
scalingAbsError = abs(trueScaling - estimatedScaling);

% Create histogram bin edges
scalingBinEdges = (linspace(min(trueScaling), max(trueScaling), nScalingBins))';
% scalingBinEdges = (linspace(0.5, 2.0, nScalingBins))';
% translationBinEdges = (linspace(min(trueTranslationMagnitude), max(trueTranslationMagnitude), nTranslationBins))';
rotBinEdges = (linspace(0, ceil(max(trueRotationDegrees)), nRotBins))';

% Create error histograms
[~, scalingBins] = histc(trueScaling, scalingBinEdges);
[~, rotBins] = histc(abs(trueRotationDegrees), rotBinEdges);

% Initialize rotation error
meanAbsErrorScaling = zeros(nRotBins, nScalingBins);
stdAbsErrorScaling = zeros(nRotBins, nScalingBins);

% Calculate the statistics for each bin
parfor x = 1 : nRotBins
    fprintf(1, [num2str(x) ' of ' num2str(nScalingBins) '\n']);
    for y = 1 : nScalingBins
        meanAbsErrorScaling(y, x) = mean(scalingAbsError( scalingBins == y & rotBins == x ));
        stdAbsErrorScaling(y, x) = std(scalingAbsError( scalingBins == y & rotBins == x ));
    end    
end

logAbsErrorScaling = log10(meanAbsErrorScaling);

% Contour levels
% contourLevels = linspace(log10(0.0005), 0.1, log10(0.1));
contourLevels = log10(0.0005): 0.1 : log10(0.1);

% Absolute error of translation
f = figure(1);

% Make contour plot
contourf( rotBinEdges(1:end), scalingBinEdges(1:end), logAbsErrorScaling(1:end, 1:end), contourLevels); 
% contourf( rotBinEdges(1:end), scalingBinEdges(1:end), logAbsErrorScaling(1:end, 1:end)); 

% Format the plot
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
% caxis([-2 1]);
caxis([log10(0.0005), log10(0.1)]);
% Make a colorbar and label it.
h = colorbar;
set(get(h, 'ylabel'), 'String', 'Scaling error', 'FontSize', fSize);

% These are the positions at which to place tick marks in the colorbar
% ticksLin = log10([0.01 0.1, 0.5, 1, 5, 10]);
ticksLin = log10([0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1]);
set(h, 'YTick', ticksLin);

% Populate the colorbar tick label strings. 
nTicks = length(ticksLin);
for k = 1 : nTicks;
    tickLabels{k} = num2str(10^ticksLin(k));
end

% Set the color bar labels
set(h, 'ytickLabel', tickLabels);

% Format the other axes.
xlabel('True rotation magnitude (degrees)', 'FontSize', fSize);
ylabel('True scaling factor', 'FontSize', fSize);
title({'Average absolute errors of scaling estimate ' ; ...
    [upper(JOBFILE.CorrelationType) ' algorithm, ' num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation region']},...
    'FontSize', fSize);
set(gca, 'FontSize', fSize);
xlim([0 180]);
ylim([0.5 2]);

% Save the figure
print(1, '-depsc', ['~/Desktop/Histograms/eps/scalingError_rotation_scaling_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.eps']);
saveas(f, ['~/Desktop/Histograms/fig/scalingError_rotation_scaling_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.fig']);

% Save the histogram data
save(['~/Desktop/Histograms/data/scalingError_rotation_scaling_hist_' caseName...
    '_' num2str(regionHeight) 'x' num2str(regionWidth)...
    '_' JOBFILE.CorrelationType '_linear.mat'], 'JOBFILE', 'logAbsErrorScaling', 'scalingBinEdges', 'rotBinEdges', 'caseName', 'contourLevels', 'regionHeight', 'regionWidth', 'fSize');

% figure(1)
% semilogy(trueTranslationMagnitude, translationAbsErrorMagnitude, '.k');
% xlabel('True translation magnitude (pix)', 'FontSize', 16)
% ylabel('Translation absolute error (pix)', 'FontSize', 16);
% title('Absolute error of translation estimate, RPC algorithm, 128x128 pix', 'FontSize', 16)
% set(gcf, 'color', 'white')
end

end


