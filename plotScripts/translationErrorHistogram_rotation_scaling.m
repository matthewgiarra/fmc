
function translationErrorHistogram_rotation_scaling(JOBLIST)

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

% True rotations
trueRotationMagnitude = abs(rad2deg(Data.trueRotation)); % Work in degrees
trueScaling = Data.trueScaling;

% True translations (components)
trueTranslationX = Data.trueTranslationX;
trueTranslationY = Data.trueTranslationY;

% Translation errors (components)
translationAbsErrorX = Data.translationAbsErrorX;
translationAbsErrorY = Data.translationAbsErrorY;

% Magnitudes of translation and translation error
trueTranslationMagnitude = sqrt(trueTranslationX.^2 + trueTranslationY.^2);
translationAbsErrorMagnitude = sqrt(translationAbsErrorX.^2 + translationAbsErrorY.^2);

% Relative error of translation estimate
translationRelErrorMagnitude = translationAbsErrorMagnitude ./ trueTranslationMagnitude;

% Create histogram bin edges
scalingBinEdges = (linspace(min(trueScaling), max(trueScaling), nScalingBins))';
% scalingBinEdges = (linspace(0.5, 1.5, nScalingBins))';
% translationBinEdges = (linspace(min(trueTranslationMagnitude), max(trueTranslationMagnitude), nTranslationBins))';
rotBinEdges = (linspace(0, ceil(max(trueRotationMagnitude)), nRotBins))';

% Create error histograms
[~, scalingBins] = histc(trueScaling, scalingBinEdges);
[~, rotBins] = histc(trueRotationMagnitude, rotBinEdges);

% Average errors
meanAbsErrorTranslation = zeros(nScalingBins, nRotBins);
meanRelErrorTranslation = zeros(nScalingBins, nRotBins);
stdAbsErrorTranslation = zeros(nScalingBins, nRotBins);

% Standard deviations of errors
stdAbsErrorTranslation = zeros(nScalingBins, nRotBins);
stdRelErrorTranslation = zeros(nScalingBins, nRotBins);

% Take the log of the translation errors
logAbsErrTranslation = log10(meanAbsErrorTranslation);
logRelErrTranslation = log10(meanRelErrorTranslation);

% histDir = fullfile(DataDir, '..', 'hist');
% if ~exist(histDir, 'dir')
%     mkdir(histDir)
% end
% 
% histName = ['translationError_rotation_scaling_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '.mat'];
% histPath = fullfile(histDir, histName);

% Calculate the statistics for each bin
parfor x = 1 : nRotBins
    fprintf(1, [num2str(x) ' of ' num2str(nScalingBins) '\n']);
    for y = 1 : nScalingBins
        meanAbsErrorTranslation(y, x) = mean(translationAbsErrorMagnitude( scalingBins == y & rotBins == x ));
         stdAbsErrorTranslation(y, x) = std(translationAbsErrorMagnitude( scalingBins == y & rotBins == x ));
%         meanRelErrorTranslation(y, x) = mean(translationRelErrorMagnitude( scalingBins == y & rotBins == x));
% 
%         stdAbsErrorTranslation(y, x) = std(translationAbsErrorMagnitude( scalingBins == y & rotBins == x));
%         stdRelErrorTranslation(y, x) = std(translationRelErrorMagnitude( scalingBins == y & rotBins == x));
    end    
end

% Calculate the log of the translation error
logAbsErrTranslation = log10(meanAbsErrorTranslation);

% Contour levels
contourLevels = -3 : 0.1 : 1;

% Absolute error of translation
f = figure(1);

% Make contour plot
contourf( rotBinEdges(1:end-1), scalingBinEdges(1:end-1), logAbsErrTranslation(1:end-1, 1:end-1), contourLevels); 

% Format the plot
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([-2 0]);
% Make a colorbar and label it.
h = colorbar;
set(get(h, 'ylabel'), 'String', 'Translation error (pixels)', 'FontSize', fSize);

% These are the positions at which to place tick marks in the colorbar
% ticksLin = log10([0.01 0.1, 0.5, 1, 5, 10]);
ticksLin = log10([0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1]);
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
title({'Absolute error of translation estimate' ; ...
    [upper(JOBFILE.CorrelationType) ' algorithm, ' num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation region']},...
    'FontSize', fSize);
set(gca, 'FontSize', fSize);
xlim([0 180]);
ylim([0.5 2]);

% Save the figure
print(1, '-depsc', ['~/Desktop/Histograms/eps/translationError_rotation_scaling_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '.eps']);
saveas(f, ['~/Desktop/Histograms/fig/translationError_rotation_scaling_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '.fig']);

% Save the histogram data
save(['~/Desktop/Histograms/data/translationError_rotation_scaling_hist_' caseName...
    '_' num2str(regionHeight) 'x' num2str(regionWidth)...
    '_' JOBFILE.CorrelationType '.mat'], 'JOBFILE', 'logAbsErrTranslation', 'scalingBinEdges', 'rotBinEdges', 'caseName', 'contourLevels', 'regionHeight', 'regionWidth', 'fSize');

% figure(1)
% semilogy(trueTranslationMagnitude, translationAbsErrorMagnitude, '.k');
% xlabel('True translation magnitude (pix)', 'FontSize', 16)
% ylabel('Translation absolute error (pix)', 'FontSize', 16);
% title('Absolute error of translation estimate, RPC algorithm, 128x128 pix', 'FontSize', 16)
% set(gcf, 'color', 'white')
end

end


