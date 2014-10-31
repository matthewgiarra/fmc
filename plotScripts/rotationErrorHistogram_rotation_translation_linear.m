
function rotationErrorHistogram_rotation_translation_linear(JOBLIST)

% Specify number of bins
nRotBins = 50;
nRotPixBins = 50;
nTranslationBins = 50;

% Plot font size
fSize = 16;

% Number of jobs to plot.
nJobs = length(JOBLIST);

for n = 1 : nJobs
    
JOBFILE = JOBLIST(n);

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
trueRotationRad = Data.trueRotation;
estimatedRotationRad = Data.estimatedRotation;

trueRotationDegrees = (rad2deg(Data.trueRotation)); % Work in degrees

% True translations (components)
trueTranslationX = Data.trueTranslationX;
trueTranslationY = Data.trueTranslationY;

% Magnitudes of translation and translation error
trueTranslationMagnitude = sqrt(trueTranslationX.^2 + trueTranslationY.^2);
rotationAbsErrorMagnitude = rad2deg(angleAbsDiff(trueRotationRad, estimatedRotationRad));



% Create histogram bin edges
rotBinEdges = (linspace(0, ceil(max(trueRotationDegrees)), nRotBins))';
% rotBinEdges = (linspace(floor(min(trueRotationDegrees)), ceil(max(trueRotationDegrees)), nRotBins))';
translationBinEdges = (linspace(0, ceil(max(trueTranslationMagnitude)), nTranslationBins))';
% translationBinEdges = (linspace(0, 16 * regionHeight / 128, nTranslationBins))';
% Create error histograms
[~, rotBins] = histc(abs(trueRotationDegrees), rotBinEdges);
[~, translationBins] = histc(trueTranslationMagnitude, translationBinEdges);

% Average errors
meanAbsErrorRotation = zeros(nRotBins, nTranslationBins);
stdAbsErrorRotation = zeros(nRotBins, nTranslationBins);

% Calculate the statistics for each bin
parfor x = 1 : nRotBins
    fprintf(1, [num2str(x) ' of ' num2str(nRotBins) '\n']);
    for y = 1 : nTranslationBins
        meanAbsErrorRotation(y, x) = mean(rotationAbsErrorMagnitude( rotBins == x & translationBins == y ));
%         stdAbsErrorRotation(y, x) = std(rotationAbsErrorMagnitude( rotBins == x & translationBins == y ));
    end
    
end

% absErrorRotationUpper = meanAbsErrorRotation + 2.575 * stdAbsErrorRotation;


% Contour levels
% contourLevels = 0 : 0.005 : 0.1;
contourLevels = linspace(0, 0.12, 20);

% Absolute error of translation
f = figure(1);
% contourf( rotBinEdges(1:end), translationBinEdges(1:end), meanAbsErrorRotation(1:end, 1:end), contourLevels); 
contourf( rotBinEdges(1:end), translationBinEdges(1:end), meanAbsErrorRotation(1:end, 1:end), contourLevels ); 
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([0 0.12]);
% Make a colorbar and label it.
h = colorbar;
set(get(h, 'ylabel'), 'String', 'Rotation error (degrees)', 'FontSize', fSize);

% These are the positions at which to place tick marks in the colorbar

% Format the other axes.
xlabel('True rotation magnitude (degrees)', 'FontSize', fSize);
ylabel('True translation matnigude (pixels)', 'FontSize', fSize);
title({'Average absolute error of rotation estimate' ; ...
    [upper(JOBFILE.CorrelationType) ' algorithm, ' num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation region']},...
    'FontSize', fSize);
set(gca, 'FontSize', fSize);
xlim([0 180]);
ylim([0 16 * regionHeight / 128]);


% g = figure(2);
% % Upper limit of confidence interval
% contourf( rotBinEdges(1:end), translationBinEdges(1:end), absErrorRotationUpper(1:end, 1:end), contourLevels); 
% axis square;
% set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
% caxis([0 0.12]);
% % Make a colorbar and label it.
% h = colorbar;
% set(get(h, 'ylabel'), 'String', 'Rotation error (degrees)', 'FontSize', fSize);
% 
% % These are the positions at which to place tick marks in the colorbar
% 
% % Format the other axes.
% xlabel('True rotation magnitude (degrees)', 'FontSize', fSize);
% ylabel('True translation matnigude (pixels)', 'FontSize', fSize);
% title({'Upper limit of 99% confidence interval of rotation error' ; ...
%     [upper(JOBFILE.CorrelationType) ' algorithm, ' num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation region']},...
%     'FontSize', fSize);
% set(gca, 'FontSize', fSize);
% xlim([0 180]);
% ylim([0 16 * regionHeight / 128]);


% Save the figure
print(f, '-depsc', ['~/Desktop/Histograms/eps/rotationError_rotation_translation_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.eps']);
saveas(f, ['~/Desktop/Histograms/fig/rotationError_rotation_translation_hist_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.fig']);

% % Save the figure
% print(g, '-depsc', ['~/Desktop/Histograms/eps/rotationError_rotation_translation_hist_upper_conf_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.eps']);
% saveas(g, ['~/Desktop/Histograms/fig/rotationError_rotation_translation_hist_upper_conf_' caseName '_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '_linear.fig']);

% Save the histogram data
save(['~/Desktop/Histograms/data/rotationError_rotation_translation_hist_' caseName...
    '_' num2str(regionHeight) 'x' num2str(regionWidth)...
    '_' JOBFILE.CorrelationType '_linear.mat'], 'JOBFILE', 'meanAbsErrorRotation', 'translationBinEdges', 'rotBinEdges', 'caseName', 'contourLevels', 'regionHeight', 'regionWidth', 'fSize');

% figure(1)
% semilogy(trueTranslationMagnitude, translationAbsErrorMagnitude, '.k');
% xlabel('True translation magnitude (pix)', 'FontSize', 16)
% ylabel('Translation absolute error (pix)', 'FontSize', 16);
% title('Absolute error of translation estimate, RPC algorithm, 128x128 pix', 'FontSize', 16)
% set(gcf, 'color', 'white')

end

end


