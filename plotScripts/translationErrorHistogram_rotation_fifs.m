
function translationErrorHistogram_rotation_fifs(JOBLIST)

% Specify number of bins
nRotBins = 50;
nRotPixBins = 50;
nFiFsBins = 50;

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
trueRotationDegrees = abs(rad2deg(Data.trueRotation)); % Work in degrees
trueRotationPixels = abs(Data.trueRotationPixels);

% True translations (components)
trueTranslationX = Data.trueTranslationX;
trueTranslationY = Data.trueTranslationY;

% Mutual information paramters.
NiFiFs = Data.NiFiFs

% Translation errors (components)
translationAbsErrorX = Data.translationAbsErrorX;
translationAbsErrorY = Data.translationAbsErrorY;

% Magnitudes of translation and translation error
trueTranslationMagnitude = sqrt(trueTranslationX.^2 + trueTranslationY.^2);
translationAbsErrorMagnitude = sqrt(translationAbsErrorX.^2 + translationAbsErrorY.^2);

% Relative error of translation estimate
translationRelErrorMagnitude = translationAbsErrorMagnitude ./ trueTranslationMagnitude;

% Create histogram bin edges
rotBinEdges = (linspace(0, ceil(max(trueRotationDegrees)), nRotBins))';
% translationBinEdges = (linspace(min(trueTranslationMagnitude), max(trueTranslationMagnitude), nTranslationBins))';
translationBinEdges = (linspace(0, floor(max(trueTranslationMagnitude)), nFiFsBins))';

% Create error histograms
[~, rotBins] = histc(trueRotationDegrees, rotBinEdges);
[~, translationBins] = histc(trueTranslationMagnitude, translationBinEdges);

% Average errors
meanAbsErrorTranslation = zeros(nRotBins, nFiFsBins);
meanRelErrorTranslation = zeros(nRotBins, nFiFsBins);

% Standard deviations of errors
stdAbsErrorTranslation = zeros(nRotBins, nFiFsBins);
stdRelErrorTranslation = zeros(nRotBins, nFiFsBins);

% Calculate the statistics for each bin
parfor x = 1 : nRotBins
    fprintf(1, [num2str(x) ' of ' num2str(nRotBins) '\n']);
    for y = 1 : nFiFsBins
        meanAbsErrorTranslation(y, x) = mean(translationAbsErrorMagnitude( rotBins == x & translationBins == y ));
        meanRelErrorTranslation(y, x) = mean(translationRelErrorMagnitude( rotBins == x & translationBins == y));

        stdAbsErrorTranslation(y, x) = std(translationAbsErrorMagnitude( rotBins == x & translationBins == y));
        stdRelErrorTranslation(y, x) = std(translationRelErrorMagnitude( rotBins == x & translationBins == y));

    end
    
end

% Take the log of the translation errors
logAbsErrTranslation = log10(meanAbsErrorTranslation);
logRelErrTranslation = log10(meanRelErrorTranslation);

% Absolute error of translation
figure(1);
contourf( rotBinEdges(1:end-1), translationBinEdges(1:end-1), logAbsErrTranslation(1:end-1, 1:end-1)); 
% surf( rotBinEdges(1:end), translationBinEdges(1:end), logAbsErrTranslation(1:end, 1:end));
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([-2 1]);
% Make a colorbar and label it.
h = colorbar;
set(get(h, 'ylabel'), 'String', 'Translation errror (pixels)', 'FontSize', fSize);

% These are the positions at which to place tick marks in the colorbar
ticksLin = log10([0.01 0.1, 0.5, 1, 5, 10]);
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
ylabel('True translation matnigude (pixels)', 'FontSize', fSize);
title({'Absolute error of translation estimate' ; ...
    [upper(JOBFILE.CorrelationType) ' algorithm, ' num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation region']},...
    'FontSize', fSize);
set(gca, 'FontSize', fSize);
xlim([0 85]);
ylim([0 20]);

% Save the figure
print(1, '-depsc', ['~/Desktop/Histograms/translationError_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '.eps']);

% savefig(1, ['~/Desktop/translationError_' num2str(regionHeight) 'x' num2str(regionWidth) '_' JOBFILE.CorrelationType '.fig']);


% figure(1)
% semilogy(trueTranslationMagnitude, translationAbsErrorMagnitude, '.k');
% xlabel('True translation magnitude (pix)', 'FontSize', 16)
% ylabel('Translation absolute error (pix)', 'FontSize', 16);
% title('Absolute error of translation estimate, RPC algorithm, 128x128 pix', 'FontSize', 16)
% set(gcf, 'color', 'white')
end

end


