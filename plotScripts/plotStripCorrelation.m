% plotStripCorrelation

NW = 1024;
NR = 64;

fileDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/mc/FMCtest_2013-10-09_smallRotations/128x128/fmc_strip/';
fileName = ['errorAnalysis_mc_fmc_strip_h128_w128_nw' num2str(NW) '_nr' num2str(NR) '_00001'];

saveDir = '~/Desktop/rotationErrorPlots';

% Load the data file
load(fullfile(fileDir, fileName));

height = 128;
width = 128;

maxRadius = ceil(min(height/2, width/2));

% Number of rings
nRings = size(estimatedRotation, 2);

% Number of image pairs
nPairs = size(estimatedRotation, 1);

% Number of wedges
nWedges = JOBFILE.Parameters.Processing.NumberOfWedges.Min;

% Estimated rotation in degrees
estimatedRotationDeg = rad2deg(estimatedRotation);

% Absolute error of rotation estimate in degrees
rotationAbsErrorDeg = rad2deg(rotationAbsError);
errDeg = reshape(rotationAbsErrorDeg, numel(rotationAbsErrorDeg), 1);

% Absolute error of displacement estimate in pix
errPix = nWedges * errDeg / 360;

% True rotation in degrees
trueRotationDeg = rad2deg(repmat(trueRotation, 1, size(estimatedRotationDeg, 2)));

% True peak displacement in pixels
trueDisplacementPixels = nWedges * trueRotationDeg / 360;

% Estimated peak displacement in pixels
peakDisplacementError = abs(peakDisplacement - trueDisplacementPixels);

% Mean error of peak displacement
meanPeakDisplacementError = mean(peakDisplacementError, 1);

% Mean estimated rotation angle (degrees)
meanRotationAbsErrDeg = mean(rotationAbsErrorDeg, 1);

% Make a vector for all the rings
ringNumber = maxRadius / nRings * repmat(1 : nRings, 1000, 1);
r = reshape(ringNumber, numel(ringNumber), 1);

% Make a new figure
figure(1);
% Plot the rotation errors
plot(r, errDeg, 'ok', 'MarkerFaceColor', 'k');
hold on;
errorbar(ringNumber(1, :), meanRotationAbsErrDeg, std(rotationAbsErrorDeg, 1), 'or', 'MarkerFaceColor', 'r');
legend('Raw', 'Mean')
hold off;
xlim([0 max(r) + 1]);
ylim([0 6]);
xlabel('Wave number', 'FontSize', 16);
ylabel('Rotation error (deg)', 'FontSize', 16);
title(['Rotation error in degrees, NW = ' num2str(nWedges) ', NR = ' num2str(nRings)], 'FontSize', 16);
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
axis square;


% Make a new figure
figure(2);
% Plot the rotation errors
plot(r, nWedges / 360 * errDeg, 'ok', 'MarkerFaceColor', 'k');
hold on;
errorbar(ringNumber(1, :), nWedges / 360 * meanRotationAbsErrDeg, nWedges / 360 * std(rotationAbsErrorDeg, 1), 'or', 'MarkerFaceColor', 'r');
legend('Raw', 'Mean')
hold off;
xlim([0 max(r) + 1]);
ylim([0 4]);
xlabel('Wave number', 'FontSize', 16);
ylabel('Rotation error (pix)', 'FontSize', 16);
title(['Rotation error in pixels, NW = ' num2str(nWedges) ', NR = ' num2str(nRings)], 'FontSize', 16);
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
axis square;

% Save the angle error plots
print(1, '-dpng', '-r200', fullfile(saveDir, ['rotationAngleError_nw' num2str(nWedges) '_nr' num2str(nRings) '.png']));
print(1, '-depsc', fullfile(saveDir, ['rotationAngleError_nw' num2str(nWedges) '_nr' num2str(nRings) '.eps']));

% Save the pixel error plots
print(2, '-dpng', '-r200', fullfile(saveDir, ['rotationPixError_nw' num2str(nWedges) '_nr' num2str(nRings) '.png']));
print(2, '-depsc', fullfile(saveDir, ['rotationPixError_nw' num2str(nWedges) '_nr' num2str(nRings) '.eps']));

