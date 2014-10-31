
function rotationErrorHistogram(JOBFILE)

% Plot font size
fSize = 16;
fmiRpcDiameter = JOBFILE.Parameters.Processing.FmiRPCDiameter;
setType = JOBFILE.SetType;
caseName = JOBFILE.CaseName;
correlationType = JOBFILE.CorrelationType;
imageType = JOBFILE.ImageType;
regionHeight = JOBFILE.Parameters.RegionHeight;
regionWidth =  JOBFILE.Parameters.RegionWidth;
startSet = JOBFILE.Parameters.Sets.Start;
endSet = JOBFILE.Parameters.Sets.End;

setDigits = 5;
setFormat = ['%0' num2str(setDigits) '.0f'];

Repository = determineLocalRepositoryPath;

DataDir = fullfile(Repository,'analysis', 'data', imageType, setType, caseName, [num2str(regionHeight) 'x' num2str(regionWidth)], 'stat');

FileName = ['errorStatistics_' setType '_' correlationType...
    '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_' ...
    num2str(startSet, setFormat) '-' num2str(endSet, setFormat)  '.mat']; 

Data = load(fullfile(DataDir, FileName));

trueRotationDegrees = rad2deg(Data.trueRotation); % Work in degrees
trueRotationPixels = abs(Data.trueRotationPixels);
rotationAbsErrorDegrees = rad2deg(Data.rotationAbsError); % Work in degrees
rotationPixelsAbsError = Data.rotationPixelsAbsError;
rotationRelError = rotationAbsErrorDegrees ./ trueRotationDegrees;
rotationRelErrorPixels = rotationPixelsAbsError ./ trueRotationPixels;
NW = Data.numWedges;
% NWmax = max(NWoriginal);
% NWmin = min(NWoriginal);
% nImages = length(NWoriginal);
% 
% NW = NWmin + (NWmax - NWmin) .* rand(nImages, 1);


logRot = log10(trueRotationDegrees);
logRotPix = log10(trueRotationPixels);
% logErr = log10(rad2deg(rotationAbsError));

nRotBins = 50;
nRotPixBins = 50;
nWedgeBins = 50;

rotBinEdges = (linspace(min(logRot), max(logRot), nRotBins))';
rotPixBinEdges = (linspace(min(logRotPix), max(logRotPix), nRotPixBins))';
wedgeBinEdges = (linspace(min(NW), max(NW), nWedgeBins))';

[~, rotBins] = histc(logRot, rotBinEdges);
[~, rotPixBins] = histc(logRotPix, rotPixBinEdges);
[~, wedgeBins] = histc(NW, wedgeBinEdges);

absErrorRotationPixels = rotationAbsErrorDegrees .* NW / 360;

meanAbsErrorRotationDegrees = zeros(nRotBins, nWedgeBins);
meanAbsErrorRotationPixels = zeros(nRotPixBins, nWedgeBins);

meanRelErrorRotation = zeros(nRotBins, nWedgeBins);
meanRelErrorRotationPixels = zeros(nRotPixBins, nWedgeBins);

stdAbsErrorRotationDegrees = zeros(nRotBins, nWedgeBins);
stdRelErrorRotation = zeros(nRotBins, nWedgeBins);

parfor x = 1 : nRotBins
    fprintf(1, [num2str(x) ' of ' num2str(nRotBins) '\n']);
    for y = 1 : nWedgeBins
        meanAbsErrorRotationDegrees(y, x) = mean(rotationAbsErrorDegrees( rotBins == x & wedgeBins == y));
        meanRelErrorRotation(y, x) = mean(rotationRelError( rotBins == x & wedgeBins == y));

        meanAbsErrorRotationPixels(y, x) = mean(absErrorRotationPixels( rotBins == x & wedgeBins == y));
        meanRelErrorRotationPixels(y, x) = mean( rotationRelErrorPixels (rotPixBins == x & wedgeBins == y) );
        
        stdAbsErrorRotationDegrees(y, x) = std(rotationAbsErrorDegrees( rotBins == x & wedgeBins == y));
        stdRelErrorRotation(y, x) = std(rotationRelError( rotBins == x & wedgeBins == y));
        stdAbsErrorRotationPixels(y, x) = std( rotationPixelsAbsError (rotPixBins == x & wedgeBins == y) );
        stdRelErrorRotationPixels(y, x) = std( rotationRelErrorPixels (rotPixBins == x & wedgeBins == y) );
%         
%         fprintf(1, [num2str(k) ', ' num2str(m) '\n']);
        
    end
    
end

logAbsErrDegrees = log10(meanAbsErrorRotationDegrees);
logRelErr = log10(meanRelErrorRotation);
logAbsErrPixels = log10(meanAbsErrorRotationPixels);
logRelErrPixels = log10(abs(meanRelErrorRotationPixels));

logStdAbsErr = log10(stdAbsErrorRotationDegrees);
logStdRelErr = log10(stdRelErrorRotation);
logStdAbsErrPix = log10(stdAbsErrorRotationPixels);
logStdRelErrPix = log10(stdRelErrorRotationPixels);

relWedgeEdges = wedgeBinEdges / Data.imageHeight(1);

% imagesc(rotBinEdges, relWedgeEdges / pi, logAbsErr); 
% axis square;
% set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
% caxis([-3 0]);
% t = colorbar('peer', gca); 
% set(get(t, 'ylabel'), 'String', 'Log_{10} of Absolute Error in Rotation (deg)', 'FontSize', fSize);
% xlabel('Log_{10} of rotation angle (deg)', 'FontSize', fSize);
% ylabel('Number of Wedges / Image Diameter (pix) / \pi', 'FontSize', fSize);
% title('Absolute rotation error (10^6 subregions)', 'FontSize', fSize);
% hold on; 
% plot([min(logRot) max(logRot)], [1 1], '--k'); 
% plot([min(logRot) max(logRot)], [1/pi 1/pi], '--k'); 
% set(gca, 'FontSize', fSize);
% hold off;

% Relative error of rotation
figure(1);
%%% DON'T FORGET TO DIVIDE BY PI AFTER DEBUGGING
% imagesc(rotBinEdges(1:end-1), relWedgeEdges(1:end-1) * regionHeight, logRelErr(1:end-1, 1:end-1)); 
contourf(rotBinEdges(1:end-1), relWedgeEdges(1:end-1) * regionHeight, logRelErr(1:end-1, 1:end-1)); 
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([-3 1]);
xlim([-2 2]);
% ylim([0 1.5])
t = colorbar('peer', gca); 
set(get(t, 'ylabel'), 'String', 'Log_{10} of Relative Error in Rotation', 'FontSize', fSize);
xlabel('Log_{10} of rotation angle (deg)', 'FontSize', fSize);
ylabel('Number of Wedges', 'FontSize', fSize);
title({['Relative rotation error, D_{RPC} = ' num2str(fmiRpcDiameter, '%0.1f')], ...
    [num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation regions']}, 'FontSize', fSize);
hold on; 
plot([min(logRot) max(logRot)], [1 1], '--k'); 
plot([min(logRot) max(logRot)], [1/pi 1/pi], '--k'); 
set(gca, 'FontSize', fSize);
hold off;

figure(2);
% imagesc(rotBinEdges(1:end-1), relWedgeEdges(1:end-1) * regionHeight, logAbsErrDegrees(1:end-1, 1:end-1)); 
contourf(rotBinEdges(1:end-1), relWedgeEdges(1:end-1) * regionHeight, logAbsErrDegrees(1:end-1, 1:end-1)); 
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([-3 1]);
xlim([-2 2]);
% ylim([0 1.5])
t = colorbar('peer', gca); 
set(get(t, 'ylabel'), 'String', 'Log_{10} of Absolute Error in Rotation (degrees)', 'FontSize', fSize);
xlabel('Log_{10} of rotation angle (deg)', 'FontSize', fSize);
ylabel('Number of Wedges', 'FontSize', fSize);
title({['Rotation error (Deg), D_{RPC} = ' num2str(fmiRpcDiameter, '%0.1f') ], ...
    [num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation regions']}, 'FontSize', fSize);
hold on; 
plot([min(logRot) max(logRot)], [1 1], '--k'); 
plot([min(logRot) max(logRot)], [1/pi 1/pi], '--k'); 
set(gca, 'FontSize', fSize);
hold off;


figure(3);
contourf(rotBinEdges(1:end-1), relWedgeEdges(1:end-1) * regionHeight, logAbsErrPixels(1:end-1, 1:end-1)); 
% contourf(relWedgeEdges(1:end-1), rotBinEdges(1:end-1), logRelErr(1:end-1, 1:end-1)); 
axis square;
set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
caxis([-3 0]);
xlim([-2 2]);
% ylim([0 1.5])
t = colorbar('peer', gca); 
set(get(t, 'ylabel'), 'String', 'Absolute Error in Rotation (pixels)', 'FontSize', fSize);
xlabel('Log_{10} of rotation angle (deg)', 'FontSize', fSize);
ylabel('Number of Wedges', 'FontSize', fSize);
title({['Rotation error (Pix), D_{RPC} = ' num2str(fmiRpcDiameter, '%0.1f')] , ...
    [num2str(regionHeight) ' x ' num2str(regionWidth) ' pixel interrogation regions']}, 'FontSize', fSize);
hold on; 
plot([min(logRot) max(logRot)], [1 1], '--k'); 
plot([min(logRot) max(logRot)], [1/pi 1/pi], '--k'); 
set(gca, 'FontSize', fSize);
hold off;

% figure;
% imagesc(rotPixBinEdges, relWedgeEdges / pi, logRelErrPix); 
% axis square;
% set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
% caxis([-3 0]);
% t = colorbar('peer', gca); 
% set(get(t, 'ylabel'), 'String', 'Log_{10} of Average Rel Error in Rotation (pix)', 'FontSize', fSize);
% xlabel('Log_{10} of rotation (pix)', 'FontSize', fSize);
% ylabel('Number of Wedges / Image Diameter (pix) / \pi', 'FontSize', fSize);
% title('Average rotation error (pixels) (10^6 subregions)', 'FontSize', fSize);
% hold on; 
% plot([min(logRotPix) max(logRotPix)], [1 1], '--w'); 
% plot([min(logRotPix) max(logRotPix)], [1/pi 1/pi], '--w'); 
% set(gca, 'FontSize', fSize);
% hold off;

% figure;
% imagesc(rotPixBinEdges, relWedgeEdges / pi, logStdAbsErrPix); 
% axis square;
% set(gca, 'ydir', 'normal');set(gcf, 'color', 'white')
% caxis([-3 0]);
% t = colorbar('peer', gca); 
% set(get(t, 'ylabel'), 'String', 'Log_{10} Std of Rel Error in Rotation (pix)', 'FontSize', fSize);
% xlabel('Log_{10} of rotation (pix)', 'FontSize', fSize);
% ylabel('Number of Wedges / Image Diameter (pix) / \pi', 'FontSize', fSize);
% title('Std rotation error (pixels) (10^6 subregions)', 'FontSize', fSize);
% hold on; 
% plot([min(logRotPix) max(logRotPix)], [1 1], '--w'); 
% plot([min(logRotPix) max(logRotPix)], [1/pi 1/pi], '--w'); 
% set(gca, 'FontSize', fSize);
% hold off;

% saveName = ['MCHist_' num2str(startSet, setFormat) '-' num2str(endSet, setFormat) '.mat'];
% 
% save(fullfile(DataDir, saveName));
    
end



