% close all;

% Repository = '~/Desktop/FMCtest_2013-01-24';
testDir = 'FMCtest_2013-01-24';


% Repository = '~/Desktop/FMCtest_2013-01-30_singleSet';
Repository = fullfile(pwd, '..', '..', testDir);
%

load(fullfile(Repository, 'results', '128x128', 'mc', 'statistics', 'errorStatistics_mc_128x128.mat') );

rpcData = load(fullfile(Repository, 'results', '128x128', 'mc', 'statistics', 'errorStatistics_mc_128x128_rpcOnly.mat') );

NiFiFsRpc = rpcData.NiFiFs;
trueRotationRpc = rpcData.trueRotation;

fsize = 28 ;
markerSize = 20;
lineWidth = 2;

% Peak ratio threshold for detectability
peakRatioThreshold = 1.2;
rotationErrorThreshold = 0.02;
scalingErrorThreshold = 0.02;
translationErrorThreshold = 0.15;

plotMarker = 'ok';
% plotMarkerRpc = 'ok';


% Number of bins for the histograms
NB = 20;

% Edges of the histogram of NiFi
histEdges = linspace(0, 30, NB);

rotationHistEdges = linspace(0, 2*pi, NB);

% Build the histogram of NiFi. bin is the bin into which each element of NiFi falls. 
[ n, bin ] = histc(NiFiFs, histEdges);
[nRpc, binRpc] = histc(NiFiFsRpc, histEdges);

[nRot, binRot] = histc(trueRotation, rotationHistEdges);
[nRotRpc, binRotRpc] = histc(trueRotationRpc, rotationHistEdges);

nBins = length(n); % Number of bins
nBinsRpc = length(nRpc);

% Binary vector specifying whether or not the measurement exceeded the
% correlation plane threshold
Detected = peakRatio >= peakRatioThreshold;
correctRotation = rotationAbsError <= rotationErrorThreshold;
correctScaling = scalingAbsError <= scalingErrorThreshold;
correctTranslationX = translationAbsErrorX <= translationErrorThreshold;
correctTranslationY = translationAbsErrorY <= translationErrorThreshold;

translationAbsErrorXrpc = rpcData.translationAbsErrorX;
translationAbsErrorYrpc = rpcData.translationAbsErrorY;

correctTranslationXrpc = translationAbsErrorXrpc <= translationErrorThreshold;
correctTranslationYrpc = translationAbsErrorYrpc <= translationErrorThreshold;

% Initialize vectors to hold statistics
nDetected = zeros(nBins, 1);
nCorrectRotation = zeros(nBins, 1);
nCorrectScaling = zeros(nBins, 1);
nCorrectTranslationX = zeros(nBins, 1);
nCorrectTranslationY = zeros(nBins, 1);
nCorrectTranslationXrotation = zeros(nBins, 1);
nCorrectTranslationYrotation = zeros(nBins, 1);

nCorrectTranslationXrpc = zeros(nBinsRpc, 1);
nCorrectTranslationYrpc = zeros(nBinsRpc, 1);
nCorrectTranslationXrotationRpc = zeros(nBinsRpc, 1);
nCorrectTranslationYrotationRpc = zeros(nBinsRpc, 1);

pDetect = zeros(nBins, 1);
pCorrectRotation = zeros(nBins, 1);
pCorrectScaling = zeros(nBins, 1);
pCorrectTranslationX = zeros(nBins, 1);
pCorrectTranslationY = zeros(nBins, 1);
pCorrectTranslationXrotation = zeros(nBins, 1);
pCorrectTranslationYrotation = zeros(nBins, 1);


pCorrectTranslationXrpc = zeros(nBinsRpc, 1);
pCorrectTranslationYrpc = zeros(nBinsRpc, 1);
pCorrectTranslationXrotationRpc = zeros(nBinsRpc, 1);
pCorrectTranslationYrotationRpc = zeros(nBinsRpc, 1);

meanErrorTranslationX = zeros(nBins, 1);
meanErrorTranslationY = zeros(nBins, 1);

meanErrorTranslationXrpc = zeros(nBinsRpc, 1);
meanErrorTranslationYrpc = zeros(nBinsRpc, 1);

meanErrorRotation = zeros(nBins, 1);
meanErrorScaling = zeros(nBins, 1);


for k = 1:nBins
    nDetected(k) = sum(Detected(bin == k)); % Fraction of measurements whose peak ratios exceeded the specified threshold
    nCorrectRotation(k) = sum(correctRotation(bin == k));
    nCorrectScaling(k) = sum(correctScaling(bin == k));
    nCorrectTranslationX(k) = sum(correctTranslationX( bin == k));
    nCorrectTranslationY(k) = sum(correctTranslationY( bin == k));
    

    
    pDetect(k) = 100 * nDetected(k) ./ n(k); % Fraction of vectors that were detected
    pCorrectRotation(k) = 100 * nCorrectRotation(k) ./ n(k);
    pCorrectScaling(k) = 100 * nCorrectScaling(k) ./ n(k);
    pCorrectTranslationX(k) = 100 * nCorrectTranslationX(k) ./ n(k);
    pCorrectTranslationY(k) = 100 * nCorrectTranslationY(k) ./ n(k);
    

    
    meanErrorTranslationX(k) = mean(translationAbsErrorX(bin == k));
    meanErrorTranslationY(k) = mean(translationAbsErrorY(bin == k));
%     

    
    meanErrorRotation(k) = mean(rotationAbsError(bin == k));
    meanErrorScaling(k) = mean(scalingAbsError(bin == k));
    
end


for k = 1 : nBinsRpc
    nCorrectTranslationXrpc(k) = sum(correctTranslationXrpc(binRpc == k));
    nCorrectTranslationYrpc(k) = sum(correctTranslationYrpc(binRpc == k));
    
    pCorrectTranslationXrpc(k) = 100 * nCorrectTranslationXrpc(k) ./ n(k);
    pCorrectTranslationYrpc(k) = 100 * nCorrectTranslationYrpc(k) ./ n(k);

    meanErrorTranslationXrpc(k) = mean(translationAbsErrorXrpc(binRpc == k));
    meanErrorTranslationYrpc(k) = mean(translationAbsErrorYrpc(binRpc == k));

end


% figure(4);
% plot(histEdges, pDetect, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth);
% xlabel('NiFi', 'fontsize', 2 * fsize);
% ylabel('Probability of valid detection',  'fontsize', 2 * fsize);
% title('Probability of valid peak detection', 'fontsize', 2 * fsize);
% set(gcf, 'Color', [1 1 1]);
% ylim([0 105]);
% set(gca, 'fontsize', 2* fsize)
% hold on;

% hold on;
% 
% figure(1);
% plot(histEdges, pCorrectTranslationX, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth);
% xlabel('NiFi', 'fontsize', fsize);
% ylabel('Probability',  'fontsize', fsize);
% title({'Probability of correct translation estimate'; ['Error < ' num2str(translationErrorThreshold, '%0.3f') ' pixels']},  'fontsize', fsize);
% set(gcf, 'Color', [1 1 1]);
% ylim([0 105]);
% set(gca, 'fontsize', fsize);

% hold on;

% 
figure(2);
subplot(2, 2, 1)
semilogy(histEdges, meanErrorTranslationX, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'g');
hold on;
semilogy(histEdges, meanErrorTranslationXrpc, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'r');
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Average Error (Pixels)',  'fontsize', fsize);
title('Average Error (Horizontal Displacement)',  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0.01 100]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;
legend('FMC', 'RPC', 'FontSize', fsize)


subplot(2, 2, 2)
semilogy(histEdges, meanErrorTranslationY, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth,  'MarkerFaceColor', 'g');
hold on;
semilogy(histEdges, meanErrorTranslationYrpc, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'r');
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Average Error (Pixels)',  'fontsize', fsize);
title('Average Error (Vertical Displacement)',  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0.01 100]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;
legend('FMC', 'RPC', 'FontSize', fsize)

subplot(2, 2, 3)
semilogy(histEdges, meanErrorRotation, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth,  'MarkerFaceColor', 'g');
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Average Error (Radians)',  'fontsize', fsize)
title('Average Error (Rotation)',  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0.001 1]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;

subplot(2, 2, 4)
semilogy(histEdges, meanErrorScaling, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth,  'MarkerFaceColor', 'g');
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Average Error (Unitless)',  'fontsize', fsize)
title('Average Error (Scaling)',  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0.001 1]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;


figure(3)
subplot(2, 2, 1)
plot(histEdges, pCorrectTranslationX, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'g');
hold on;
plot(histEdges, pCorrectTranslationXrpc, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'r');
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Probability',  'fontsize', fsize);
title(['Probability of detecting horizontal translation (error < ' num2str(translationErrorThreshold, '%3.2f')  ' pix )'],  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0 105]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;
legend('FMC', 'RPC', 'FontSize', fsize);
hold on;
plot([0 20], [100 100], '--k', 'LineWidth', lineWidth);
hold off;


subplot(2, 2, 2)
plot(histEdges, pCorrectTranslationY, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'g');
hold on;
plot(histEdges, pCorrectTranslationYrpc, plotMarker, 'markerSize', markerSize, 'lineWidth', lineWidth , 'MarkerFaceColor', 'r');
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Probability',  'fontsize', fsize);
title(['Probability of detecting vertical translation (error < ' num2str(translationErrorThreshold, '%3.2f')  ' pix )'],  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0 105]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;
legend('FMC', 'RPC',  'FontSize', fsize);
hold on;
plot([0 20], [100 100], '--k', 'LineWidth', lineWidth);
hold off;


subplot(2, 2, 3)
plot(histEdges, pCorrectRotation, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth,  'MarkerFaceColor', 'g');
hold on;
plot([0 20], [100 100], '--k', 'LineWidth', lineWidth);
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Probability',  'fontsize', fsize);
title(['Probability of detecting rotation (error < ' num2str(rotationErrorThreshold, '%3.3f')  ' radians )'],  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0 105]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;

subplot(2, 2, 4)
plot(histEdges, pCorrectScaling, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth,  'MarkerFaceColor', 'g');
hold on;
plot([0 20], [100 100], '--k', 'LineWidth', lineWidth);
hold off
xlabel('N_iF_iF_s', 'fontsize', fsize);
ylabel('Probability',  'fontsize', fsize);
title(['Probability of detecting scaling ( error < ' num2str(scalingErrorThreshold, '%3.3f') ' )'],  'fontsize', fsize);
set(gcf, 'Color', [1 1 1]);
ylim([0 105]);
xlim([0 20]);
set(gca, 'fontsize', fsize);
grid on;


% 
% figure(3);
% plot(histEdges, pDetect, plotMarker,  'markerSize', markerSize, 'lineWidth', lineWidth);
% xlabel('NiFi', 'fontsize', 2 * fsize);
% ylabel('Probability of valid detection',  'fontsize', 2 * fsize);
% title(['Probability of generating a correlation plane whose peak ratio is greater than ' num2str(peakRatioThreshold, '%3.1f')], 'fontsize', 2 * fsize);
% set(gcf, 'Color', [1 1 1]);
% ylim([0 105]);
% set(gca, 'fontsize',2 *  fsize)



