spatialWindowFraction = JobFile.Parameters.Processing.SpatialWindowFraction;
fftSize = JobFile.Parameters.Processing.FFTSize;
fftHeight = fftSize(1); fftWidth = fftSize(2);
numRings = JobFile.Parameters.Processing.NumberOfRings(1).Min;
numWedges = JobFile.Parameters.Processing.NumberOfWedges(1).Min;
rMin = JobFile.Parameters.Processing.MinimumRadius;

figure(1);
plot(fxTheory, fxTheory - fmiTranslationX(:, 1), 'sk', 'Markerfacecolor', 'red', 'MarkerSize', 3);
% plot(fxTheory, fxTheory - fmiTranslationX(:, 1), 'ok', 'Markerfacecolor', 'green');

xlim([-25 25]);
ylim([-1 1])
axis square
grid on
set(gcf, 'Color', 'White');
set(gca, 'FontSize', 16);
% title2Eval = ['title({''$\Delta x_{theory} = (1 - n_r) \ln (\sigma) / \ln( \frac{R}{r_o})$''', ...
%     ', ''$  ' num2str(imageHeight) ' \times ' num2str(imageWidth) '$ region, ' num2str(100 * spatialWindowFraction(1)) '\% spatial window'' '...
%     ', ''$  ' num2str(fftHeight) ' \times ' num2str(fftWidth) '$ FFT, $n_r = $ ' num2str(numRings) ', $n_w$ = '  num2str(numWedges) ', $r_{min} = $ ' num2str(rMin) ' '' }' , ...
%     ', ''interpreter'', ''latex'', ''fontsize'', 12)'];
% eval(title2Eval);

title2Eval = ['title({''$\Delta x_{theory} = (1 - n_r) \ln (\sigma) / \ln( \frac{R}{r_o})$''', ...
    ', ''Effective win. res. $51 \times 51$ pixels'' '...
    ', ''$  ' num2str(fftHeight) ' \times ' num2str(fftWidth) '$ FFT, $n_r = $ ' num2str(numRings) ', $n_w$ = '  num2str(numWedges) ', $r_{min} = $ ' num2str(rMin) ' '' }' , ...
    ', ''interpreter'', ''latex'', ''fontsize'', 12)'];
eval(title2Eval);


xlabel('$\Delta x_{theory}$ (pixels)',...
    'interpreter', 'latex', 'FontSize', 16)
ylabel('$\Delta x_{theory} - \Delta x_{measured} $    (pixels) ',...
    'interpreter', 'latex', 'fontsize', 16);

plotName = ['fxTheoryError_'...
    num2str(imageHeight) 'x' num2str(imageWidth)...
    '_win' num2str(100 * spatialWindowFraction(1))...
    '_fft' num2str(fftHeight) 'x' num2str(fftWidth)...
    '_nr' num2str(numRings) '_nw' num2str(numWedges)...
    '_rMin' num2str(rMin) '.eps'];

plotDir = '~/Desktop/scalingErrorPlots';

plotPath = fullfile(plotDir, plotName);

% print(1, '-depsc', plotPath);