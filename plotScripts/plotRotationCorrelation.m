imDir = '~/Desktop/combined_transformations';
imBase = 'combined_transformations__';
ext = '.tiff';
plotDir = '~/Desktop/fmcplots3';

nDigits = 6;
startImage = 1;
endImage = 300;
fSize = 40;


numberFormat = ['%0' num2str(nDigits) '.0f'];
imageNumbers = startImage : endImage;
nImages = length(imageNumbers);

nWedges = 128;
nRings = 128;

spatialWindow = gaussianWindowFilter([128, 128], [0.5, 0.5], 'fraction');
fmiWindow = gaussianWindowFilter([nWedges, nRings], [1, 0.8], 'fraction');
imageSpectralFilter = spectralEnergyFilter(128, 128, 2.8);
fmiSpectralFilter = spectralEnergyFilter(nWedges, nRings, 3.3);

for k = 1 : nImages
   imagePaths(k, :) = fullfile(imDir, [imBase num2str(k, numberFormat) ext]);
    
end

% First image
im1 = double(imread(imagePaths(k, :)));
im1win = spatialWindow .* im1;
fmi1 = im2fmi(im1win, nWedges, nRings);
fmiWin1 = fmiWindow .* fmi1;

figure('visible', 'off');
for k = 2 : nImages
    disp(['Processing image ' num2str(k) ' of ' num2str(nImages) ]);
    im2 = double(imread(imagePaths(k, :)));
    im2win = spatialWindow .* im2;
    
    [fmi2, mag2] = im2fmi(im2win, nWedges, nRings);
    magFilt = mag2 .* fmiSpectralFilter;
    fmiWin2 = fmiWindow .* fmi2; 
    
    % Do the correlation
    [~, ~, fmiTy, fmiTx, fmiCorr] = FMIRPC(fmiWin1, fmiWin2, fmiSpectralFilter);
    
    subplot(2, 2, 1);
    imagesc(im2);
%     title('Original Image', 'FontSize', fSize);
    axis image; axis off; colormap gray;  freezeColors;
    
    subplot(2, 2, 2);
    imagesc(mag2);
%     title('Spectral Magnitude', 'FontSize', fSize);
    axis image; axis off; colormap jet; caxis([0 2e6]); freezeColors;
    
    subplot(2, 2, 3);
    imagesc(fmi2);
%     title({'Fourier-Mellin transform (FMT)', '(Log-polar transform of magnitude)'}, 'FontSize', fSize);
    axis image; axis off; colormap jet; caxis([0 2e6]); freezeColors;
    
    subplot(2, 2, 4);
    mesh(fmiCorr ./ max(fmiCorr(:)), 'EdgeColor', 'black'); 
%     title({'Cross Correlation', 'of the pair of FMTs' }, 'FontSize', fSize);
    axis off; freezeColors; xlim([1 nWedges]); ylim([1 nRings]);
    
    set(gcf, 'color', 'white')
    set(gcf, 'OuterPosition', [-1919         241        1920        1200]);

    
    fileName = ['plot_' num2str(k-1, numberFormat) '.png'];
    filePath = fullfile(plotDir, fileName);
    export_fig(filePath);

%     pause(0.01)
    hold off
    
end

hold off