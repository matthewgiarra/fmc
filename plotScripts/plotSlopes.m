function [a, windowFractions] = findConstant()
JobList = fmcJobList_mc;
JobFile = JobList(end);

windowFractions = (0.1 : 0.1 : 1.0)';
a = zeros(length(windowFractions), 1);

plotDir = '~/Desktop/scalingErrorPlots';

for k = 1:length(windowFractions)
   JobFile.Parameters.Processing.SpatialWindowFraction ...
       = windowFractions(k) * [1 1];
   runFmcJobFileMc(JobFile);
   load('/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/mc/FMCtest_2014-03-31_scaling_only/128x128/fmc/errorAnalysis_mc_fmc_h128_w128_00001.mat');
   
   fmiX = fmiTranslationX(:, 1);
   
   % Determine the constants.
   a(k) = const_finder(trueScaling, fmiX, rMin, rMax, numRings(1));
   
   
   fxTheory = 1 * (1 - numRings(1)) * log(trueScaling) / log(rMax / rMin);

   % Plot scaling error
    figure(1); 
    plot(trueScaling, 10^3 * (trueScaling - estimatedScaling), 'ok');
    xlim([0.6 1.8]);
    ylim([-30 40])
    axis square
    grid on
    set(gcf, 'Color', 'white');
    set(gca, 'FontSize', 16)
    title1Eval = ['title({''$\sigma_{theory} =  (\frac{R}{r_o}) ^' ...
        '{\Delta x / (1 - n_r) }$''', ...
        ', ''$128 \times 128$ with ' num2str(windowFractions(k) * 100, '%02.0f') '\% window, $d_p = ' ...
        num2str(particleDiameter, '%0.2f') '$''}' , ...
        ', ''interpreter'', ''latex'', ''fontsize'', 16)'];
    eval(title1Eval);
    xlabel('True scaling factor $\sigma_{true}$',...
        'interpreter', 'latex', 'FontSize', 16)

    ylabel('$(\sigma_{true} - \sigma_{theory}) \times 10^3$',...
        'interpreter', 'latex', 'FontSize', 16);
    
    % Plot tx error
    figure(2);
    plot(fxTheory, fxTheory - fmiTranslationX(:, 1), 'ok', 'MarkerFaceColor', 'green');
    xlim([-25 25]);
    ylim([-2 2])
    axis square
    grid on
    set(gcf, 'Color', 'White');
    set(gca, 'FontSize', 16);
    title2Eval = ['title({''$\Delta x_{theory} = (1 - n_r) \ln (\sigma) / \ln( \frac{R}{r_o})$''', ...
        ', ''$128 \times 128$ with ' num2str(windowFractions(k) * 100, '%02.0f') '\% window, $d_p = ' ...
        num2str(particleDiameter, '%0.2f') '$''}' , ...
        ', ''interpreter'', ''latex'', ''fontsize'', 16)'];
    eval(title2Eval);
    xlabel('$\Delta x_{theory}$ (pixels)',...
        'interpreter', 'latex', 'FontSize', 16)
    ylabel('$\Delta x_{theory} - \Delta x_{measured} $ (pixels)',...
        'interpreter', 'latex', 'fontsize', 16);

     plotName1 = ['scalingError_win' ...
    num2str(imageHeight) 'x' num2str(imageWidth)...
    '_win' num2str(windowFractions(k), '%0.2f')...
    '_dp' num2str(particleDiameter, '%0.2f') '.eps'];

    plotName2 = ['txError_win' ...
    num2str(imageHeight) 'x' num2str(imageWidth)...
    '_win' num2str(windowFractions(k), '%0.2f')...
    '_dp' num2str(particleDiameter, '%0.2f') '.eps'];

    print(1, '-depsc', fullfile(plotDir, 'scaling', plotName1));
    print(2, '-depsc', fullfile(plotDir, 'tx', plotName2));
     
    
end

    varsName = ['constantVals_win' ...
    num2str(imageHeight) 'x' num2str(imageWidth)...
    '_dp' num2str(particleDiameter, '%0.2f') '.mat'];

    save(fullfile(plotDir, varsName), 'a', 'windowFractions');

end


function a = const_finder(trueScaling, fmiTranslationX, rMin, rMax, nr);

fxTheory_init =  (1 - nr) * log(trueScaling) / log(rMax/rMin);

val = abs(fxTheory_init) < 10;

% a = lsqnonlin(@const_function, 1.2, [], [], [], ...
%     trueScaling(val), nr, rMax, rMin, fmiTranslationX(val));

% a = lsqnonlin(@const_function, 1.2, [], [], [], ...
%     trueScaling(val), nr, rMax, rMin, fmiTranslationX(val));

a = lsqnonlin(@const_function, 1, [], [], [], ...
    trueScaling(val), nr, rMax, rMin, fmiTranslationX(val));

end


function F = const_function(a, s, nr, rMax, rMin, fx)

% F = (a - nr) * log(s) / log(rMax/rMin) - fx;

% F = exp(log(rMax/rMin) * fx / (a - nr) ) - s;
F = 1 * exp(log(rMax/rMin) * fx / (a - nr) ) - s;

% F = 1 ./ exp(log(rMax/rMin) * fx / (a - nr)) - s;


end










