vectDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/stat';
cSteps = 1:8;

corrTypes = {'rpc', 'fmc'};
startSet = 1;
endSet = 100;

nSteps = length(cSteps);

Pass = 1;

meanPeakHeightRatio = zeros(nSteps, 2);
stdPeakHeightRatio = zeros(nSteps, 2);
for c = 1 : 2
    for k = 1 : nSteps
       fileName =  ['errorStats_' corrTypes{c} '_lambvortex_h1024_w1024_c_0.0250_cstep_' num2str(cSteps(k), '%02.0f') '_sets_' num2str(startSet, '%05.0f') '-' num2str(endSet, '%05.0f') '.mat'];
       filePath = fullfile(vectDir, fileName);
       load(filePath);
       meanPeakHeightRatio(k, c) = mean(spatialPeakRatioInside{Pass}(:));
       stdPeakHeightRatio(k, c) = std(spatialPeakRatioInside{Pass}(:)); 
    end
end

hold off
plot(cSteps, meanPeakHeightRatio(:, 1),...
    '-ok', 'markerfacecolor', 'black', 'markersize', 10);
hold on
plot(cSteps, meanPeakHeightRatio(:, 2),...
    '-sk', 'markerfacecolor', 'white', 'markersize', 10);
plot([0 cSteps(end)+1], [1 1], '--k');
hold off
grid on
axis square

set(gca, 'ytick', [0 1 5 10 15 20 25]);
