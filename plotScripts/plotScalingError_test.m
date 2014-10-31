function plotScalingError()

JobList = fmcJobList_mc;

JobFile = JobList(1);

windowSizes = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,...
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

nJobs = length(windowSizes);

saveDir = '~/Desktop/scalingErrorTests/mat';

meanError = zeros(nJobs, 1);
stdError = zeros(nJobs, 1);
windowFraction = zeros(nJobs, 1);

imageHeight = JobFile.Parameters.RegionHeight;
imageWidth = JobFile.Parameters.RegionWidth;

for n = 1 : nJobs
  
   fileName = ['scalingPeakError_region' num2str(imageHeight) 'x' num2str(imageWidth) '_win' num2str(windowSizes(n), '%0.2f') 'x' num2str(windowSizes(n), '%0.2f')];
   load(fullfile(saveDir, [fileName '.mat']));
   meanError(n) = mean(abs(fxTheoryError));
   stdError(n) = std(fxTheoryError);
   windowFraction(n) = spatialWindowFraction(1);
    
end

% errorbar(windowFraction, meanError, 1.96 * stdError, 'ok');
loglog(windowFraction, meanError, 'ok')
hold on;
loglog(windowFraction, meanError, '--k');
hold off


end
