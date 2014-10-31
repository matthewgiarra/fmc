% plotVectorFiles

% Font size
fSize = 16;

vectDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/stat';

% RPC DWO
vectorFiles{1} = 'errorStats_rpc_prana_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat'; 

% RPC prana deform
% vectorFiles{2} = 'errorStats_rpc_converged_prana_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';
vectorFiles{2} = 'errorStats_rpc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% FMC DWO
vectorFiles{3} = 'errorStats_fmc_converged_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% FMC deform
% vectorFiles{4} = 'errorStats_fmc_deform_disparity_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';
vectorFiles{4} = 'errorStats_fmc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% Count vector files
nFiles = length(vectorFiles);

% Create the paths to the vector files
for k = 1 : nFiles
   filePaths{k} = fullfile(vectDir, vectorFiles{k});   
end

% Initialize error metrics
mean_err = zeros(nFiles, 3);
std_err = zeros(nFiles, 3);

% Calculate the error statistics.
for k = 1 : nFiles
    
   % Load the vector file.
   load(filePaths{k});   
   
   % Count the number of passes
   nPasses = length(X);
   
   % Calculate the mean and std dev error for each pass
   for p = 1 : nPasses
      mean_err(k, p) =  mean(velocityErrorMagnitudeInside{p}(:));
      std_err(k, p) = std(velocityErrorMagnitudeInside{p}(:));
   end
   
end

rpc_conv = find(abs(diff(mean_err(2, :))) < 0.01, 1);
fmc_conv = find(abs(diff(mean_err(4, :))) < 0.01, 1);

% Force first passes of FMC DWO and IID to be equivalent because
% they must be, 
mean_err(3, 1) = mean_err(4, 1);
std_err(3, 1) = std_err(4, 1);

% Markers
marks = {'ok'; 'ok'; 'sk'; 'sk'};
markColors = {'black'; 'white'; 'black'; 'white'};

% Plot the error statistics
figure(1);
hold on
for k = 1 : nFiles
    errorbar((k-1)*4 + 1 : (k-1)*4 + 3, mean_err(k, 1:3), std_err(k, 1:3), marks{k}, 'markerfacecolor', markColors{k}, 'MarkerSize', 10);
    
end

% Plot converged errors
% plot([0 15], mean_err(2, 10) * [1 1], '--k');
% plot([0 15], mean_err(4, 7) * [1 1], '--b');



% patch([-1, 16, 16, -1], mean_err(2, rpc_conv) + [-std_err(2,rpc_conv), -std_err(2,rpc_conv), std_err(2,rpc_conv), std_err(2,rpc_conv)], 180/255 * [1 1 1]);
% patch([-1, 16, 16, -1], mean_err(4, fmc_conv) + [-std_err(4, fmc_conv), -std_err(4, fmc_conv), std_err(4, fmc_conv), std_err(4, fmc_conv)], 100/255 * [ 1 1 1]);


hold off
axis square
grid on

% Set plot limites
xlim([0 16]);
ylim([-5 20]);

% Set plot ticks
xt = [1:3, 5:7, 9:11, 13:15];
xtl = {'1' '2' '3' '1' '2' '3' '1' '2' '3' '1' '2' '3'};
set(gca, 'xtick', xt);
set(gca, 'xticklabel', xtl);

set(gca, 'ytick', [0, 5, 10, 15, 20]);

% % Create plot legend
% h = legend('RPC DWO', 'RPC IID', 'FMC DWO', 'FMC IID', 'RPC IID (converged)', 'FMC IID (converged)'); 
% Create plot legend
h = legend('RPC DWO', 'RPC IID', 'FMC DWO', 'FMC IID'); 

% Labels
xlabel('Iteration number', 'FontSize', fSize);
ylabel('Displacement error magnitude (pixels)', 'FontSize', fSize);
title({'Average displacement error magnitude', 'near vortex cores (r < r_c)'}, 'Interpreter', 'tex', 'FontSize', fSize);
set(gca, 'FontSize', fSize);
% Enable box
box on;

figure(2);
errorbar(1, mean_err(2, 3), std_err(2, 3), marks{2}, 'markerfacecolor', markColors{2}, 'markersize', 10);
hold on
errorbar(2, mean_err(2, 8), std_err(2, 8), marks{2}, 'markerfacecolor', markColors{2},'MarkerSize', 10);
errorbar(3, mean_err(3, 3), std_err(3,3), marks{3}, 'markerfacecolor', markColors{3}, 'MarkerSize', 10);
errorbar(4, mean_err(4, 3), std_err(4, 3), marks{4}, 'markerfacecolor', markColors{4},'MarkerSize', 10);
set(gca, 'ytick', 0:0.25:2.0);
set(gca, 'xtick', 1:4);

g = get(gca, 'ytick');

for k = 1 : length(g)
    ytl{k} = num2str(g(k), '%0.2f');
end
set(gca, 'yticklabel', ytl);
set(gca, 'xticklabel', {'3', '8(c)', '3', '3(c)'})
set(gca, 'plotboxaspectratio', [0.3 1 1]);
xlim([0.5 4.5]);
ylim([-1.4 1.9])
grid on
box on;
set(gca, 'fontsize', fSize);








