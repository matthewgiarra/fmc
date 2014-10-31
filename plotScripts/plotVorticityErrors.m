% plotVectorFiles

% Font size
fSize = 16;

vectDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/stat';

% RPC DWO Gradient
vectorFiles{1} = 'errorStats_rpc_prana_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat'; 

% RPC Deform Gradient
vectorFiles{2} = 'errorStats_rpc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% FMC DWO Gradient
vectorFiles{3} = 'errorStats_fmc_converged_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% FMC Deform Gradient
vectorFiles{4} = 'errorStats_fmc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% FMC DWO Direct
vectorFiles{5} = 'errorStats_fmc_converged_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% Count vector files
nFiles = length(vectorFiles);

% Create the paths to the vector files
for k = 1 : nFiles
   filePaths{k} = fullfile(vectDir, vectorFiles{k});   
end

% Initialize error metrics
mean_err = zeros(nFiles, 3);
std_err = zeros(nFiles, 3);
converged_err_mean = zeros(nFiles, 1);

% Calculate the error statistics.
for k = 1 : nFiles
    
   % Load the vector file.
   load(filePaths{k});   
      
   if k ~=5
   % Calculate the mean and std dev error for each pass
   for p = 1 : length(X)
      mean_err(k, p) =  mean(rErrorGradInside{p}(:));
      std_err(k, p) = std(rErrorGradInside{p}(:));
      converged_err_mean(k) = mean(rErrorGradInside{end}(:)); 
   end
   
   else
       for p = 1 : length(X)
           mean_err(k, p) =  mean(rErrorDirInside{p}(:));
           std_err(k, p) = std(rErrorDirInside{p}(:));
           converged_err_mean(k) = mean(rErrorDirInside{end}(:)); 
       end
   end
   
end

% Markers
marks = {'ok'; 'ok'; 'sk'; 'sk'; 'vk'};
markColors = {'black'; 'white'; 'black'; 'white'; 'black'};

% Plot the error statistics
figure(1);
hold on
for k = 1 : nFiles
    errorbar((k-1)*4 + 1 : (k-1)*4 + 3, mean_err(k, 1:3), std_err(k, 1:3), marks{k}, 'markerfacecolor', markColors{k}, 'MarkerSize', 10);    
end

% % Plot converged errors
% plot([0 20], converged_err_mean(2) * [1 1], '--k');
% plot([0 20], converged_err_mean(4) * [1 1], '--b');

hold off
axis square
grid on

% Set plot limites
xlim([0 20]);
ylim([-0.1 0.6]);

% Set plot ticks
xt = [1:3, 5:7, 9:11, 13:15, 17:19];
xtl = {'1' '2' '3' '1' '2' '3' '1' '2' '3' '1' '2' '3' '1' '2' '3'};
set(gca, 'xtick', xt);
set(gca, 'xticklabel', xtl);

set(gca, 'ytick', 0:0.1:0.5);
g = get(gca, 'ytick');
for k = 1 : length(g)
   ytl_01{k} = num2str(g(k), '%0.1f');
end
set(gca, 'yticklabel', ytl_01);

% set(gca, 'ytick', [0, 5, 10, 15, 20]);

% Create plot legend
h = legend('RPC DWO Gradient', 'RPC IID Gradient', 'FMC DWO Gradient', 'FMC IID Gradient', 'FMC DWO Direct'); 
set(h, 'FontSize', 14);
% Labels
xlabel('Iteration number', 'FontSize', fSize);
ylabel('Vorticity error magnitude (1 / frames)', 'FontSize', fSize);
title({'Average vorticity error magnitude', 'near vortex cores (\omega > 1/2 \omega_{max})'}, 'Interpreter', 'tex', 'FontSize', fSize);
set(gca, 'FontSize', fSize);
% Enable box
box on;

% % Flip the plot order
% c = get(gca,'Children'); %Get the handles for the child objects from the current axes
% set(gca,'Children',flipud(c)) %Invert the order of the objects

figure(2);
errorbar(1, mean_err(2, 3), std_err(2, 3), marks{2}, 'markerfacecolor', markColors{2}, 'MarkerSize', 10);
hold on
errorbar(2, mean_err(2, 8), std_err(2, 8), marks{2}, 'markerfacecolor', markColors{2},'MarkerSize', 10);
errorbar(3, mean_err(3, 3), std_err(3, 3), marks{3}, 'markerfacecolor', markColors{3},'MarkerSize', 10);
errorbar(4, mean_err(4, 3), std_err(4, 3), marks{4}, 'markerfacecolor', markColors{4},'MarkerSize', 10);
set(gca, 'ytick', 0:0.01:0.07);
set(gca, 'xtick', 1:4);

g = get(gca, 'ytick');

for k = 1 : length(g)
    ytl_02{k} = num2str(g(k), '%0.2f');
end
set(gca, 'yticklabel', ytl_02);
set(gca, 'xticklabel', {'3', '8(c)', '3', '3(c)'});
set(gca, 'plotboxaspectratio', [0.30 1 1]);
xlim([0.5 4.5]);
ylim([-0.045 0.070])
grid on
box on;
set(gca, 'fontsize', fSize);




