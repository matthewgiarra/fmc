% plotVectorFiles
close all
% Font size
fSize = 16;

vectDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/stat';

% RPC prana deform
% vectorFiles{1} = 'errorStats_rpc_converged_prana_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';
% vectorFiles{1} = 'errorStats_rpc_deform_disparity_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';
vectorFiles{1} = 'errorStats_rpc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';


% FMC deform
% vectorFiles{2} = 'errorStats_fmc_deform_disparity_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';
vectorFiles{2} = 'errorStats_fmc_deform_lambvortex_h1024_w1024_c_0.0250_cstep_05_sets_00001-00100.mat';

% Count vector files
nFiles = length(vectorFiles);

% Create the paths to the vector files
for k = 1 : nFiles
   filePaths{k} = fullfile(vectDir, vectorFiles{k});   
end

% Line styles
colors_01 = flipud({'black'; 'red'; 'blue'; 0.5 * [1 1 1]; 'black'; 'red'});
styles_01 = flipud({'-'; '-'; '-'; '-'; '--'; '--'});

colors_02 = flipud({'black'; 'red'; 'blue'; 0.5 * [1 1 1]});
styles_02 = flipud({'-'; '-'; '-'; '-'});

% Pass numbers for each processing method
rpc_pass_numbers = [1, 2, 3, 4, 7, 8];
fmc_pass_numbers = [1, 2, 3, 4];

% Initialize legend entry string structures
rpc_legend_entries{length(rpc_pass_numbers)} = '';
fmc_legend_entries{length(fmc_pass_numbers)} = '';

% Load the first data file
rpc_data = load(filePaths{1});
figure(1);
ha = tight_subplot(1, 2);
axes(ha(1));
for k = 1 : length(rpc_pass_numbers)
    rpc_legend_entries{k} = ['Pass ' num2str(rpc_pass_numbers(k))];
    cdfplot(rpc_data.velocityErrorMagnitudeInside{rpc_pass_numbers(k)}(:));
    hold on
end

hold off
axis square
xlim([0 1]);

% Flip children for labeling
h = get(gca, 'children');

% Specify the line colors and styles
for k = 1 : length(rpc_pass_numbers)
    set(h(k), 'color', colors_01{k});
    set(h(k), 'linestyle', styles_01{k});
    set(h(k), 'LineWidth', 2);
    
end

% Make the rpc legend
g = legend(rpc_legend_entries, 'Location', 'Southeast');
set(gca, 'xtick', 0 : 0.2 : 0.8);
set(gca, 'ytick', 0 : 0.2 : 1);

%%%% FMC

% Load the second data file
fmc_data = load(filePaths{2});
axes(ha(2));
for k = 1 : length(fmc_pass_numbers)
    fmc_legend_entries{k} = ['Pass ' num2str(fmc_pass_numbers(k))];
    cdfplot(fmc_data.velocityErrorMagnitudeInside{fmc_pass_numbers(k)}(:));
    hold on
end

hold off
axis square
xlim([0 1]);

% Flip children for labeling
h = get(gca, 'children');

% Specify the line colors and styles
for k = 1 : length(fmc_pass_numbers)
    set(h(k), 'color', colors_02{k});
    set(h(k), 'linestyle', styles_02{k});
    set(h(k), 'LineWidth', 2);
end

% set(gca,'Children',flipud(h)) %Invert the order of the objects
% Create legend
g = legend(fmc_legend_entries, 'Location', 'Southeast');
set(gca, 'xtick', 0 : 0.2 : 1);
set(gca, 'ytick', 0 : 0.2 : 1);

% Plot labels
title('');
ylabel('');
set(gca, 'yticklabel', '');
xlabel('');
set(gca, 'FontSize', fSize);

% Plot labels
axes(ha(1));
title('');
ylabel('Cumulative probability', 'FontSize', fSize);
xlabel('');
set(gca, 'FontSize', fSize);

a_xlabel = annotation('textbox', ...
    'position', [0.25 0.05 0.5 0.1],...
    'String', 'Displacement error magnitude (pix)', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'None', ...
    'FontSize', 16);

a_rpc_label = annotation('textbox', ...
    'position', [0.35, 0.705, 0.13, 0.06],...
    'String', 'RPC IID', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'Black', ...
    'backgroundcolor', 'white', ...
    'FontSize', 16);

a_fmc_label = annotation('textbox', ...
    'position', [0.81, 0.705, 0.13, 0.06],...
    'String', 'FMC IID', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'Black', ...
    'backgroundcolor', 'white', ...
    'FontSize', 16);

a_title = annotation('textbox', ...
    'position', [0.15, 0.76, 0.7, 0.1],...
    'String', 'Displacement error magnitude near vortex cores', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'None', ...
    'FontSize', 16);


