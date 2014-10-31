% plotVectorFiles

% Close figures
close all

c_step = 5;

start_image = 860;

% Case name
case_name = 'vortexring_d03_f60_t06_';

% Base vector directory
vect_dir = '~/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/';

% RPC DWO Gradient
paths{1} = fullfile(vect_dir, 'rpc', '128x128', ['cstep_' num2str(c_step, '%02.0f')],  [case_name 'rpc_grid32x32_region128x128_' num2str(start_image, '%06.0f') '_' num2str(start_image + c_step, '%06.0f') '.mat']);

% RPC IID Gradient
paths{2} = fullfile(vect_dir, 'rpc_deform', '128x128', ['cstep_' num2str(c_step, '%02.0f')],  [case_name 'rpc_deform_grid32x32_region128x128_' num2str(start_image, '%06.0f') '_' num2str(start_image + c_step, '%06.0f') '.mat']);

% FMC DWO Gradient
paths{3} = fullfile(vect_dir, 'fmc', '128x128', ['cstep_' num2str(c_step, '%02.0f')],  [case_name 'fmc_grid32x32_region128x128_' num2str(start_image, '%06.0f') '_' num2str(start_image + c_step, '%06.0f') '.mat']);

% FMC IID Gradient
paths{4} = fullfile(vect_dir, 'fmc_deform', '128x128', ['cstep_' num2str(c_step, '%02.0f')],  [case_name 'fmc_deform_grid32x32_region128x128_' num2str(start_image, '%06.0f') '_' num2str(start_image + c_step, '%06.0f') '.mat']);

% FMC DWO Direct
paths{5} = paths{3};


% Number of vector files.
nFiles = length(paths);

% Load color map
load('~/Documents/MATLAB/mycmap.mat');

% Which passes to plot.
Pass = [3, 10, 3, 7, 3];

% Vector scaling
Scale = 2 * ones(1, nFiles); 

% Vector skips
Skip = 1;

line_width = 0.25;

% Subplot axis limits
xl = 480 + [-100, 100];
yl1 = 416 + [-100, 100];
yl2 = 816 + [-100, 100];

% Plot parameters
gap_h = 0.001;
gap_w = 0.1;
marg_lower = 0.1;
marg_upper = 0.1;
marg_left = 0.1;
marg_right = 0.1;

% Make the subplot figure.
ha = tight_subplot(2, nFiles, [gap_h, gap_w], [marg_lower, marg_upper], [marg_left, marg_right]);

% Vorticity methods
kinds = {'grad', 'grad', 'grad', 'grad', 'dir'};
titles = {'RPC DWO Gradient'; 'RPC IID Gradient'; 'FMC DWO Gradient'; 'FMC IID Gradient'; 'FMC DWO Direct'};

% Plot colors
q_color = 0 * [1/255, 1/255, 1/255];
c_axis = [-0.5 0.5];

% Horizontal gap
g_h = 0.085;

% Horiziontal positions of plots
xpos = 0.1 : g_h : 0.1 + g_h * (nFiles-1);

quiver_skip = 2;
color_skip = 1;
fSize = 4.5;

for k = 1 : nFiles
    load(paths{k});
    x{k} = X{Pass(k)};
    y{k} = Y{Pass(k)};
    u{k} = U{Pass(k)};
    v{k} = V{Pass(k)};
    r{k} = R{Pass(k)};
    uVal{k} = UVAL{Pass(k)};
    vVal{k} = VVAL{Pass(k)};
    rVal{k} = RVAL{Pass(k)};
    
    % Gradient based vorticity
    rGrad{k} = calculateVorticity_socdiff(x{k}, y{k}, uVal{k}, vVal{k}, rVal{k}, 0);
    
    % Direct vorticity
    rDir{k}  = 2 * r{k};
    
end


for k = 1 : nFiles
    
    isGrad = ~isempty(regexpi(kinds{k}, 'grad'));
    
    axes(ha(k));
%     subplot(2, 3, k)
    
    if isGrad;
        im{k} = imagesc(x{k}(:), y{k}(:), rGrad{k});     
    else
         im{k} = imagesc(x{k}(:), y{k}(:), 2*r{k});
    end
    
    axis image;
    colormap(mycmap);
    caxis(c_axis);
    xlim(xl);
    ylim(yl1);
    hold on
    
    q{k} = quiver(x{k}(1:quiver_skip:end, 1:quiver_skip:end), ...
    y{k}(1:quiver_skip:end, 1:quiver_skip:end), ...
    Scale(k) * u{k}(1:quiver_skip:end, 1:quiver_skip:end), ...
    Scale(k) * v{k}(1:quiver_skip:end, 1:quiver_skip:end), 0);
    
    title(titles{k}, 'FontSize', fSize);

%     if isGrad
%         q{k} = quiver(x{k}(1:Skip:end, 1:Skip:end), ...
%         y{k}(1:Skip:end, 1:Skip:end), ...
%         Scale(k) * uVal{k}(1:Skip:end, 1:Skip:end), ...
%         Scale(k) * vVal{k}(1:Skip:end, 1:Skip:end), 0);
%     else
%         q{k} = quiver(x{k}(1:Skip:end, 1:Skip:end), ...
%         y{k}(1:Skip:end, 1:Skip:end), ...
%         Scale(k) * u{k}(1:Skip:end, 1:Skip:end), ...
%         Scale(k) * v{k}(1:Skip:end, 1:Skip:end), 0);
%     end
    
    set(q{k}, 'Color', q_color);
    set(q{k}, 'LineWidth', line_width);
    set(gca, 'xticklabel', []);
    if k > 1
        set(gca, 'yticklabel', []);
    end
    
    
    p = get(gca, 'Position');
    p(1) = xpos(k);
    p(2) = 0.6;
    set(gca, 'Position', p);
    
    set(gca, 'FontSize', fSize);
    
    
end

for k = nFiles + 1 : 2 * nFiles
    
    axes(ha(k));
%     subplot(2, 3, k);

    isGrad = ~isempty(regexpi(kinds{k - nFiles}, 'grad'));

    if isGrad
         im{k} = imagesc(x{k - nFiles}(:), y{k - nFiles}(:), rGrad{k - nFiles});        
    else
         im{k} = imagesc(x{k - nFiles}(:), y{k - nFiles}(:), 2*r{k - nFiles});
    end
    
    axis image;
    colormap(mycmap);
    caxis(c_axis);
    xlim(xl);
    ylim(yl2);
    hold on
    
    if isGrad
    q{k} = quiver(x{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        y{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        Scale(k - nFiles) * uVal{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        Scale(k - nFiles) * vVal{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), 0);
    else
        q{k} = quiver(x{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        y{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        Scale(k - nFiles) * u{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), ...
        Scale(k - nFiles) * v{k - nFiles}(1:quiver_skip:end, 1:quiver_skip:end), 0);
    end
    
    set(q{k}, 'Color', q_color);
    set(q{k}, 'LineWidth', line_width);
    
    if k > nFiles + 1
        set(gca, 'yticklabel', []);
    end
    
    p = get(gca, 'Position');
    p(1) = xpos(k - nFiles);
    p(2) = 0.488;
    set(gca, 'Position', p);
    xlabel(' ');
    ylabel(' ');
    
    set(gca, 'FontSize', fSize);
    
    set(gca, 'xtick', 350 : 50 : 550);
    
end


print(1, '-depsc', '~/Desktop/test.eps');

!open ~/Desktop/test.eps
% print(1, '-depsc', '/Users/matthewgiarra/Documents/School/VT/Research/Manuscripts/FMIRPC/text/latex/figures/experimental_vorticity_comparison_plot_raw.eps');











