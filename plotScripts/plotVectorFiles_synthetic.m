% plotVectorFiles

% Close figures
close all

paths{1} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/rpc_prana/128x128/cstep_05/lambvortex_h1024_w1024_rpc_prana_grid32x32_region128x128_000010_000015.mat';

% RPC IID Gradient
paths{2} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/rpc_converged_prana_deform/128x128/cstep_05/lambvortex_h1024_w1024_rpc_converged_prana_deform_grid32x32_region128x128_000010_000015.mat';

% FMC DWO Gradient
paths{3} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/fmc/128x128/cstep_05/lambvortex_h1024_w1024_fmc_grid32x32_region128x128_000010_000015.mat';

% FMC DWO Direct
paths{4} = paths{3};

% FMC IID Gradient
paths{5} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/fmc_deform_disparity/128x128/cstep_05/lambvortex_h1024_w1024_fmc_deform_disparity_grid32x32_region128x128_000010_000015.mat';

% Number of vector files.
nFiles = length(paths);

% Load color map
load('~/Documents/MATLAB/mycmap.mat');

% Which passes to plot.
Pass = [3, 10, 3, 3, 4];

% Vector scaling
Scale = ones(1, nFiles); 

% Vector skips
Skip = 1;

line_width = 0.1;

% Subplot axis limits
xl = 512 + [-200, 200];
yl1 = 296 + [-200, 200];
yl2 = 700 + [-200, 200];

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
kinds = {'grad', 'grad', 'grad', 'dir', 'grad'};

% Plot colors
q_color = 0 * [1/255, 1/255, 1/255];
c_axis = [-1 1];

% Horizontal gap
g = 0.085;

% Horiziontal positions of plots
xpos = 0.1 : g : 0.1 + g * (nFiles-1);

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
    p(2) = 0.49;
    set(gca, 'Position', p);
    xlabel(' ');
    ylabel(' ');
    
    set(gca, 'FontSize', fSize);
    
    set(gca, 'xtick', 400 : 100 : 700);
    
end


print(1, '-depsc', '~/Desktop/test.eps');

!open ~/Desktop/test.eps
% print(1, '-depsc', '/Users/matthewgiarra/Documents/School/VT/Research/Manuscripts/FMIRPC/text/latex/figures/experimental_vorticity_comparison_plot_raw.eps');











