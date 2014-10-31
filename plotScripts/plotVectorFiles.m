% plotVectorFiles

load('~/Documents/MATLAB/mycmap.mat');
startImage = 850;
endImage = 949;
imSkip = 1;
nDigits = 6;
numFormat = ['%0' num2str(nDigits) '.0f'];

correlationStep = 5;
corrMethod = 'rpc';
regionHeight = 128;
regionWidth = 128;

imageNumbers = startImage : imSkip : endImage;
nImages = length(imageNumbers);

regionStr = [num2str(regionHeight) 'x' num2str(regionWidth)];
corrString = ['cstep_' num2str(correlationStep, '%02.0f')];

vectBaseDir = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect';
vectDir = fullfile(vectBaseDir, corrMethod, regionStr, corrString);

vectBase = 'vortexring_d03_f60_t06_';
vectFullBase = [vectBase corrMethod '_grid32x32_region128x128_'];

Pass = 3;
% Scale = 2;
yStart = 40;
yEnd = 5;
xStart = 40;
xEnd = 80;

% xl(1, :) = 512 + [-200 200];
% yl(1, :) = 408 + [-200 200];
% 
% xl(2, :) = xl(1, :);
% yl(2, :) = yl(1, :);
% 
% yl(3, :) = 408 + [-75 75];
% xl(3, :) = 535 + [-75 75];

% Scale = [3.5, 1, 1.5];
Scale = [1, 1, 1];

xl = [430 630];
yl1 = [330 530];
yl2 = [700, 900];

% xl = 512 + [-200, 200];
% yl1 = 296 + [-200, 200];
% yl2 = 700 + [-200, 200];


gap_h = 0.001;
gap_w = 0.1;
marg_lower = 0.1;
marg_upper = 0.1;
marg_left = 0.1;
marg_right = 0.1;

close all
ha = tight_subplot(2, 3, [gap_h, gap_w], [marg_lower, marg_upper], [marg_left, marg_right]);
% load('~/Desktop/paths.mat');
% paths{1} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/fmc/128x128/cstep_06/lambvortex_h1024_w1024_fmc_grid32x32_region128x128_000001_000007.mat';
% paths{2} = paths{1};
% paths{3} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/synthetic/vortex/lambvortex_2014-05-22_centralDifference/c_0.0250/lambvortex_h1024_w1024_00001/vect/rpc/128x128/cstep_06/lambvortex_h1024_w1024_rpc_grid32x32_region128x128_000001_000007.mat';


paths{1} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/fmc/128x128/cstep_10/vortexring_d03_f60_t06_fmc_grid32x32_region128x128_000913_000923.mat';
paths{2} = paths{1};
paths{3} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/rpc/128x128/cstep_10/vortexring_d03_f60_t06_rpc_grid32x32_region128x128_000913_000923.mat';



for k = 1 : 3
    load(paths{k});
    x{k} = X{Pass};
    y{k} = Y{Pass};
    u{k} = U{Pass};
    v{k} = V{Pass};
    r{k} = R{Pass};
    uVal{k} = UVAL{Pass};
    vVal{k} = VVAL{Pass};
    rVal{k} = RVAL{Pass};
    
    if k == 1
        rGrad{k} = calculateVorticity_socdiff(x{k}, y{k}, uVal{k}, vVal{k}, rVal{k}, 0);
    else
        rGrad{k} = calculateVorticity_socdiff(x{k}, y{k}, uVal{k}, vVal{k}, rVal{k}, 0);
    end
    
    rDir{k}  = 2 * r{k};
    
end

kinds = {'grad', 'dir', 'grad'};

Skip = 1;
% Scale = 0.8;

q_color = 0 * [1/255, 1/255, 1/255];
c_axis = [-1 1];

xpos = [0.1, 0.32, 0.54];

for k = 1:3
    
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
    
    if isGrad
        q{k} = quiver(x{k}(1:Skip:end, 1:Skip:end), ...
        y{k}(1:Skip:end, 1:Skip:end), ...
        Scale(k) * uVal{k}(1:Skip:end, 1:Skip:end), ...
        Scale(k) * vVal{k}(1:Skip:end, 1:Skip:end), 0);
    else
        q{k} = quiver(x{k}(1:Skip:end, 1:Skip:end), ...
        y{k}(1:Skip:end, 1:Skip:end), ...
        Scale(k) * u{k}(1:Skip:end, 1:Skip:end), ...
        Scale(k) * v{k}(1:Skip:end, 1:Skip:end), 0);
    end
    
    set(q{k}, 'Color', q_color);
    set(gca, 'xticklabel', []);
    if k > 1
        set(gca, 'yticklabel', []);
    end
    
    
    p = get(gca, 'Position');
    p(1) = xpos(k);
    p(2) = 0.39;
    
    set(gca, 'Position', p);
    
end

for k = 4:6
    
    axes(ha(k));
%     subplot(2, 3, k);

    isGrad = ~isempty(regexpi(kinds{k-3}, 'grad'));

    if isGrad
         im{k} = imagesc(x{k-3}(:), y{k-3}(:), rGrad{k-3});        
    else
         im{k} = imagesc(x{k-3}(:), y{k-3}(:), 2*r{k-3});
    end
    
    axis image;
    colormap(mycmap);
    caxis(c_axis);
    xlim(xl);
    ylim(yl2);
    hold on
    
    if isGrad
    q{k} = quiver(x{k-3}(1:Skip:end, 1:Skip:end), ...
        y{k-3}(1:Skip:end, 1:Skip:end), ...
        Scale(k-3) * uVal{k-3}(1:Skip:end, 1:Skip:end), ...
        Scale(k-3) * vVal{k-3}(1:Skip:end, 1:Skip:end), 0);
    else
        q{k} = quiver(x{k-3}(1:Skip:end, 1:Skip:end), ...
        y{k-3}(1:Skip:end, 1:Skip:end), ...
        Scale(k-3) * u{k-3}(1:Skip:end, 1:Skip:end), ...
        Scale(k-3) * v{k-3}(1:Skip:end, 1:Skip:end), 0);
    end
    
    set(q{k}, 'Color', q_color);
    
    if k > 4
        set(gca, 'yticklabel', []);
    end
    
    p = get(gca, 'Position');
    p(1) = xpos(k-3);
    set(gca, 'Position', p);
    xlabel(' ');
    ylabel(' ');
    
end

% print(1, '-depsc', '/Users/matthewgiarra/Documents/School/VT/Research/Manuscripts/FMIRPC/text/latex/figures/experimental_vorticity_comparison_plot_raw.eps');











