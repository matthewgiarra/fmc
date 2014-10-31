% plotVectorFiles

load('~/Documents/MATLAB/mycmap.mat');

Skip = 1;

dx = 100;
dy = 100;
% xl_01 = 512 + [-200, 200];
xl = 525 + dx * [-1, 1];
% yl_01 = 296 + [-200, 200];
% yl_01 = 736 + [-200, 200];
yl_02 = 800 + dy * [-1, 1];
yl_01 = 425 + dy * [-1, 1];

% xl_02 = 512 + [-75, 75];
% xl_02 = 525 + [-75, 75];
% yl_02 = 296 + [-75, 75];
%  yl_02 = 736 + [-75, 75];
%  yl_02 = 816 + [-75, 75];
% yl_02 = 416 + [-75, 75];
 
gap_h = 0.01;
gap_w = 0.01;
marg_lower = 0.2; 
marg_upper = 0.1;
marg_left = 0.2;
marg_right = 0.2;

px = 0.43;
py = 0.03;
fSize = 15;

Scale = 3;
Skip = 2;
close all
ha = tight_subplot(2, 4, [gap_h, gap_w], [marg_lower, marg_upper], [marg_left, marg_right]);
q_color = 'black';

c_lim = 0.5;

% RPC DWO
paths{1} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/rpc/128x128/cstep_05/vortexring_d03_f60_t06_rpc_grid32x32_region128x128_000913_000918.mat';

% RPC Deform
paths{2} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/rpc_prana_deform/128x128/cstep_05/vortexring_d03_f60_t06_rpc_prana_deform_grid32x32_region128x128_000913_000918.mat';

% FMC
paths{3} = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/FMISPOMF/analysis/data/experimental/vortex/vortexring_2013-11-12_d03_f60_t06/vect/fmc/128x128/cstep_05/vortexring_d03_f60_t06_fmc_grid32x32_region128x128_000913_000918.mat';


% RPC DWO
load(paths{1});
rGrad_rpc = calculateVorticity_socdiff(X{3}, Y{3}, UVAL{3}, VVAL{3}, RVAL{3}, 0);

axes(ha(1));
imagesc(X{3}(:), Y{3}(:), rGrad_rpc);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_01);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
title('RPC DWO Gradient', 'FontSize', 8);
caxis(c_lim*[-1 1]);

axes(ha(5));
imagesc(X{3}(:), Y{3}(:), rGrad_rpc);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_02);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
p = get(gca, 'position');
p(2) = 0.35;
set(gca, 'position', p);
caxis(c_lim*[-1 1]);

% RPC Deform
load(paths{2});
rGrad_rpc_deform = calculateVorticity_socdiff(X{3}, Y{3}, UVAL{3}, VVAL{3}, RVAL{3}, 0);

axes(ha(2));
imagesc(X{3}(:), Y{3}(:), rGrad_rpc_deform);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_01);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
title('RPC CWO Gradient', 'FontSize', 8);
caxis(c_lim*[-1 1]);

axes(ha(6));
imagesc(X{3}(:), Y{3}(:), rGrad_rpc_deform);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_02);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
p = get(gca, 'position');
p(2) = 0.35;
set(gca, 'position', p);
caxis(c_lim*[-1 1]);


% FMC DWO
load(paths{3});
rGrad_fmc = calculateVorticity_socdiff(X{3}, Y{3}, UVAL{3}, VVAL{3}, RVAL{3}, 0);

axes(ha(3));
imagesc(X{3}(:), Y{3}(:), rGrad_fmc);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_01);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
title('FMC DWO Gradient', 'FontSize', 8);
caxis(c_lim*[-1 1]);

axes(ha(7));
imagesc(X{3}(:), Y{3}(:), rGrad_fmc);
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * UVAL{3}(1:Skip:end, 1:Skip:end), Scale * VVAL{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_02);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
p = get(gca, 'position');
p(2) = 0.35;
set(gca, 'position', p);
caxis(c_lim*[-1 1]);

axes(ha(4));
imagesc(X{3}(:), Y{3}(:), 2*R{3});
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * U{3}(1:Skip:end, 1:Skip:end), Scale * V{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_01);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
title('FMC Direct', 'FontSize', 8);
caxis(c_lim*[-1 1]);


axes(ha(8));
imagesc(X{3}(:), Y{3}(:), 2*R{3});
axis image;
colormap(mycmap);
hold on;
quiver(X{3}(1:Skip:end, 1:Skip:end), Y{3}(1:Skip:end, 1:Skip:end), Scale * U{3}(1:Skip:end, 1:Skip:end), Scale * V{3}(1:Skip:end, 1:Skip:end), 0, 'black'); axis image;
xlim(xl);
ylim(yl_02);
% plot([xl_02(1) xl_02(2) xl_02(2) xl_02(1) xl_02(1)], ...
%     [yl_02(1) yl_02(1) yl_02(2) yl_02(2) yl_02(1)], ...
%     '--k', 'linewidth', 2);
hold off
set(gca, 'FontSize', fSize);
xlabel(' '); ylabel(' ');
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
p = get(gca, 'position');
p(2) = 0.35;
set(gca, 'position', p);
caxis(c_lim*[-1 1]);


% print(1, '-depsc', '~/Desktop/vorticity_compare_test.eps');




