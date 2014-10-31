nPasses = length(allErrInsideFmc);
saveDir = '~/Desktop/cdfs';


for p = 1 : nPasses
hold off 
figure(1);
cdfplot(allErrInsideFmc{p});
hold on
cdfplot(allErrInsideRpc{p});

axis square
xlim([0 1]);
ylim([0 1]);

h = get(gca, 'children');

set(h(2), 'color', 'black');
set(h(6), 'color', 'black');
set(h(6), 'linestyle', '--');
legend('FMC', 'RPC');
xlabel('Velocity error magnitude (pix');
ylabel('Cumulative probability');
title(['Velocity error near core, Pass ' num2str(p)]);

saveName = ['error_cdf_inside_pass_' num2str(p, '%02.0f') '.eps'];

savePath = fullfile(saveDir, saveName);

print(1, '-depsc', savePath);

hold off

end
