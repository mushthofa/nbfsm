function plotAccuracyDist(savedir, drug, shufAcc, nofsAcc, topfAcc, fsAcc, diffAcc)
% Plot the RF accuracies using all features, top statistical features,
% network-based feature selection and diffusion against the background
% distribution of accuracies when using shuffled data

[~, mudd, sigdd] = zscore(shufAcc);

figure;
histfit(shufAcc, 20);
hold on;
g=gcf;
ax=g.CurrentAxes;
ylim = get(ax, 'Ylim' ); 

plot([nofsAcc nofsAcc], ylim, '--r', 'linewidth', 2);
plot([topfAcc topfAcc], ylim, '--b', 'linewidth', 2);
plot([fsAcc fsAcc], ylim,'--g', 'linewidth', 2);
plot([diffAcc diffAcc], ylim,'--y', 'linewidth', 2);
set(ax, 'XLim', [0 1]);
text(0.1, 0.9*ylim(2), sprintf('Shuffled data RF accuracy \\mu=%.2f,\\sigma=%.2f', mudd, sigdd));
title(sprintf('%s - Accuracies of random forests (RFs) against data RF with shuffled data', drug));
legend('RF Accuracies from shuffled data (100 repeats)', 'Normal fit', 'Random forest with all features',  'Random forest with top features', 'Random forest with network feature selection', 'Random forest with feature selection+diffusion');
g.Position = [0 0 1920 1080];

filename = sprintf('%saccshufdist/%s.png', savedir, matlab.lang.makeValidName(drug));
print(filename, '-dpng');
close all;
    
end