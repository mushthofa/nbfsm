% Create plots from CV results

% Assuming savedir is already defined
savedir = 'cvres1/';

% Folder to hold all accuracy distribution plots
mkdir(savedir, 'accshufdist'); 
% Folder to hold all predictor matrix plots
mkdir(savedir, 'predictormat');

% Load drug data
load(sprintf('%sdrugdata.mat', savedir));

% Load Accuracy results from the RF with shuffled features ~ drug response
load('shufcv.mat');
shuf_p = shufCVResult(:,:,1)+shufCVResult(:,:,4);
shuf_n = shufCVResult(:,:,2)+shufCVResult(:,:,3);
shufacc = (shufCVResult(:,:,1) + shufCVResult(:, :, 3))./(shuf_p + shuf_n);

% Load cv results
load(sprintf('%scvresult.mat', savedir));


% Plot results of accuracies compared to shuffled data
nofsacc = zeros(height(thaatab), 1);
topfacc = zeros(height(thaatab), 1);
fsacc = zeros(height(thaatab), 1);
diffacc = zeros(height(thaatab), 1);

nofsrec = zeros(height(thaatab), 1);
topfrec = zeros(height(thaatab), 1);
fsrec = zeros(height(thaatab), 1);
diffrec = zeros(height(thaatab), 1);

nofspre = zeros(height(thaatab), 1);
topfpre = zeros(height(thaatab), 1);
fspre = zeros(height(thaatab), 1);
diffpre = zeros(height(thaatab), 1);

pvalsNOFS = nan(height(thaatab), 1);
pvalsTOPF = nan(height(thaatab), 1);
pvalsFS = nan(height(thaatab), 1);

fprintf('Creating plots in the folder %s...\n', savedir);
for i=1:height(thaatab)
    drug = thaatab.DRUG{i};
    
    % Computing accuracies
    [nofsacc(i), nofsrec(i), nofspre(i)] = acc(predictionTable.NOFS(i, :), referenceTable(i, :));
    [topfacc(i), topfrec(i), topfpre(i)] = acc(predictionTable.TOPF(i, :), referenceTable(i, :));
    [fsacc(i), fsrec(i), fspre(i)] = acc(predictionTable.FS(i, :), referenceTable(i, :));
    [diffacc(i), diffrec(i), diffpre(i)] = acc(predictionTable.DIFF(i, :), referenceTable(i, :));
    
    % Plot and compare the accuracies against the shuffled data accuracies
    dd = shufacc(i, :);

    [zz, mudd, sigdd] = zscore(dd);
    pvalsNOFS(i) = normcdf(nofsacc(i), mudd, sigdd, 'upper');
    pvalsTOPF(i) = normcdf(topfacc(i), mudd, sigdd, 'upper');
    pvalsFS(i) = normcdf(fsacc(i), mudd, sigdd, 'upper');
    
    plotAccuracyDist(savedir, drug, shufacc(i, :), nofsacc(i), topfacc(i), fsacc(i), diffacc(i));
    
    % Create Heatmap plots to visualize the top predictors
    plotFeatsHeatmap(savedir, drug, gdscnet, thaatab, allfeatsfs);
    
end

idxtakegood = pvalsNOFS < 0.05 ; % NOFS is significantly better than random (good data)
idxtakebad = pvalsNOFS >= 0.05; % NOFS is not significantly better than random (good data)

fsBetter = fsacc> topfacc;
diffBetter = diffacc > topfacc;

% Bar plots of accuracies, split into 3 cases/categories: 
% (1) bad signal in data (NOFS not significant), but network selection can still perform well
% (2) good signal in data (NOFS significant) but network selection does not
% performs better than random
% (3) Good signal and network selection performs better than random

% Category 1)
idxtake = idxtakebad & pvalsFS < 0.05;
druglabel = thaatab.DRUG(idxtake);
bardata = [topfacc(idxtake) fsacc(idxtake) diffacc(idxtake)];
plotname = 'cat1_TOPFvsFSvsDIFF';
legendlabel = {'Random forest with top 30 statistical features', 'Random Forest with network-selected features', 'Network selection + diffusion'};
plottitle = 'Comparison of the accuracies between RF with top 30 features, network-selected features (Category 1) and diffused features';
plotAccBar(savedir, plotname, bardata, druglabel, legendlabel, plottitle);

fprintf('Case 1, fs better = %d/%d\n', sum(idxtake & fsBetter), sum(idxtake));
fprintf('Case 1, diff better = %d/%d\n', sum(idxtake & diffBetter), sum(idxtake));

% Category 2)
idxtake = idxtakegood &  pvalsFS >= 0.05;  
druglabel = thaatab.DRUG(idxtake);

bardata = [topfacc(idxtake) fsacc(idxtake) diffacc(idxtake)];
plotname = 'cat2_TOPFvsFSvsDIFF';
legendlabel = {'Random forest with top 30 statistical features', 'Random Forest with network-selected features', 'Network selection + diffusion'};
plottitle = 'Comparison of the accuracies between RF with top 30 features, network-selected features (Category 1) and diffused features';
plotAccBar(savedir, plotname, bardata, druglabel, legendlabel, plottitle);

fprintf('Case 2, fs better = %d/%d\n', sum(idxtake & fsBetter), sum(idxtake));
fprintf('Case 2, diff better = %d/%d\n', sum(idxtake & diffBetter), sum(idxtake));

% Category 3)
idxtake = idxtakegood & pvalsFS < 0.05;
druglabel = thaatab.DRUG(idxtake);
% a) Network feature selection vs statistical baseline (top 30 features) 

bardata = [topfacc(idxtake) fsacc(idxtake) diffacc(idxtake)];
plotname = 'cat3_TOPFvsFSvsDIFF';
legendlabel = {'Random forest with top 30 statistical features', 'Random Forest with network-selected features', 'Network selection + diffusion'};
plottitle = 'Comparison of the accuracies between RF with top 30 features, network-selected features (Category 1) and diffused features';
plotAccBar(savedir, plotname, bardata, druglabel, legendlabel, plottitle, 1);

fprintf('Case 3, fs better = %d/%d\n', sum(idxtake & fsBetter), sum(idxtake));
fprintf('Case 3, diff better = %d/%d\n', sum(idxtake & diffBetter), sum(idxtake));