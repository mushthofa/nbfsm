%% Script for performing cross-validation on ALL dataset

rng('default');

% Put all the result in the following dir
% Make sure this folder GETS CREATED before running the script!
% All results contained already there will be overwritten!!

savedir = 'cvres1/'; % directory to save the result


% Parameters
nfold = 5;
thpvtop = 0.05;

topn = 10; % Top n genes to select using the network scores

% Network kernel parameter
kernelmethod = 'LEXP';
kernelparam = 1e-2; 

% Network diffusion parameters

diffopt.maxiter = 50; % Maximum number of diffusion iterations
diffopt.eps = 1e-8; % Diffusion convergence tolerance
alpha = 0.5:0.05:0.95; % Range of alpha values to try

% Random forest/tree settings

minleafsize = 10; % Minimum leaf size for trees 
maxnumsplit = 5; % Maximum split on each tree
tt = templateTree('MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
mdlopt.tt = tt;
mdlopt.method = 'AdaBoostM1';
mdlopt.numtrees = 50;
mdlopt.maxnumsplit = maxnumsplit;
mdlopt.minleafsize = minleafsize;

% Top features random forest, this should be 3 * number of selected nodes
% from the network feature selection, since there are 3 types of features:
% MUT, CNV, NET. So, if topn=10, we potentially take at most 30 features

ntopfeat = 3*topn;

% How much percentage of total predictive power should we get from the top features 
varimpprop = 0.95;


% Prepare for data load
% load gene entrez id to name mapping

disp('Loading Entrez - Gene name mapping...');
entr2name = readtable('data/entrez2name.csv');
entrez2name = containers.Map('keytype', 'double', 'valuetype', 'char');
name2entrez = containers.Map('keytype', 'char', 'valuetype', 'double');
for i=1:height(entr2name)
    entrez2name(entr2name.entrezid(i)) = entr2name.genename{i};
    name2entrez(entr2name.genename{i}) = entr2name.entrezid(i);
end

% Some extra aliases for gene names are needed becauses of inconsistencies
% in gene naming between datasets

name2entrez('ABL') = 25;
name2entrez('MLL') = 4297;
name2entrez('EWRS1') = 2130;
name2entrez('FTSJD1') = 55783;
name2entrez('MLL2') = 8085;
name2entrez('MLL3') = 58508;

entrezmap.name2entrez = name2entrez;
entrezmap.entrez2name = entrez2name;


% load GDSC (Iorio et al, 2016)
gdsc = loadGDSC();

% Load gene network information
disp('Loading network...');
net = loadNet('net/kegg-acsn-regnetwork.csv');

% Restrict data to nodes in the network
gdscnet = restrictNet(gdsc, net);

disp('Computing marginals...');
gdscnet = marginalExp(gdscnet, 0.9, 3, 1);

% for top features random forest
pnames_topf = strcat(gdsc.mutGenes, '-MUT');
pnames_topf = [pnames_topf; strcat(gdsc.cnvGenes, '-CNV')];
pmat_topf = [gdsc.mutMat; gdsc.cnvMat];

% Data on selected drugs: name, threshold or negative/positive, number of
% cell lines
thaatab = readtable('thaa/quantilethaa_0.30_0.70.csv', 'Delimiter', ',');

save(sprintf('%sdrugdata.mat', savedir), 'gdscnet','net', 'entrez2name', 'thaatab', '-v7.3');

% Set up prediction result holders
% 4 tables for each of NOFS, TOPF, FS and DIFF
clear predictionTable allfeatsnofs allfeatsfs allfeatsdiff;
predictionTable.NOFS = nan(height(thaatab), length(gdscnet.cellNames));
predictionTable.TOPF = nan(height(thaatab), length(gdscnet.cellNames));
predictionTable.FS = nan(height(thaatab), length(gdscnet.cellNames));
predictionTable.DIFF = nan(height(thaatab), length(gdscnet.cellNames));

% Feature tables
allfeatsfs.DRUG = [];
allfeatsfs.FEATS = [];
allfeatsfs.RANK = [];
allfeatsfs.IMP = [];
allfeatsnofs.DRUG = [];
allfeatsnofs.FEATS = [];
allfeatsnofs.RANK = [];
allfeatsnofs.IMP = [];
allfeatsdiff.DRUG = [];
allfeatsdiff.FEATS = [];
allfeatsdiff.RANK = [];
allfeatsdiff.IMP = [];

% For every selected drugs
disp('Doing cross validation per drug ...');

for i=1:height(thaatab)
    
    % Fix RNG on every drug CV run for reproducibility
    rng('default');
    
    % Get positive and negative samples
    drug = thaatab.DRUG{i};
    thaa = [thaatab.THLO(i); thaatab.THHI(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    idxrest = ~idxneg & ~idxpos;
    
    %  Do 5-fold cross validation partition
    idxAll = [find(idxneg)'; find(idxpos)'];
    classAll = [zeros(sum(idxneg), 1); ones(sum(idxpos), 1)];
    cvpat = cvpartition(classAll, 'KFold', nfold);
    
    % First, do a cross-fold validation on the extreme data
    % Compare accuracy, positive recall and negative recall
  
    for j=1:nfold
        
        % Get current training data
        trainidxpos = idxAll(cvpat.training(j) & classAll==1);
        trainidxneg = idxAll(cvpat.training(j) & classAll==0);
        
        % Training target
        classTarget = zeros(length(trainidxneg), 1);
        classTarget = [classTarget; ones(length(trainidxpos), 1)];
        
        % Get current test data
        testidxpos = idxAll(cvpat.test(j) & classAll==1);
        testidxneg = idxAll(cvpat.test(j) & classAll==0);
        
        % Train/test without feature selection
        [fmdl1, ~]  = trainModel(gdscnet, net, trainidxneg, trainidxpos, [], varimpprop, 'NOFS', [], mdlopt);
        
        % Get test features for random forest model
        testFeat = predictorMatrix(gdscnet, fmdl1.PredictorNames);
        testFeat = testFeat(:, [testidxneg; testidxpos]);
        
        % Save predictions
        predictClass = predict(fmdl1, testFeat');
        predictionTable.NOFS(i, [testidxneg; testidxpos]) = predictClass;
        
         % Train/test restricted to top features only 
        predImp = predictorImportance(fmdl1);
        featMat_topf = pmat_topf(:,[trainidxneg; trainidxpos])';
        [~, idximp] = sort(predImp, 'descend');
        featMatImp = featMat_topf(:, idximp(1:ntopfeat));
        pnamesImp = pnames_topf(idximp(1:ntopfeat));
        fmdltop = fitcensemble(featMatImp, classTarget, 'Method', mdlopt.method,  'NumLearningCycles', mdlopt.numtrees, 'Learners',  mdlopt.tt, 'PredictorNames', pnamesImp, 'CategoricalPredictors', 'all');
        
        % Get test features and do prediction
        testFeat = pmat_topf(:, [testidxneg; testidxpos])';
        testFeatImp = testFeat(:, idximp(1:ntopfeat));
        
         % Save predictions
        predictClass = predict(fmdltop, testFeatImp);
        predictionTable.TOPF(i, [testidxneg; testidxpos]) = predictClass;
        
        % Compute network scores from training samples
        drugnet{j} = drugGeneNetAA(gdscnet, idxd, net, trainidxpos, trainidxneg, kernelmethod, kernelparam);
        
        % Network-based feature selection
        selectedGenes = topGenes(drugnet{j}, gdscnet, net, topn);
            
         % Train/test WITH feature selection
        [fmdl2, ~]  = trainModel(gdscnet, net, trainidxneg, trainidxpos, selectedGenes, varimpprop, 'FS', [], mdlopt);
        
        % Get test features
        testFeat = predictorMatrix(gdscnet, fmdl2.PredictorNames);
        testFeat = testFeat(:, [testidxneg; testidxpos]);
        
        % Save predictions
        predictClass = predict(fmdl2, testFeat');
        predictionTable.FS(i, [testidxneg; testidxpos]) = predictClass;
    end
    
    max_acc = -Inf;
    for cur_alpha = alpha
        diffopt.alpha = cur_alpha;
        diffpred = nan(1, length(gdscnet.cellNames));
        referenceClass = nan(1, length(gdscnet.cellNames));
        referenceClass(idxneg) = 0;
        referenceClass(idxpos) = 1;
        referenceClass(isnan(referenceClass)) = [];
        
        for j=1:nfold
             % Get current training data
            trainidxpos = idxAll(cvpat.training(j) & classAll==1);
            trainidxneg = idxAll(cvpat.training(j) & classAll==0);
            
            % training target
            classTarget = zeros(length(trainidxneg), 1);
            classTarget = [classTarget; ones(length(trainidxpos), 1)];

            % Get current test data
            testidxpos = idxAll(cvpat.test(j) & classAll==1);
            testidxneg = idxAll(cvpat.test(j) & classAll==0);
            testClass = zeros(length(testidxneg), 1);
            testClass = [testClass; ones(length(testidxpos), 1)];

            selectedGenes = topGenes(drugnet{j}, gdscnet, net, topn);
            
            % Train WITH feature selection
            [fmdl3, ~]  = trainModel(gdscnet, net, trainidxneg, trainidxpos, selectedGenes, varimpprop, 'FS', [], mdlopt);

            % Get test features, use diffusion on the test feature
            testFeat = predictorMatrix(gdscnet, fmdl3.PredictorNames, net, diffopt);
            testFeat = testFeat(:, [testidxneg; testidxpos]);
            predictClass = predict(fmdl3, testFeat');
            diffpred([testidxneg; testidxpos]) = predictClass;
            
        end
        % Optimize for accuracy on the random forests performance on
        % diffused data
        
        curpred = diffpred;
        curpred(isnan(diffpred)) = [];
        cur_acc = acc(curpred, referenceClass);
        
        if(cur_acc > max_acc)
            max_acc = cur_acc;
            best_alpha = cur_alpha;
            predictionTable.DIFF(i, : ) = diffpred;
        end
    end
    
    % Compute accuracies
    
    referenceClass = nan(1, length(gdscnet.cellNames));
    referenceClass(idxneg) = 0;
    referenceClass(idxpos) = 1;
    referenceClass(isnan(referenceClass)) = [];
    
    nofsPred = predictionTable.NOFS(i, :);
    fsPred = predictionTable.FS(i, :);
    topfPred = predictionTable.TOPF(i, :);
    diffPred = predictionTable.DIFF(i, :);
    
    nofsPred(isnan(nofsPred)) = [];
    fsPred(isnan(fsPred)) = [];
    topfPred(isnan(topfPred)) = [];
    diffPred(isnan(diffPred)) = [];
    
    nofs_acc = acc(nofsPred, referenceClass);
    topf_acc = acc(topfPred, referenceClass);
    fs_acc = acc(fsPred, referenceClass);
    diff_acc = acc(diffPred, referenceClass);
    
    fprintf('%d, %s, %d, %d, %.2f, %.2f, %.2f, %.2f\n', i, drug, sum(idxneg), sum(idxpos), nofs_acc, topf_acc, fs_acc, diff_acc);

%   Now, use all the extreme samples as training, and use them to derive
%   features

    [~, ip, ipvals]  = trainModel(gdscnet, net, find(idxneg)', find(idxpos)', [], varimpprop, 'NOFS', [], mdlopt);
    if(~isempty(ip))
        for j=1:length(ip)
            allfeatsnofs.DRUG = [allfeatsnofs.DRUG; {drug}];
            allfeatsnofs.FEATS = [allfeatsnofs.FEATS; {ip{j}}];
            allfeatsnofs.RANK = [allfeatsnofs.RANK; j];
            allfeatsnofs.IMP = [allfeatsnofs.IMP; ipvals(j)];
        end
    end

%   Train/test WITH feature selection
%   Feature selection 
    gdsc_drugnet{i} = drugGeneNetAA(gdscnet, idxd, net, idxpos, idxneg, kernelmethod, kernelparam);
    netGenes{i} = topGenes(gdsc_drugnet{i}, gdscnet, net, topn);
    
    [~, ip, ipvals]  = trainModel(gdscnet, net, find(idxneg)', find(idxpos)', netGenes{i}, varimpprop, 'FS', [], mdlopt);
    
    if(~isempty(ip))
        for j=1:length(ip)
            allfeatsfs.DRUG = [allfeatsfs.DRUG; {drug}];
            allfeatsfs.FEATS = [allfeatsfs.FEATS; {ip{j}}];
            allfeatsfs.RANK = [allfeatsfs.RANK; j];
            allfeatsfs.IMP = [allfeatsfs.IMP; ipvals(j)];
        end
    end
   
     % Train/test WITH feature selection AND diffusion
    [~, ip, ipvals]  = trainModel(gdscnet, net, find(idxneg)', find(idxpos)', netGenes{i}, varimpprop, 'DIFF', diffopt, mdlopt);

    if(~isempty(ip))
        for j=1:length(ip)
            allfeatsdiff.DRUG = [allfeatsdiff.DRUG; {drug}];
            allfeatsdiff.FEATS = [allfeatsdiff.FEATS; {ip{j}}];
            allfeatsdiff.RANK = [allfeatsdiff.RANK; j];
            allfeatsdiff.IMP = [allfeatsdiff.IMP; ipvals(j)];
        end
    end
    
end

% Finished...
% Write results to output

% Predictor tables
allfeatsfs = struct2table(allfeatsfs);
allfeatsnofs = struct2table(allfeatsnofs);
allfeatsdiff = struct2table(allfeatsdiff);

% Reference class
referenceTable = nan(height(thaatab), length(gdscnet.cellNames));
for i=1:height(thaatab)
    drug = thaatab.DRUG{i};
    thaa = [thaatab.THLO(i); thaatab.THHI(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    referenceTable(i, idxneg) = 0;
    referenceTable(i, idxpos) = 1;
end
writePrediction(sprintf('%sreference_class.csv', savedir), gdscnet, thaatab, referenceTable);

save(sprintf('%scvresult.mat', savedir), 'predictionTable', 'referenceTable', 'gdsc_drugnet', ...
           'netGenes', 'allfeatsfs', 'allfeatsnofs', 'allfeatsdiff', '-v7.3');

writetable(allfeatsfs, sprintf('%sallfeatsfs_topn_%d.csv', savedir,  topn));
writetable(allfeatsnofs, sprintf('%sallfeatsnofs_topn_%d.csv', savedir,  topn));
writetable(allfeatsdiff, sprintf('%sallfeatsdiff_topn_%d.csv', savedir,  topn));


% Prediction output and reference
writePrediction(sprintf('%sprediction_file_nofs.csv', savedir), gdscnet, thaatab, predictionTable.NOFS);
writePrediction(sprintf('%sprediction_file_topf.csv', savedir), gdscnet, thaatab, predictionTable.TOPF);
writePrediction(sprintf('%sprediction_file_fs.csv', savedir), gdscnet, thaatab, predictionTable.FS);
writePrediction(sprintf('%sprediction_file_diff.csv', savedir), gdscnet, thaatab, predictionTable.DIFF);





