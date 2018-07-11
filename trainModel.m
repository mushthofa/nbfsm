function [mdl, ip, ipvals]  = trainModel(ds, net, idxneg, idxpos, selectedGenes, varimpprop, method, diffopt, mdlopt)

%% Train a random forest model for predicting drug response
%   either without network-based feature selection, with feature selection
%   or with feature selection + diffusion
%
%   [mdl, predMat]  = trainModel(ds, net, idxneg, idxpos, selectedGenes, topn, method, diffopt, mdlopt)

%   
%   method = 'NOFS' | 'FS' | 'DIFF'
%%

% Diffusion options
doDiffuse = 0;
if(~isempty(diffopt))
    doDiffuse = 1;
    if(~isempty(diffopt))
        alpha = diffopt.alpha;
        maxiter = diffopt.maxiter;
        eps = diffopt.eps;
    else
        alpha = 0.7;
        maxiter = 50;
        eps = 1e-5;
    end
    
end


if(~isempty(selectedGenes))
    mutgenes = selectedGenes.MUT;
    cnvgenes = selectedGenes.CNV;
    netgenes = selectedGenes.NET;
 

    idxnode = ismember(netgenes, ds.mutGenes); 
    mutgenes = union(mutgenes, netgenes(idxnode), 'stable');
    idxnode = ismember(netgenes, ds.cnvGenes);
    cnvgenes = union(cnvgenes, netgenes(idxnode), 'stable');
 

    % Annoying problem with the cell shape, if the gene list is of length 1,
    % then union will reshape it into a row cell array.
    % we make sure it is a column

    mutgenes = reshape(mutgenes, length(mutgenes), 1);
    cnvgenes = reshape(cnvgenes, length(cnvgenes), 1);
    
   
%     % set of genes for diffusion
    diffgenes = union(mutgenes, cnvgenes);


    idxmut = ismember(ds.mutGenes, mutgenes);
    idxcnv = ismember(ds.cnvGenes, cnvgenes);
end


% construct the predictor matrix from the given training data
% either using feature selection or not

featMat = [];
if(strcmpi(method, 'NOFS'))
    featMat = [featMat; ds.mutMat(:, [idxneg; idxpos])];
    featMat = [featMat; ds.cnvMat(:, [idxneg; idxpos])];
else
    featMat = [featMat; ds.mutMat(idxmut, [idxneg; idxpos])];
    featMat = [featMat; ds.cnvMat(idxcnv, [idxneg; idxpos])];
end


if(strcmpi(method,'DIFF'))
    % Diffusion method is selected, diffuse training data first
    
    % get subnetwork of selected genes
    idxdiffnet = cell2mat(values(net.name2node, diffgenes));
    adjmat = net.mat(idxdiffnet, idxdiffnet);
    adjmat = single((adjmat + adjmat')>0);
    
    % Diffuse mutations
   
    % get mutation data
    idxmut = cellfun(@(x) find(ismember(ds.mutGenes, x)), mutgenes);
    mutmat = ds.mutMat(idxmut, [idxneg; idxpos]);
    
    % create diffusion matrix holding mutation data
    diffmutMat = zeros(length(diffgenes), size(mutmat, 2));
    idxmutdiff = cellfun(@(x) find(ismember(diffgenes, x)), mutgenes);
    diffmutMat(idxmutdiff, :) = mutmat;
    
    diffmutMat = netSmooth(diffmutMat', adjmat, alpha, eps, maxiter)';
    diffmutMat = double(diffmutMat >= (1-alpha));
    diffmutMat = diffmutMat(idxmutdiff, :);
    
    % Diffuse Amplifications
    
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, [idxneg; idxpos]) > 0); % Find only 1's
    
       
    diffAmfeatMat = zeros(length(diffgenes), size(cnvmat, 2));
    idxcnvdiff = cellfun(@(x) find(ismember(diffgenes, x)), cnvgenes);
    diffAmfeatMat(idxcnvdiff, :) = cnvmat;
    
    diffAmfeatMat = netSmooth(diffAmfeatMat', adjmat, alpha, eps, maxiter)';
    diffAmfeatMat = diffAmfeatMat(idxcnvdiff, :);
    
    % Diffuse Deletions
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, [idxneg; idxpos]) < 0); % Find only -1's
    
       
    diffDelfeatMat = zeros(length(diffgenes), size(cnvmat, 2));
    idxcnvdiff = cellfun(@(x) find(ismember(diffgenes, x)), cnvgenes);
    diffDelfeatMat(idxcnvdiff, :) = cnvmat;
    
    diffDelfeatMat = netSmooth(diffDelfeatMat', adjmat, alpha, eps, maxiter)';
    diffDelfeatMat = diffDelfeatMat(idxcnvdiff, :);
    
    % Combine diffused amplifications and deletion data
    
    diffCNVMat = double(abs(diffAmfeatMat - diffDelfeatMat) >= (1-alpha)) .* sign((diffAmfeatMat - diffDelfeatMat));
    
    featMat = [];
    featMat = [featMat; diffmutMat];
    featMat = [featMat; diffCNVMat];
    
end

% Train the model

% set up predictor names
if(strcmpi(method, 'DIFF') || strcmpi(method, 'FS'))
    pnames = strcat(ds.mutGenes(idxmut), '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes(idxcnv), '-CNV')];
else
    pnames = strcat(ds.mutGenes, '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes, '-CNV')];
end

    % Training target
    classTarget = zeros(length(idxneg), 1);
    classTarget = [classTarget; ones(length(idxpos), 1)];

    % Do training
    featMat = featMat';
    mdl = fitcensemble(featMat, classTarget, 'Method', mdlopt.method, 'NumLearningCycles', mdlopt.numtrees, 'Learners', mdlopt.tt, 'PredictorNames', pnames, 'CategoricalPredictors', 'all');
        
    % Get train/test features for single tree model
    % Get important predictors from RF
    
    ii = predictorImportance(mdl);
    ip = mdl.PredictorNames; 
    % Take only predictors with non-zero importance
    ip = ip(ii>0);
    ii = ii(ii>0);
    
    % Get most important predictors (totalling at most varimpprop % of predictor
    % importance)
    [~, idx] = sort(ii, 'descend');
    topn = find((cumsum(ii)./sum(ii))>=varimpprop); 
    ip = ip(idx);
    ip = ip(1:topn); % Take top-n predictors
    ipvals = ii./sum(ii)*100;
    ipvals = ipvals(idx);
    ipvals = ipvals(1:topn);
    
end

