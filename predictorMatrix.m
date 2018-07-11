function [pmat, mutgenes, cnvgenes, gexgenes, tissfeats] = predictorMatrix(ds, feats, varargin)
%% Extracting a matrix of feature values from data set based on a set of feature names
% [pmat, mutgenes, cnvgenes, gexgenes, tissfeats] = predictorMatrix(ds, feats, varargin)
%


% If network + diffusion information are supplied, then also do diffusion
doDiffuse = 0;
if(length(varargin)>=2 && ~isempty(varargin{2}))
    net = varargin{1};
    diffopt = varargin{2};
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

doIgnore = 0;
if(length(varargin)>=3)
    doIgnore = varargin{3};
end


pmat = [];
cnvgenes = cell(0);
mutgenes = cell(0);
%metgenes = cell(0);
gexgenes = cell(0);
tissfeats = cell(0);




idxmutfeat = [];
%idxmetfeat = [];
idxcnvfeat = [];


for i=1:length(feats)
    
    ff = strsplit(feats{i}, '-');
    genename = strjoin(ff(1:end-1), '-');
    ptype = ff{end};
    
    switch(ptype)
        case 'MUT'
            idxg = strcmpi(ds.mutGenes, genename);
            if(sum(idxg)~=1)
                if(doIgnore==0)
                    error('Feature %s not found\n', feats{i});
                end
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.mutMat(idxg, :)];
                mutgenes = [mutgenes; genename];
                idxmutfeat = [idxmutfeat; i];
            end
        case 'CNV'
            idxg = strcmpi(ds.cnvGenes, genename);
            if(sum(idxg)~=1)
                if(doIgnore==0)
                    error('Feature %s not found\n', feats{i});
                end
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.cnvMat(idxg, :)];
                cnvgenes = [cnvgenes; genename];
                idxcnvfeat = [idxcnvfeat; i];
            end
        case 'GEX'
            idxg = strcmpi(ds.gexGenes, genename);
            if(sum(idxg)~=1)
                if(doIgnore==0)
                    error('Feature %s not found\n', feats{i});
                end
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.gexMat(idxg, :)];
                gexgenes = [gexgenes; genename];
                idxgexfeat = [idxgexfeat; i];
            end
            
%         case 'MET'
%             idxg = strcmpi(ds.metGenes, genename);
%             if(sum(idxg)~=1)
%                 if(doIgnore==0)
%                     error('Feature %s not found\n', feats{i});
%                 end
%                 pmat = [pmat; zeros(1, length(ds.cellNames))];
%             else
%                 pmat = [pmat; ds.metMat(idxg, :)];
%                 metgenes = [metgenes; genename];
%                 idxmetfeat = [idxmetfeat; i];
%             end
        otherwise
            tissfeat = feats{i};
            idxt = length('tissue_is_') + 1;
            tissname = tissfeat(idxt:end);
            idxtiss = strcmpi(alltissues, tissname);
            if(sum(idxtiss)~=1)
                
                if(doIgnore == 0)
                    error('Unknown predictor type %s', ptype);
                end
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; tissmat(idxtiss, :)];
                tissfeats = [tissfeats; tissname];
                idxtissfeat = [idxtissfeat; i];
            end
            
    end
    
end

% IF need to do diffusion

if(doDiffuse == 1)
    % Do not forget to diffuse test data as well
    diffgenes = union(mutgenes, cnvgenes);
    
    % get subnetwork of selected genes
    idxdiffnet = cell2mat(values(net.name2node, diffgenes));
    adjmat = net.mat(idxdiffnet, idxdiffnet);
    adjmat = single((adjmat + adjmat')>0);
    
    % Diffuse mutations
   
    % get mutation data
    idxmut = cellfun(@(x) find(ismember(ds.mutGenes, x)), mutgenes);
    mutmat = ds.mutMat(idxmut, :);
    
    % create diffusion matrix holding mutation data
    diffmutMat = zeros(length(diffgenes), size(mutmat, 2));
    idxmutdiff = cellfun(@(x) find(ismember(diffgenes, x)), mutgenes);
    diffmutMat(idxmutdiff, :) = mutmat;
    
    diffmutMat = netSmooth(diffmutMat', adjmat, alpha, eps, maxiter)';
    diffmutMat = double(diffmutMat >= (1-alpha));
    diffmutMat = diffmutMat(idxmutdiff, :);
    
    % Diffuse Amplifications
    
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, :) > 0); % Find only 1's
    
       
    diffAmfeatMat = zeros(length(diffgenes), size(cnvmat, 2));
    idxcnvdiff = cellfun(@(x) find(ismember(diffgenes, x)), cnvgenes);
    diffAmfeatMat(idxcnvdiff, :) = cnvmat;
    
    diffAmfeatMat = netSmooth(diffAmfeatMat', adjmat, alpha, eps, maxiter)';
    diffAmfeatMat = diffAmfeatMat(idxcnvdiff, :);
    
    % Diffuse Deletions
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, :) < 0); % Find only -1's
    
       
    diffDelfeatMat = zeros(length(diffgenes), size(cnvmat, 2));
    idxcnvdiff = cellfun(@(x) find(ismember(diffgenes, x)), cnvgenes);
    diffDelfeatMat(idxcnvdiff, :) = cnvmat;
    
    diffDelfeatMat = netSmooth(diffDelfeatMat', adjmat, alpha, eps, maxiter)';
    diffDelfeatMat = diffDelfeatMat(idxcnvdiff, :);
    
    % Combine diffused amplifications and deletion data
    
    diffCNVMat = double(abs(diffAmfeatMat - diffDelfeatMat) >= (1-alpha)) .* sign((diffAmfeatMat - diffDelfeatMat));
    
    % Diffuse methylation
   
%     % get methylation data
%     idxmet = cellfun(@(x) find(ismember(ds.metGenes, x)), metgenes);
%     metmat = ds.metMat(idxmet, :);
%     
%     % create diffusion matrix holding mutation data
%     diffmetMat = zeros(length(diffgenes), size(metmat, 2));
%     idxmetdiff = cellfun(@(x) find(ismember(diffgenes, x)), metgenes);
%     diffmetMat(idxmetdiff, :) = metmat;
%     
%     diffmetMat = netSmooth(diffmetMat', adjmat, alpha, eps, maxiter)';
%     diffmetMat = double(diffmetMat >= (1-alpha));
%     diffmetMat = diffmetMat(idxmetdiff, :);
%     
    
    featMat = zeros(size(pmat));
    featMat(idxmutfeat, :) = diffmutMat;
    featMat(idxcnvfeat, :) = diffCNVMat;
    %featMat(idxmetfeat, :) = diffmetMat;
    
%     % GEX and Tissue features are taken as-is
%     idx = cellfun(@(x) find(strcmpi(ds.gexGenes, x)), gexgenes);
%     featgex = ds.gexMat(idx, :);
%     featMat(idxgexfeat, :) = featgex;
%     
%     idx = cellfun(@(x) find(strcmpi(alltissues, x)), tissfeats);
%     feattiss = tissmat(idx, :);
%     featMat(idxtissfeat, :) = feattiss;
    pmat = featMat;
end


end