function glob = buildNetMat(prof, mutmat, muteid, cnvmat, cnveid, gexmat, gexeid, net)


% Sizes and indexes
nsample = size(mutmat, 2);


nmut = size(mutmat, 1);
ncnv = size(cnvmat, 1);
%nmet = size(metmat, 1);
ngex = size(gexmat, 1);

glob.size = net.size + nsample + nmut + ncnv + ngex;
glob.mat = single(zeros(glob.size, glob.size));

glob.startSample = net.size + 1;
glob.endSample = glob.startSample + nsample - 1;

glob.startMut = glob.endSample + 1;
glob.endMut = glob.startMut + nmut -1;

glob.startCNV = glob.endMut + 1;
glob.endCNV = glob.startCNV + ncnv - 1;

%glob.startMET = glob.endCNV + 1;
%glob.endMET = glob.startMET + nmet - 1;

glob.startGEX = glob.endCNV + 1;
glob.endGEX = glob.startGEX + ngex - 1;

% Fill up data nodes
% Give weights for connections from samples to data nodes using normalised profile
glob.mat(1:net.size, 1:net.size) = net.mat;
glob.mat(glob.startSample:glob.endSample, glob.startMut:glob.endMut) = repmat(prof, 1, size(mutmat, 1)).*mutmat';
glob.mat(glob.startSample:glob.endSample, glob.startCNV:glob.endCNV) = repmat(prof, 1, size(cnvmat, 1)).*cnvmat';
%glob.mat(glob.startSample:glob.endSample, glob.startMET:glob.endMET) = repmat(prof, 1, size(metmat, 1)).*metmat';
glob.mat(glob.startSample:glob.endSample, glob.startGEX:glob.endGEX) = repmat(prof, 1, size(gexmat, 1)).*gexmat';


% Fill up connection nodes
mutnodes = single(zeros(net.size, nmut));
idxrow = cell2mat(values(net.name2node, muteid));
idxcol = (1:nmut)';
sidx = sub2ind(size(mutnodes), idxrow, idxcol);
mutnodes(sidx) = 1;
glob.mat(1:net.size, glob.startMut:glob.endMut) = mutnodes;

cnvnodes = single(zeros(net.size, ncnv));
idxrow = cell2mat(values(net.name2node, cnveid));
idxcol = (1:ncnv)';
sidx = sub2ind(size(cnvnodes), idxrow, idxcol);
cnvnodes(sidx) = 1;
glob.mat(1:net.size, glob.startCNV:glob.endCNV) = cnvnodes;

% metnodes = single(zeros(net.size, nmet));
% idxrow = cell2mat(values(net.name2node, meteid));
% idxcol = (1:nmet)';
% sidx = sub2ind(size(metnodes), idxrow, idxcol);
% metnodes(sidx) = 1;
% glob.mat(1:net.size, glob.startMET:glob.endMET) = metnodes;

gexnodes = single(zeros(net.size, ngex));
idxrow = cell2mat(values(net.name2node, gexeid));
idxcol = (1:ngex)';
sidx = sub2ind(size(gexnodes), idxrow, idxcol);
gexnodes(sidx) = 1;
glob.mat(1:net.size, glob.startGEX:glob.endGEX) = gexnodes;

% 
glob.mat = single((glob.mat + glob.mat') > 0);
glob.mat(1:glob.size + 1: end) = 1;

