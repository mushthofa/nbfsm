function drugnet = drugGeneNetAA(ds, idxdrug, net, idxpos, idxneg, kernelmethod, kernelparam)
    % Perform network integration and kernel calculation for responsive and
    % resistant sampels of a certain drug


    posprof = ds.AAMat(idxdrug, idxpos);
    negprof = ds.AAMat(idxdrug, idxneg);

    posprof = (posprof-min(posprof))./(max(posprof)-min(posprof));
    negprof = (negprof-min(negprof))./(max(negprof)-min(negprof));
    negprof = 1 - negprof;
   
    drugnet.posprof = posprof;
    drugnet.negprof = negprof;

    % Build positive global network & Computing positive kernel

    globpos = buildNetMat(drugnet.posprof', ds.mutMat(:, idxpos), ds.mutGenes, ds.cnvMat(:, idxpos)~=0, ds.cnvGenes, ds.dgexMat(:, idxpos)~=0, ds.dgexGenes, net);
    kernelmat = kernel(globpos.mat, kernelmethod, kernelparam);
    kernelmat = kernelNorm(kernelmat);
    simpos = mean(kernelmat(globpos.startSample:globpos.endSample, 1:end));
    globpos = rmfield(globpos, 'mat'); % Save memory by removing unnecessary field

    % Build negative global network & Computing negative kernel
    globneg = buildNetMat(drugnet.negprof', ds.mutMat(:, idxneg), ds.mutGenes, ds.cnvMat(:, idxneg)~=0, ds.cnvGenes,  ds.dgexMat(:, idxneg)~=0, ds.dgexGenes, net);
    kernelmat = kernel(globneg.mat, kernelmethod, kernelparam);
    kernelmat = kernelNorm(kernelmat);
    simneg = mean(kernelmat(globneg.startSample:globneg.endSample, 1:end));
    globneg = rmfield(globneg, 'mat');

    drugnet.simpos = simpos;
    drugnet.simneg = simneg;
    drugnet.globpos = globpos;
    drugnet.globneg = globneg;
end







