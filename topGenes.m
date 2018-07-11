function toptbl = topGenes(dgnet, ds, net, topn)
% Analysing kernel results and obtain a list of top selected features

    startmutp = dgnet.globpos.startMut;
    endmutp = dgnet.globpos.endMut;
    startcnvp = dgnet.globpos.startCNV;
    endcnvp = dgnet.globpos.endCNV;

    
    startmutn = dgnet.globneg.startMut;
    endmutn = dgnet.globneg.endMut;
    startcnvn = dgnet.globneg.startCNV;
    endcnvn = dgnet.globneg.endCNV;
    
%     startmetp = dgnet.globpos.startMET;
%     endmetp = dgnet.globpos.endMET;
%     startmetn = dgnet.globneg.startMET;
%     endmetn = dgnet.globneg.endMET;
    
    
    startgexn = dgnet.globneg.startGEX;
    endgexn = dgnet.globneg.endGEX;
    
    if(length(dgnet.simpos) < dgnet.globpos.size)
        dgnet.simpos = zeros(1, dgnet.globpos.size);
    end
    
    if(length(dgnet.simneg) < dgnet.globneg.size)
        dgnet.simneg = zeros(1, dgnet.globneg.size);
    end
    
    dgnet.simpos(isnan(dgnet.simpos)) = 0;
    dgnet.simneg(isnan(dgnet.simneg)) = 0;
    
    simdrugpos = dgnet.simpos*10e5;
    simdrugneg = dgnet.simneg*10e5;
	
    smut = [simdrugpos(startmutp:endmutp)' simdrugneg(startmutn:endmutn)'];
    scnv = [simdrugpos(startcnvp:endcnvp)' simdrugneg(startcnvn:endcnvn)'];
%     smet = [simdrugpos(startmetp:endmetp)' simdrugneg(startmetn:endmetn)'];
    
    
    smut = quantnorm(smut, 'median');
    scnv = quantnorm(scnv, 'median');
%     smet = quantnorm(smet, 'median');

    
    smut = smut(:, 1) - smut(:, 2);
%     smet = smet(:, 1) - smet(:, 2);
    scnv = scnv(:, 1) - scnv(:, 2);
    

    simnetp = dgnet.simpos(1:net.size)*10e5;
    simnetn = dgnet.simneg(1:net.size)*10e5;
    
    simnet = quantnorm([simnetp' simnetn'], 'median');
    simnet = simnet(:, 1) - simnet(:, 2);

    [~, topmutidx] = sort(abs(smut), 'descend');
    topmutgenes = ds.mutGenes(topmutidx(1:topn));
    
    
    [~, topidx] = sort(abs(scnv), 'descend');
    topcnvgenes = ds.cnvGenes(topidx(1:topn));
    
%     [~, topidx] = sort(abs(smet), 'descend');
%     topmetgenes = ds.metGenes(topidx(1:topn));
    
    
    [~, idx] = sort(abs(simnet), 'descend');
    topnetgenes = net.nodes(idx(1:topn));
    
   
    
    maxl = max([length(topmutgenes); length(topcnvgenes); length(topnetgenes);]);
    
    
    topmutgenes = reshape(topmutgenes, maxl, 1);
    
    topcnvgenes = reshape(topcnvgenes, maxl, 1);
    
    topnetgenes = reshape(topnetgenes, maxl, 1);
    
    toptbl = table(topmutgenes,topcnvgenes, topnetgenes, ... '
        'VariableNames', {'MUT', 'CNV', 'NET'});
    
end