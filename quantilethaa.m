quantileth = [0.3 0.7];
filesaveth = sprintf('thaa/quantilethaa_%.2f_%.2f.csv', quantileth(1), quantileth(2));
ff = fopen(filesaveth, 'w');
fprintf(ff, 'DRUG, THLO, THHI, NUMALL, NUMNEG, NUMPOS\n');
passdrugs = {};
passidx = 0;
for i=1:length(gdsc.allDrugs)
    drug = gdsc.allDrugs{i};
    aa = gdsc.AAMat(i, :)';
    aa(isnan(aa)) = [];
    n_all = length(aa);
    thaa = quantile(aa, quantileth);
    thlo = min(thaa(1), 0.1);
    n_neg = sum(aa<=thlo);
    thhi = max(thaa(2), 0.1);
    if(thhi<=thlo)
        error('Overlaps!');
    end
    n_pos = sum(aa>=thhi);
    
    % Special treatment for Mitomycin C
    pass = 0;
    if(strcmpi(drug, 'mitomycin c'))
        thaa = quantile(aa, quantileth);
        thlo = min(thaa(1), 0.2);
        n_neg = sum(aa<=thlo);
        thhi = max(thaa(2), 0.1);
        n_pos = sum(aa>=thhi);
        pass = 1;
    end
    
    if((n_neg>=30 && n_pos>=30) || pass==1)
        fprintf(ff, '%s, %.4f, %.4f, %d, %d, %d\n', drug, thlo, thhi, n_all, n_neg, n_pos);
        idxlo = aa<=thlo;
        idxhi = aa>=thhi;
        edges=0:0.05:1;
        hold on;
        blo = histc(aa(idxlo), edges);
        bhi = histc(aa(idxhi), edges);
        bmi = histc(aa(~idxlo & ~idxhi), edges);
        bar(edges, blo, 'r');
        bar(edges, bmi, 'b');
        bar(edges, bhi, 'g');
        cf=gcf;
        ax=cf.CurrentAxes;
        ax.XLim = [-0.1 1.1];
        cf.Position = [0 0 800 600];
        title(sprintf('Activity Area threshold for drug %s (total samples = %d)', drug, n_all));
        legend(sprintf('Negatives = %d', n_neg),'Middle (unused)',sprintf('Positives = %d', n_pos));
        print(sprintf('thaa/%s_%.2f_%.2f.png', matlab.lang.makeValidName(drug), quantileth(1), quantileth(2)), '-dpng');
        close all;
    end
end
fclose(ff);