function writePrediction(outfile, ds, thtab, predtable)
fpred = fopen(outfile, 'w+');
fprintf(fpred, 'DRUG');
for i=1:length(ds.cellNames)
    fprintf(fpred, ', %s', ds.cellNames{i});
end
fprintf(fpred, '\n');

for i=1:height(thtab)
    drug = thtab.DRUG{i};
    idxd = strcmpi(ds.allDrugs,drug);
    fprintf(fpred, '%s', ds.allDrugs{idxd});
    fprintf(fpred, num2str(predtable(i, :), ',%d'));
    fprintf(fpred,'\n');
end
fclose(fpred);

end