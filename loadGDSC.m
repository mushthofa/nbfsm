function ds = loadGDSC()
% Loading and preprocessing GDSC data

    % Loading from the pan-cancer data set
    disp('Processing CFE table...');
    cfe = readtable('data/PANCAN_simple_MOBEM.rdata.txt', 'Delimiter', '\t', 'ReadVariableNames', 0);
    gdsc_cells = readtable('data/GDSC_Cells_MappedTissue.csv');
    cfe_mat = table2array(cfe(2:end, 2:end));
    cfe_events = cfe.Var1(2:end);
    cfe_cellids = cfe(1, 2:end);
    cfe_cellids = table2array(cfe_cellids);
    
    cfe_cellnameidx = arrayfun(@(x) find(gdsc_cells.CELL_ID==x), cfe_cellids, 'UniformOutput', false);
    cfe_cellnameidx = cell2mat(cfe_cellnameidx);
    cfe_cellnames = gdsc_cells.CELL_NAME(cfe_cellnameidx);
    cfe_celltissues = gdsc_cells.MAPPEDTISS(cfe_cellnameidx);
    cfe_oritissues = gdsc_cells.TISSUE2(cfe_cellnameidx);
    
    mut_genes = {};
    mut_idx = [];

    cnvu_genes = {};
    cnvu_idx = [];

    cnvd_genes = {};
    cnvd_idx = [];

    meth_genes = {};
    meth_idx = [];


    for i=1:length(cfe_events)
        if(strcmp(cfe_events{i}(end-3:end), '_mut'))
            genename = cfe_events{i}(1:end-4);
            
            mut_genes = [mut_genes; genename];
            mut_idx = [mut_idx; i];
            
        elseif(strcmp(cfe_events{i}(1:5), 'gain:'))
            genenamesIdx =  regexp(cfe_events{i}, '\(([^()]*)\)');
            if(isempty(genenamesIdx))
                continue;
            end

            genenames = cfe_events{i}(genenamesIdx + 1:end-1);
            genenames = strtrim(strsplit(genenames, ','));
            
            

            cnvu_genes = [cnvu_genes; genenames'];
            cnvu_idx = [cnvu_idx; repmat(i, length(genenames), 1)];

        elseif(strcmp(cfe_events{i}(1:5), 'loss:'))
            genenamesIdx =  regexp(cfe_events{i}, '\(([^()]*)\)');
            if(isempty(genenamesIdx))
                continue;
            end

            genenames = cfe_events{i}(genenamesIdx + 1:end-1);
            genenames = strtrim(strsplit(genenames, ','));
            if(isempty(genenames))
                continue;
            end
           

            cnvd_genes = [cnvd_genes; genenames'];
            cnvd_idx = [cnvd_idx; repmat(i, length(genenames), 1)];

        elseif(strcmp(cfe_events{i}(1:3), 'chr'))
            genenamesIdx =  regexp(cfe_events{i}, '\(([^()]*)\)');
            if(isempty(genenamesIdx))
                continue;
            end

            genenames = cfe_events{i}(genenamesIdx + 1:end-8);
            genenames = strtrim(strsplit(genenames, ','));
            if(isempty(genenames) || (length(genenames)==1 && isempty(genenames{1})))
                continue;
            end
            meth_genes = [meth_genes; genenames'];
            meth_idx = [meth_idx; repmat(i, length(genenames), 1)];
        end
    end
    
    % Process MUT, CNV and MET matrices
    
    mutmat = zeros(size(mut_genes, 1), length(cfe_cellids));
    for i=1:size(mut_genes, 1)
        mutmat(i, :) = cfe_mat(mut_idx(i), :);
    end

    cnv_genes = union(cnvu_genes, cnvd_genes);
    cnvmat = zeros(size(cnv_genes, 1), length(cfe_cellids));
    for i=1:size(cnv_genes, 1)
        idxu = strcmp(cnvu_genes, cnv_genes{i});
        idxd = strcmp(cnvd_genes, cnv_genes{i});

        % Average out repeated gene data
        if(sum(idxu)>0)
            arru = double(mean(cfe_mat(cnvu_idx(idxu), :), 1)>0);
        else
            arru = zeros(1, length(cfe_cellids));
        end

        if(sum(idxd)>0)
            arrd = double(mean(cfe_mat(cnvd_idx(idxd), :), 1)>0);
        else
            arrd = zeros(1, length(cfe_cellids));
        end

        cnvmat(i, :) = arru-arrd;
    end

    methmat = zeros(length(meth_genes), length(cfe_cellids));
    toDel = [];
    for i=1:length(meth_genes)
        idxu = strcmp(meth_genes,  meth_genes{i});
        % Max-out repeated gene data
        if(sum(idxu)>1)
            arru = max(cfe_mat(meth_idx(idxu), :));
        else
            arru = cfe_mat(meth_idx(idxu), :);
        end

        fi = find(idxu);
        methmat(fi(1), :) = arru;
        toDel = [toDel; fi(2:end)];
    end
    
    methmat(toDel, :) = [];
    meth_genes(toDel) = [];

    % Consider only samples with minimum number of tissues
    minsample = 10;
    tt=tabulate(cfe_celltissues);
    tisslow = tt(cell2mat(tt(:, 2))<minsample, 1);
    
    idxlow = ismember(cfe_celltissues, tisslow);
    
    % Delete tissues with low number of samples
    cfe_cellnames(idxlow) = [];
    cfe_celltissues(idxlow) = [];
    cfe_oritissues(idxlow) = [];
    
    mutmat(:, idxlow) = [];
    cnvmat(:, idxlow) = [];
    methmat(:, idxlow) = [];
    
    
    % Read GEX data
    disp('Reading gene expression data...');
    genexpTab = readtable('data/Cell_line_RMA_proc_basalExp.txt', 'Delimiter', '\t');
    genes_exp = genexpTab.GENE_SYMBOLS;
    
    expmat = table2array(genexpTab(:, 3:end));
    cellids = genexpTab.Properties.VariableNames(3:end);
    expsids = cellfun(@(x) str2double(x(6:end)), cellids, 'UniformOutput', 0);
    
    % Delete the sample ID's that cannot be converted to cell names
    idxnanemp = cellfun(@(x) isempty(x)|isnan(x), expsids);
    expmat(:, idxnanemp) = [];
    expsids(idxnanemp) = [];
    expsids = cell2mat(expsids);
    
    % Map Cell IDs to standardized cell names
    idxnav = ~ismember(expsids, gdsc_cells.CELL_ID);
    expsids(idxnav) = [];
    expmat(:, idxnav) = [];
    idxcell = arrayfun(@(x) find(gdsc_cells.CELL_ID==x), expsids);
    expcellnames = gdsc_cells.CELL_NAME(idxcell);
       
    
    icl = intersect(cfe_cellnames, expcellnames);
    idxc1 = ismember(cfe_cellnames, icl);
    idxc2 = ismember(expcellnames, icl);
    
    % Consider only cell lines in the intersection
    mutmat(:, ~idxc1) = [];
    cnvmat(:, ~idxc1) = [];
    methmat(:, ~idxc1) = [];
    cfe_cellnames(~idxc1) = [];
    cfe_celltissues(~idxc1) = [];
    cfe_oritissues(~idxc1) = [];
    
    expmat(:, ~idxc2) = [];
    expcellnames(~idxc2) = [];
    
    % Rearrange sample labels according to the intersection
    idxc1 = cellfun(@(x) find(strcmp(x, cfe_cellnames)), icl);
    idxc2 = cellfun(@(x) find(strcmp(x, expcellnames)), icl);
    
    mutmat = mutmat(:, idxc1);
    cnvmat = cnvmat(:, idxc1);
    methmat = methmat(:, idxc1);
    cfe_cellnames = cfe_cellnames(idxc1);
    cfe_celltissues = cfe_celltissues(idxc1);
    cfe_oritissues = cfe_oritissues(idxc1);
    
    expmat = expmat(:, idxc2);
    
    % Save the data
    
    ds.mutMat = mutmat;
    ds.mutGenes = mut_genes;
    ds.cnvMat = cnvmat;
    ds.cnvGenes = cnv_genes;
    ds.metMat = methmat;
    ds.metGenes = meth_genes;
    ds.gexMat = expmat;
    ds.gexGenes = genes_exp;
    ds.cellNames = cfe_cellnames;
    ds.cellTissues = cfe_celltissues;
    ds.cellOriTissues = cfe_oritissues;
    
    idxemptyname = strcmpi(ds.gexGenes, '');
    ds.gexGenes(idxemptyname) = [];
    ds.gexMat(idxemptyname, :) = [];
    
    % Read drug response data
    disp('Reading drug response data...');
    drugIDmap = containers.Map('KeyType','int32', 'ValueType', 'char');
    gdscdrugs = readtable('data/GDSC_DRUGS.csv');
    for i=1:height(gdscdrugs)
        drugIDmap(gdscdrugs.DRUG_ID(i)) = gdscdrugs.DRUG_NAME{i};
    end

    drugresponse = readtable('data/v17_fitted_dose_response.csv');
    drugresponse = drugresponse(ismember(drugresponse.DRUG_ID, cell2mat(keys(drugIDmap))), :);
    drugnames = values(drugIDmap, num2cell(drugresponse.DRUG_ID));
    drugresponse.DRUG_NAME = upper(drugnames);
    drugresponse.AA = 1 - drugresponse.AUC;
    
    
    % Remove profiles for cell lines that cannot be matched with the cell
    % line data
    idxnav = ~ismember(drugresponse.COSMIC_ID, gdsc_cells.CELL_ID);
    drugresponse(idxnav, :) = [];
    
    % Translate CellLineName in Profile Table to standardized name
    idxcells = arrayfun(@(x) find(gdsc_cells.CELL_ID==x), drugresponse.COSMIC_ID);
    drugresponse.Cell_Name = gdsc_cells.CELL_NAME(idxcells);
    
    % Remove profiles for cell lines without molecular data
    idxnav = ~ismember(drugresponse.Cell_Name, ds.cellNames);
    drugresponse(idxnav, :) = [];
    
    idxcells = cell2mat(cellfun(@(x) find(strcmp(x, ds.cellNames)), drugresponse.Cell_Name, 'UniformOutput', 0));
    
    alldrugs = unique(upper(drugresponse.DRUG_NAME));
    IC50Mat = nan(length(alldrugs), length(ds.cellNames));
    AAMat = nan(length(alldrugs), length(ds.cellNames));
    for i=1:length(alldrugs)
        drug = alldrugs{i};
        ic50sel = drugresponse.LN_IC50(strcmpi(drugresponse.DRUG_NAME, drug));
        aasel = drugresponse.AA(strcmpi(drugresponse.DRUG_NAME, drug));
        idxcellsel = idxcells(strcmpi(drugresponse.DRUG_NAME, drug));
        IC50Mat(i, idxcellsel) = ic50sel';
        AAMat(i, idxcellsel) = aasel';
    end
    
    ds.allDrugs = alldrugs;
    ds.IC50Mat = IC50Mat;
    ds.AAMat = AAMat;
end
