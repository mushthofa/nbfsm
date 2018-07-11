function stdname = stdCellName(name)
    % Standardize cell names
    % 1. Remove all punctuations : space, _, -, (, ), /
    % 2. Capitalise

    %Black listed
    if(ismember(name, {'KM-H2';  'KMH-2'; 'T-T'; }))
        stdname = name;
        return;
    end

    if(strcmp(name, 'T.T'))
        stdname = 'T-T';
        return;
    end

    splittedname = strsplit(name,  {' ', '-', '_', '(', ')', '/', '.', ','}, 'CollapseDelimiters', 1);
    stdname = '';
    for i=1:length(splittedname) - 1
        stdname = strcat(stdname, splittedname{i});
    end
    stdname = strcat(stdname, splittedname{end});
    stdname = upper(stdname);
end









