function plotAccBar(savedir, plotname, bardata, druglabel, legendlabel, plottitle, split)
% Plot accuracy bar for all drugs to compare between our methods and
% baseline 


    if(nargin<7)
        split = 0;
    end
    
    % Sort based on the difference between second and first 
    difdata = bardata(:, 2) - bardata(:, 1);
    [~, idxsort] = sort(difdata);
    bar(bardata(idxsort, :));
    druglabel = druglabel(idxsort);
    
    % IF split, just call plotAccBar twice
    if(split==1)
        idx1 = 1:floor(length(druglabel)/2);
        idx2 = floor(length(druglabel)/2)+1:length(druglabel);
        plotname1 = strcat(plotname,'_1');
        plotname2 = strcat(plotname,'_2');
        plotAccBar(savedir, plotname1, bardata(idx1, :), druglabel(idx1), legendlabel, plottitle, 0);
        plotAccBar(savedir, plotname2, bardata(idx2, :), druglabel(idx2), legendlabel, plottitle, 0);
        return;
    end
    
% 
    
    set(gca, 'XTickLabel', druglabel); set(gca, 'XTick', 1:length(druglabel)); xtickangle(45);
    legend(legendlabel);
    g=gcf;
    g.Position = [0 0 1920 1080];
    filename = sprintf('%s%s.png', savedir, plotname);
    print(filename, '-dpng');
    title(plottitle);
    close all;
end
