function [ quant] = quantnorm( raw_data, option )

% This function normalizes data from large scale data set
% columns contain the observations
% rows contain the different features
%   Detailed explanation goes here
% the function performs as a default a median normalization
% other options are listed below
%
% option={'median','mean','trimean','geomean','harmean'};
% 
% How it works?
% Quantile normalization is simple and effective
% 1 Step: sort all columns amd store the sorted data and the sort index 
% 2 Step: calculate the row 'median','mean','trimean','geomean','harmean'
% of the sorted data
% 3 Step: replace the original raw value with the matching quantile value
% 
% Ex: Example with and 6 x 4 array
% Input array
% 
% 668	1515  469	 772
% 1148	1631 1500	 593
% 439	931	  449	 441
% 1473	236	 1711	1041
% 500	254	  645	 733
% 992	438	  362	689

% Step 1: sorted raw data
% 439	 236	 362	 441
% 500	 254	 449	 593
% 668	 438	 469	 689
% 992	 931	 645	 733
% 1148	1515	1500	 772
% 1473	1631	1711	1041

% Step 2: median value for each quantile
%  400.5
%  474.5
%  568.5
%  832
% 1324
% 1552

% Step 3: replace the original raw data value with the matching quantile
% value. This is done with the resort index
% 3	5	3	5		568.5	1324	568.5	1324
% 5	6	5	2		1324	1552	1324	474.5
% 1	4	2	1		400.5	832	474.5	400.5
% 6	1	6	6		1552	400.5	1552	1552
% 2	2	4	4		474.5	474.5	832	832
% 4	3	1	3		832	568.5	400.5	568.5

% Job done

try
% Get the number of dimensions
cols=size(raw_data,2);


% Now normalize the data

[quant,quant_idx]=sort(raw_data);
[~,resort_idx]=sort(quant_idx);

% Check which normalization options was requested
if nargin==1
    normalizer=median(raw_data);
else
    switch lower(option)
        case 'median'
            normalizer=median(quant,2);
        case 'mean'
            normalizer=mean(quant,2);
        case 'trimmean'
            normalizer=trimmean(quant',10)';
            
        case 'harmmean'
            normalizer=harmmean(quant,2);
            
    end
end

dummy=repmat(normalizer,1,cols);
quant=dummy(resort_idx);


catch
    disp('Usage: provide rawdata')
    disp('Optionally provide the normalization option')
    disp('Options: median,mean,trimean,harmean')
end

