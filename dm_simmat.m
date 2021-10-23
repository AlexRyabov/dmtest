function [S, D, Sfull] = dm_simmat(X, varargin)
%Alexey Ryabov, 2021
%Find similarity matrix beween columns (variables) of data matrix X and secles 10 largest entries in each column. 
%Rows of X correspond to observations and columns to variables.
%
% [S] = dm_simty(X) returns the cleaned matrix S where
% each column contains k_max largest entries, by default k_max = 10
% for the N by P data matrix X. Rows of X correspond to observations and columns to
% variables.
%
% [S, Sfull] = dm_simmat(X) also returns the uncleaned similarity matrix Sfull
%
% [...] = dm_simmat(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%     parameter name/value pairs to control the computation and handling of
%     special data types. Parameters are:
%
% 'Metric' -   Metric of the similarity matrix. Choices are:
%             'euclidean'; 'cityblock'; 'chebychev'; 'mahalanobis'; 'minkowski';
%             'cosine';  'spearman'; 'pearson'; 'kendall'; 'jaccard';'gaussian'
%             see pdist for their meaning. 
%             When metric is distance based, e.g. 'euclidean' 'cityblock' 
%             then similarity matrix S = 1/D, 
%             for correlation based metric  ('spearman'; 'pearson'; 'kendall') similarity S = (D + 1)/2,
%             for cosine, jaccard and gaussian based distance the similarity S = 1-D
% 'Norm' -    Normalization of input varaibles. Choices are:
%     'none' no normalization
%     'normalized'  Normilized along columns, so that sum(X.^2, 1) = 1
%     'zscore' each column standartized to have mean = 0 and std = 1, X = (X-mu)/sigma
% 'k_min' -   : the minimal number of entries in the cleaned similarity
% matrix, default k_min = 10
%
%   Example:
%  [Sclned, Sfull] = dm_simmat(X);
%   [Sclned] = dm_simmat(r, 'k_min', 10,  'Metric', 'cosine', 'Norm', 'normalized');
%
%  See also dm_dmit, dm_plot, pdist, zscore
% References:
%   [1]

%get paramters
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[n, p] = size(X);
if n == 0 || p == 0
    error('Input matrix X should not be empty');
end
internal.stats.checkSupportedNumeric('X',X); % not complex

% Parse arguments and check if parameter/value pairs are valid
%Norm  none/NORMALIZED/standardized
%Distance — Distance metric


paramNames = {'Metric',   'Norm', 'k_min'};
defaults   = {'spearman', 'normalized', 10};

[Metric, Normalization, k_min, sf, rest]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

% Validate String value for  Distance value
MetricNames = {'euclidean'; 'cityblock'; 'chebychev'; ...
    'mahalanobis'; 'minkowski'; 'cosine';  ...
    'spearman'; 'pearson'; 'kendall'; 'jaccard';'gaussian'};
Metric = internal.stats.getParamVal(Metric,MetricNames,...
    '''Similarity metric''');

% Validate String value for Normalization value
Normalizations = {'none'; 'normalized'; 'zscore'};
Normalization = internal.stats.getParamVal(Normalization,Normalizations,...
    '''Norm''');

% Validate integer value for 'k_min' option
if k_min > p
    warning('''k_min (%i)'' should be less or equal the number of columns in X', k_min)
end

% Validate the number of components option 'NumComponents'
if isempty(X)
    error(message('diffusion map:Input matrix X is empty'));
end

%Normalize
switch Normalization
    case 'normalized'
        fNorm = sqrt(sum((X).^2, 1));
        fNorm(fNorm==0) = 1;
        X = X./ fNorm;
    case 'standardized'
        X = zscore(X);
end

%Find the similarity matrix
switch Metric
    case {'cosine', 'jaccard'}
        D = squareform(pdist(X',Metric));
        Sfull = 1-D;
    case 'gaussian'  %Gaussian distance
        D = squareform(pdist(X','euclidean'));
        Sfull = exp(-D.^2);
    case {'pearson', 'spearman', 'kendall'}  %correlation distance.
        Sfull = (corr(X,'Type', Metric)+1)/2;
        D = 1-Sfull;
    otherwise %'euclidean', 'cityblock' ... 
        D = squareform(pdist(X',Metric));
        Sfull = 1./D;
end

%Cleaned similarity matrix
S = Sfull;

%set diagonal elements to zero
for iS = 1:size(S, 1)
    S(iS, iS) = 0;
end

%Remove infinite similatiry entiries
ind_inf = S == Inf;
if sum(ind_inf(:)) > 0
    [i_inf, j_inf] = find(ind_inf == 1);
    warning('The following pairs of variables have infinite similarity. We recommend either combining these variables or using a correlation, cosine, or gaussian similarity')
    display([i_inf, j_inf]);
    MaxSim = max(S(S(:) < Inf));
    S(S == Inf) = MaxSim;
    warning('Infinite similarity is set to max observed similarity of %f', MaxSim)
end

%clean the similarity matrix (leave k_min largest elements in each colum
B =(maxk(S,k_min));
S(S<repmat(B(end,:), size(S, 1), 1)) = 0;
S = max(S, S');
D(S==0) = 0;

end
