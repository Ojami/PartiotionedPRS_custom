function [dups, dupCounts, count, Sidx] = duplicates(S)
% detects duplicates in a string or cellstr vector.
% INPUT:
%       S: a string/cellstr vector.
% OUTPUTS:
%       dups:      duplicated values.
%       dupCounts: duplicate counts for each value in dups.
% Oveis Jamialahmadi, GU. Oct. 2020.

if nargin < 1
    error('input is missing!')
elseif ~any(size(S) == 1)
    error('only row/column vector!')
elseif ~isstring(S) && ~iscellstr(S)
    S = string(S);
end

[Sunique, Sidx, Suidx] = unique(S, 'stable');
count = histcounts(Suidx,[1:numel(Sunique), inf])';
dups = Sunique(count > 1);
dupCounts = count(count > 1);

end % END