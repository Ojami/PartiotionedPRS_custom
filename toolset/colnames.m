function [cols, idx] = colnames(tab, opts)
% returns column names of table.
% Oveis Jamialahmadi, University of Gothenburg, February 2023.
% 
% @18MAY2023: 'index' option was added to return only the columns
%             correponding to a set of index values. Note that this option
%             does not affect 'find' option, because 'find' is first taken
%             care of.
arguments
    tab
    opts.find {mustBeText, mustBeVector}
    opts.index {mustBeVector} % index of columns to be read
end

if ~istable(tab) && ~istall(tab)
    error("colnames:input tab must be table or a tall table!")
end

idx = [];
cols = string(tab.Properties.VariableNames);

if isfield(opts, 'find')
    idx = ismember(cols, opts.find);
end

if isfield(opts, "index")
    cols = cols(opts.index);
end

end % END