function str = truncateStr(str, n, opts)
% truncates string str so that max length of characters be equal or less
% than 'n'. str and n can be vector

% @07SEP2023: 'dot' option was added to add "..." to end of truncated
% strings.

arguments
    str {mustBeText}
    n {mustBeNonempty, mustBeNumeric} = 31 % max length of str vector
    opts.dot (1,1) logical = false
end

str = string(str);
if numel(str) > 1 && numel(n) == 1 
    n = repmat(n, numel(str), 1);
end

if numel(str) ~= numel(n)
    error('truncateStr::str and n must have the same size!')
end

for k = 1:numel(str)
    slen = str(k).strlength;
    ntmp = min(n(k), slen);
    str(k) = str(k).extractBefore(ntmp + 1);

    if n(k) < slen && opts.dot
        str(k) = str(k) + "...";
    end
end

end % END