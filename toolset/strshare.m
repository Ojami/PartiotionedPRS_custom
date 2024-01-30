function shared = strshare(str, opts)
% shared part of two strings.
% example: str = ["this.is.first", "this.is.second"]
%           shared = "this.is."

% Oveis Jamialahmadi, Sahlgrenska Academy, Jan 2022.

arguments
    str {mustBeText}
    opts.pattern (1,1) logical = false % create a pattern (e.g. file%d_ext)
    opts.patternstr {mustBeTextScalar} = '\d+'; % custom pattern to match
end

if numel(str) < 2
    error('strshare:input string must have 2 members')
end

str = string(str);
str = str(1:2);
if isrow(str); str = str.'; end

strc = char(str);

if opts.pattern
    num1 = regexp(str(1), opts.patternstr, 'match');
    num2 = regexp(str(2), opts.patternstr, 'match');
    if numel(num1) ~= numel(num2)
        error('strshare: input strings don''t have the same pattern!')
    end
    diffIdx = num1 ~= num2;
    shared = char(str(1));
    [start, stop] = regexp(shared, "("+opts.patternstr+")", 'start', 'end');
    for i = 1:numel(start)
        if diffIdx(i)
            offset = numel(start(i):stop(i)) - 2;
            if offset < 0
                idx1 = start(i); idx2 = stop(i);
                offset = -offset;
                start = start + offset; stop = stop + offset;
                shared = [shared(1:idx2), blanks(offset), shared(idx2+1:end)];
                shared(idx1:stop(i)) = '%d';
            else
                shared(start(i):stop(i)) = ['%d', blanks(offset)];
            end
        end
    end
    shared = erase(string(shared), whitespacePattern);
    return
    % Deprecated @03/13/2022
    if numel(unique(str.strlength)) > 1 % both strings much have the same length
        error('strshare: both strings must have the same length!')
    end
    idx = unique([0, find((strc(1, :) == strc(2, :)) == 0)]);
    if idx(end) == size(strc, 2)
        addLast = true;
    else
        addLast = false;
        idx(end + 1) = size(strc, 2) + 1;
    end
    shared = cell(1, numel(idx) - 1);
    for i = 1:numel(idx)-1
        if addLast
            shared{i} = [strc(1, (idx(i)+1):(idx(i+1)-1)), '%d'];
        else
            if i == (numel(idx) - 1)
                shared{i} = strc(1, (idx(i)+1):(idx(i+1)-1));
            else
                shared{i} = [strc(1, (idx(i)+1):(idx(i+1)-1)), '%d'];
            end
        end
    end
    shared = string(horzcat(shared{:}));
else
    sharedidx = find(diff((strc(1, :) == strc(2, :))) == -1, 1, 'first');
    shared = string(str{1}(1:sharedidx));
end

end % END