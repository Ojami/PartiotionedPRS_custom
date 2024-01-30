function pth = makeWSLpath(pth, toWin)
% converts win path to WSL path (use ispc alongside to decide to use this
% func) on mounted drives (mnt).
% Oveis Jamialahmadi, Jan 2020.
% 
% @03/19/2022: 'toWin' flag was added to convert back a wsl path to windows
%               path
% 
% TODO: WSL now has wslpath to convert win<->linux paths.

arguments
    pth {mustBeTextScalar}
    toWin (1,1) logical = false % to convert back a wsl to win path 
end

if ~ispc; return; end % already in unix/mac
pth = string(pth);
if fileparts(pth) == ""; return; end % exists in pwd

% note that the objective of this func is not to check if pth is a valid
% isfile/isfolder!

if toWin % wsl -> win
    if contains(pth, "\") % this is already a win path!
        return % should be silent?
    end
    pth = regexprep(regexprep(pth, '^/mnt/', ''), '^.*?(?=/)', ...
    '${upper($0)+":"}');
    pth = replace(pth, "/", "\");
else % win -> wsl
    if contains(pth, "/") % this is already a unix path!
        return % should be silent?
    end
    pth = regexprep(pth, '^.*?[:]', '${lower($0)}', 'once');
    pth = '/mnt/'+ strrep(pth, '\', '/'); % only works with wsl mounted drives
    pth = strrep(pth, ':/', '/');
    pth = strrep(pth, ":", "/"); % to deal with z:XXX instead of z:/XXXX
end

end % END