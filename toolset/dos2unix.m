function dos2unix(file, opts)
% a wrapper around dos2unix module (needs WSL)
% input: 
%   file: full path to file. The path can be either in unix (e.g.
%         /mnt/e/file) or windows format (e.g. E:\file).
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, March 2022.

arguments
    file {mustBeUnixFile}
    opts.verbose (1,1) logical = false
end

if ~ispc
    return
end

if checkCR(file) % check if carriage return is present
    if opts.verbose
        system("wsl sudo dos2unix " + makeWSLpath(file));
    else
        [~, ~] = system("wsl sudo dos2unix " + makeWSLpath(file));
    end
end
end % END

%% subfunctions ===========================================================
function mustBeUnixFile(file)
% validation function for unix file (equivalent to mustBeFile) inside
% windows

file = makeWSLpath(file, true); % convert back to win path
mustBeFile(file)

% directly check a win path instead of system cmd to wsl
% if ispc && contains(file, filesep) % should be converted to WSL path first
%     file = makeWSLpath(file);
% end
% 
% [doesntExist, ~] = system('wsl [ -f "' + string(file) + '" ]');
% if doesntExist
%     eid = 'WSL:notFound';
%     msg = 'File doesn''t exist!';
%     throwAsCaller(MException(eid, msg))
% end
end

%% ------------------------------------------------------------------------
function out = checkCR(file)
% first convert the WSL to Windows path

file = makeWSLpath(file, true); % convert back to win path
fid = fopen(file);
out = double(fgets(fid));
fclose(fid);

out = any(out == 13); % regexp(out, '\r')
end