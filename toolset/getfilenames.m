function filenames = getfilenames(tdir, extension, opts)
% get file names of desired file types within target directory
% Oveis Jamialahmadi, University of Gothenburg, March 2022.
% 
% @14NOV2022: 'dir' option was added to return folder names in 'tdir'.

arguments
    tdir {mustBeFolder}
    extension {mustBeText, mustBeVector} = "*"
    opts.fullpath (1,1) logical = false
    opts.dir (1,1) logical = false % to return folder names. If true, ignores 'extension'
end

if opts.dir 
    % only return folder names
    filenames = struct2table(dir(tdir));
    filenames(~filenames.isdir, :) = [];
    filenames(ismember(filenames.name, {'.', '..'}), :) = [];
    if opts.fullpath
        filenames = fullfile(filenames.folder, filenames.name);
    else
        filenames = filenames.name;
    end
    filenames = string(filenames);
    return
end

extension = string(extension);
extension = regexprep(extension, "^\.", "");
ext = "*." + extension;

filenames = cell(numel(extension), 1);
for i = 1:numel(extension)
    tmpdir = struct2table(dir(fullfile(tdir, ext(i))));
    tmpdir(tmpdir.isdir, :) = [];
    filenames{i} = string(tmpdir.name);
    if ~isempty(filenames{i})
        filenames{i} = natsort(filenames{i});
        if opts.fullpath
            filenames{i} = fullfile(tdir, filenames{i});
        end
    end
end

emptyIdx = cellfun(@isempty, filenames);
filenames(emptyIdx) = [];
emptyExtension = extension(emptyIdx);
extension(emptyIdx) = [];

if ~isempty(extension)  
    filenames = cell2struct(filenames, matlab.lang.makeValidName(extension));
end

if ~isempty(emptyExtension)
    for i = 1:numel(emptyExtension)
        filenames.(matlab.lang.makeValidName(emptyExtension(i))) = [];
    end
end


end % END