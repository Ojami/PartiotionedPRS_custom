function makeFUMAfiles
clc

% calls ss2fuma function to generate summary stat files to be uploaded on
% FUMA web server

% get summary stat paths and files ----------------------------------------
pth = fileparts(fileparts(pwd));
dir = getfilenames(pth, dir=true);
dir(~startsWith(dir.lower, ["liverironcorrected", "protondensity"])) = [];

ct = 1;
for k = 1:numel(dir)
    tmp = getfilenames(fullfile(pth, dir(k)), "txt", fullpath=true).txt;
    tmp(~endsWith(tmp, "sumstat.step2.txt")) = [];
    if isempty(tmp), continue; end % no summary stat file

    files.path(ct, 1) = tmp;
    files.tag(ct, 1) = dir(k);
    ct = ct + 1;
end

% for lead independent variants
cojo = load(fullfile(fileparts(pwd), '\COJO\cojo_files_10mb\cojo_merged.mat')).cojo;

% call ss2fuma
for k = 1:numel(files.path)
    tab = cojo(cojo.Trait == files.tag(k), :);
    tab(tab.p > 5e-8, :) = [];
    ss2fuma(files.path(k), output=files.tag(k), leadsnps=tab)
end

end % END