function ss2fuma(sumstat, opts)
% converts summary statistics file for FUMA web service. 
% The output of this function should be uploaded on https://fuma.ctglab.nl/
% 
% Oveis Jamialahmadi, University of Gothenburg, March 2023.

arguments
    sumstat {mustBeFile}
    opts.output {mustBeTextScalar}
    opts.leadsnps {mustBeA(opts.leadsnps, 'table')} % lead snp file
    opts.map {mustBeA(opts.map, 'table')} % mapping between old and new variable names with two columns: 1-old and 2-new
end

t1 = tic;
tmpname = getRandomName("fuma", 5);
sumstat = string(sumstat);
[pth, name] = fileparts(sumstat);
if pth == "", pth = pwd; end

if isfield(opts, 'output')
    [pth, name, ext] = fileparts(string(opts.output));

    if pth == "", pth = pwd; end

    if ~isfolder(pth)
        mkdir(pth)
    end
    output = fullfile(pth, name + ext);

else
    output = fullfile(pth, name + "_fuma");
end

if isfield(opts, 'leadsnps')
    tab = opts.leadsnps;
    h = colnames(tab);
    snpidx = find(ismember(h.lower, ["snp", "varid", "rsid", "rsids", "id"]), 1);
    chridx = find(contains(h.lower, "chr"), 1);
    posidx = find(ismember(h.lower, ["bp", "position", "pos"]) | contains(h.lower, "genpos"), 1);
    tab = tab(:, [snpidx, chridx, posidx]);
    writetable(tab, fullfile(pth, name + "_leadsnps.txt"), Delimiter="\t");
end

hdr = bfilereader(sumstat, "header", true, "summary", "firstline");

if any(contains(hdr.lower, "log10p")) % regenie
    pcol = find(contains(hdr.lower, "log10p"));

     % convert log10 P -> P
    runbash("awk 'NR!=1 {$"+ pcol(1) + "=10^-$" + pcol(1) + ...
        "}1' " + makeWSLpath(sumstat) + " > " + makeWSLpath(output), ...
        tmpname, "wait", true);
    sumstat = output;
end

% change column names
if isfield(opts, 'map')
    hdr_new = hdr;
    for k = 1:height(opts.map)
        [idx1, idx2] = ismember(hdr.lower, lower(opts.map.old));
        hdr_new(idx1) = opts.map.new(idx2);
    end
    hdr = hdr_new;
else
    hdr = createGWASheader(hdr);
end

if sumstat == output
    cmd = "sed -i '1s/.*/" + hdr.join(" ") + "/' '" + makeWSLpath(sumstat) + ...
        "'";
else
    cmd = "sed '1s/.*/" + hdr.join(" ") + "/' '" + makeWSLpath(sumstat) + ...
        "'  > '" + makeWSLpath(output) + "'";
end
runbash(cmd, tmpname, "wait", true);

% compress it
gzip(output)
delete(output)

fprintf('summary stat file %s is ready (done in %.0f sec)\n', output, toc(t1))

end %END