function [geno, inopts] = getbulkgeno(snpList, chrList, opts)
% reads raw genotypes for variants in string array snpList with the
% corresponding chromosome in chrList (double or string).
% OPTIONAL:
%       - dosage: true (default, use dosage) or false (hard calls).
%       - minorAllele to prune variants sharing the same variant ID but
%         have different minor allele. Must correspond to snpList
% OUTPUT: geno, a struct containing raw genotypes each stored in a distinct
% field.
% Oveis Jamialahmadi, GU, Jan 2020.
% 
% [EDITED]: infoscore flag was added to select the variants based on their
%           imputation INFO score. 05/14/2021.
% [EDITED]: parallel flag was added to support for large number of
%           variants. 05/14/2021
% [EDITED]: verbose flag was added. 30/09/2021.
% [EDITED] : 'swap' field was added back to 'geno' struct. 'swap' is aliase
%           'SWAP_FLAG' field returned from bgenreader function containing
%           information on variants with flipped (swapped) alleles. In
%           original implementation of bgenreader, variants are swapped so
%           that allele 2 always be minor allele. So, this 'swap' field can
%           be useful when comparing to original bgen files. 27/01/2022.
% [EDITED]: a bug was fixed with parallel read. 27/01/2022.
% [EDITED]: 'datatype' option was added with default as being 'double'.
%           'single' reduces memory overhead and improves reading speed.
%           27/01/2022.
% [EDITED]: 'getgeno' function was modified due to the change made to
%           bgenreader function with I_A == 1 and 'GP' flag. 16/02/2022. 
% [EDITED]: 'eid' was added for sample subsetting. 'eid' should be a vector
%           of samples for whom genotype (dosage) calls must be read.
%           17/02/2022.
% [EDITED]: a bug was fixed in sample subsetting. 03/02/2022.
% [EDITED]: 'tall' flag was added to read large genotype calls lazily to
%           fit into memory as tall arrays. 18/05/2022.
% [EDITED]: 'merge' flag was added (default is false) to merge all calls
%           over different chromosomes into a single struct. Calls for
%           different eids between chromosomes are filled with NaN values
%           in the final call matrix. 18/05/2022.
% [EDITED]: 'tag' option was added. These tags will be included as a new
%           field in the 'geno' struct correponding to variants in snpList.
%           06/07/2022.


arguments
       snpList {mustBeText}
       chrList {mustBeText}
       opts.dosage (1,1) logical = true % dosage or hard genotype calls
       opts.parallel (1,1) logical = false % use parallel toolbox to fetch variants (useful when there are large number of variants).
       opts.minorAllele {mustBeText} = ""
       opts.model {mustBeMember(opts.model, ["add", "dom", "rec"])} = "add"
       opts.home {mustBeFolder} = "D:\Imputed" % must only contains unique BGEN files for each chr
       opts.infoscore (1,1) double {mustBeInRange(opts.infoscore, 1e-6, 1)} = 0.7 % INFO score cut-off: variants < INFO score are removed.
       opts.chunk (1,1) double = 100 % number of snps read per each parfor loop (parallel mode only)
       opts.optionsOnly (1,1) logical = false % don't fetch variant calls, only information about files (for external software like regenie).
       opts.verbose (1,1) logical = true
       opts.datatype {mustBeTextScalar, mustBeMember(opts.datatype, ["single", "double"])} = "double"
       opts.eid {mustBeVector} % samples for whom genotype calls must be read
       opts.tall (1,1) logical = false % read memory friendly genotype calls as tall arrays. Only works in parallel.
       opts.merge (1,1) logical = false % to merge all chromosomes into one single struct
       opts.tag {mustBeText, mustBeVector} = "" % should have the same size as snpList
end

fclose('all');

% sample subsetting (see bgenreader help)
if ~isfield(opts, 'eid')
    opts.eid = [];
end

% variant tags
opts.tag = string(opts.tag);
opts.tag(ismissing(opts.tag) | opts.tag == "") = [];

% set variants options
inopts.imputedFlag = true;
inopts.dosage = opts.dosage;
inopts.verbose = false;

chrList = string(chrList);
chrom = natsort(unique(chrList));
geno = struct;

% find bgen files
bgenfiles = dir(fullfile(opts.home, '*.bgen'));
opts.checkbed = false;
if isempty(bgenfiles) && opts.optionsOnly % check bed files also
    bgenfiles = dir(fullfile(opts.home, '*.bed'));
    opts.checkbed = true;
end

if isempty(bgenfiles)
    error('no bgen file found within %s!', opts.home)
end
bgenfiles = natsort({bgenfiles.name}');

if ~opts.checkbed
    % further check if bgi/sample files exist
    samplefiles = getfilenames(opts.home, '.sample').sample;
    bgifiles = getfilenames(opts.home, '.bgi').bgi;
    if isempty(samplefiles) || isempty(bgifiles)
        error('no bgi/sample files found in %s!', opts.home)
    end
    bgifiles = regexprep(bgifiles, '.bgen.bgi$', '.bgen');
    samplefiles = regexprep(samplefiles, '.sample$', '.bgen');
    if ~isempty(setdiff(bgenfiles, bgifiles)) || ~isempty(setdiff(bgenfiles, samplefiles))
        error('bgi/sample files do not correspond to bgi files!')
    end
end

% construct bgen file string
if any(contains(lower(bgenfiles), ["_cx_", "_cy_", "_cxy_"])) % sex chromosomes
    chrPos = strshare(bgenfiles, 'pattern', true, 'patternstr', "(?<=_c)(.*?)(?=_)");
else
    chrPos = strshare(bgenfiles, 'pattern', true);
end
chrPos = replace(chrPos, "%d", "%s"); 

if opts.optionsOnly
    inopts.bgenfile = fullfile(opts.home, compose(chrPos, chrom))';
    inopts.chrPos = chrPos;
    inopts.chrom = chrom;
    return % don't fetch variant calls
end

if opts.verbose
    % create a loading bar 
    progressBar = floor(linspace(1, numel(chrom)+1, 11));
    if ~opts.parallel
        fprintf('getting raw genotypes [           ]')
    end
end

for i = 1:numel(chrom)
    idx = chrList == chrom(i);
    if ~any(idx); continue; end
    inopts.bgenfile = fullfile(opts.home, string(sprintf(chrPos, chrom(i))));
    
    % check for sex chromosomes. 
    sexchr = contains(inopts.bgenfile, wildcardPattern + ("X"|"Y"|"XY"|"23"|"24"|"PAR") + "_");
    opts.sexchr = false;
    if sexchr || ~isempty(opts.eid) % extract males and/or sample subsetting
        opts.N = bgenheader(inopts.bgenfile, 'sample', false, 'verbose', false);
        fid = fopen(regexprep(inopts.bgenfile, ".bgen$", ".sample"), 'r');

        if sexchr % read also last column (sex) in sample file
            samfile = textscan(fid, '%s %*s %*s %s %*[^\n]');
        else
            samfile = textscan(fid, '%s %*[^\n]');
        end

        fclose(fid);
        if opts.N > numel(samfile{1})
            error('sample/bgen fies mismatch!')
        end
        diff_N = numel(samfile{1}) - opts.N;
        if ~isempty(opts.eid) % subset samples
            opts.samples = double(string(samfile{1}(diff_N+1:numel(samfile{1}))));
            opts.samples = ismember(opts.samples, opts.eid);
            if any(ismissing(opts.samples)), opts.eid = []; end
        end
                
        if sexchr % add 'sex' information to 2nd column of 'samples'
            opts.sexchr = true;
            sexinfo = double(string(samfile{2}(diff_N+1:numel(samfile{2}))));
            
            % 1=male and 2=female in *.sample files; convert to logical
            % (true=male, false=female)
            sexinfo(sexinfo == 2) = 0;
            sexinfo = logical(sexinfo);

            if ~isfield(opts, 'samples') || any(ismissing(opts.samples))
                opts.samples = true(numel(sexinfo), 1);
            end
            opts.samples = [opts.samples, sexinfo];
            
        end
    end

    if opts.parallel
        if opts.verbose
            fprintf('getting raw genotypes for chrom: %s (of %d)\n', chrom(i), numel(chrom))
        end
        
        snpChunk = snpList(idx);
        chunkIdx = unique([1:opts.chunk:numel(snpChunk), numel(snpChunk)+1]);
        if numel(chunkIdx) < 3 % use bgenread internal parfor when chunk numel is 1 and parallel is true
            inopts.numWorkers = inf; % this doesn't affect bgenreader workers, it's only to provoke internal bgenreader parfor
            Nworkers = 0;
        else
            inopts.numWorkers = 0;
            Nworkers = gcp('nocreate');
            if isempty(Nworkers)
                Nworkers = parpool("Processes");
            end
            Nworkers = Nworkers.NumWorkers;
        end

        [callsT, bed] = deal({});
        verbose = opts.verbose;
        if opts.verbose; ptime = tic; end
        parfor (j = 1:numel(chunkIdx)-1, Nworkers)
            inchunk = snpChunk;
            if verbose
                fprintf('fetching chunk: %d (of %d)\n', j, numel(chunkIdx)-1)
            end
            [callsT{j}, bed{j}] = getgeno(inchunk(chunkIdx(j):chunkIdx(j+1)-1), inopts, opts);
        end

        if verbose
            gettime = toc(ptime);
            fprintf('\b [done in %.2f seconds]\n', gettime);
        end
        calls = struct;
        [calls.I_A, calls.SWAP_FLAG] = deal([]);
        calls.SWAP_FLAG = logical(calls.SWAP_FLAG);
        [calls.bim.a1, calls.bim.a2, calls.bim.chr, calls.bim.pos, calls.bim.snp] = deal([]);
        calls.fam = callsT{1}.fam;
        for j = 1:numel(callsT)
            calls.I_A = [calls.I_A; callsT{j}.I_A];
            calls.SWAP_FLAG = [calls.SWAP_FLAG; callsT{j}.SWAP_FLAG];
            calls.bim.a1 = [calls.bim.a1; callsT{j}.bim.a1];
            calls.bim.a2 = [calls.bim.a2; callsT{j}.bim.a2];
            calls.bim.chr = [calls.bim.chr; callsT{j}.bim.chr];
            calls.bim.pos = [calls.bim.pos; callsT{j}.bim.pos];
            calls.bim.snp = [calls.bim.snp; callsT{j}.bim.snp];
            callsT{j} = [];
        end
        bed = horzcat(bed{:});
    else
        [calls, bed] = getgeno(snpList(idx), inopts, opts);
    end

    if ~isnumeric(calls.fam) || any(isinteger(calls.fam))
        if ~all(isnan(double(calls.fam)))
            calls.fam = double(calls.fam); 
        end
    end
    calls.bim.a1 = string(calls.bim.a1);
    calls.bim.a2 = string(calls.bim.a2);
    
    if ~isempty(opts.tag)
        [tidx1, tidx2] = ismember(calls.bim.snp, snpList); tidx2(tidx2 < 1) = [];
    end

    if opts.merge % merge all calls from different chromosomes into a single struct
        geno.eid{i, 1} = calls.fam;
        geno.bed{i, 1} = bed;
        geno.I_A{i, 1} = calls.I_A;
        geno.swap{i, 1} = calls.SWAP_FLAG;  
        geno.snp{i, 1} = calls.bim.snp;
        geno.chr{i, 1} = calls.bim.chr;
        geno.pos{i, 1} = calls.bim.pos;
        geno.a1{i, 1} = calls.bim.a1;
        geno.a2{i, 1} = calls.bim.a2;
        if ~isempty(opts.tag)
            geno.tag{i, 1} = strings(numel(calls.bim.snp), 1);
            geno.tag{i, 1}(tidx1) = opts.tag(tidx2);
        end
    
    else      
        geno.("chr" + chrom(i)).eid = calls.fam;
        geno.("chr" + chrom(i)).bed = bed;
        geno.("chr" + chrom(i)).I_A = calls.I_A;
        geno.("chr" + chrom(i)).swap = calls.SWAP_FLAG;  
        geno.("chr" + chrom(i)).snp = calls.bim.snp;
        geno.("chr" + chrom(i)).chr = calls.bim.chr;
        geno.("chr" + chrom(i)).pos = calls.bim.pos;
        geno.("chr" + chrom(i)).a1 = calls.bim.a1;
        geno.("chr" + chrom(i)).a2 = calls.bim.a2;
        if ~isempty(opts.tag)
            geno.("chr" + chrom(i)).tag = strings(numel(calls.bim.snp), 1);
            geno.("chr" + chrom(i)).tag(tidx1) = opts.tag(tidx2);
        end
    end

    clear calls bed
    
    if ~opts.parallel && opts.verbose
        progressTxt = [repmat('=', 1, sum(i >= progressBar)),'>',...
                repmat(' ', 1, 10-sum(i >= progressBar))];
        fprintf(repmat('\b', 1, 12))
        fprintf('%s]', progressTxt)
    end
end

if opts.merge
    eid = unique(vertcat(geno.eid{:})); % union of all eids
    [idx, idx2] = cellfun(@(x)ismember(eid, x), geno.eid, 'uni', false);
    geno.eid = eid;
    idx = horzcat(idx{:});
    idx2 = horzcat(idx2{:});
    for i = 1:numel(geno.bed)
        tmp = nan(numel(eid), size(geno.bed{i}, 2), opts.datatype);
        tmp(idx(:, i), :) = geno.bed{i}(idx2(idx(:, i), i), :); 
        geno.bed{i} = tmp;
    end
    geno.bed = horzcat(geno.bed{:});
    fis = string(setdiff(fieldnames(geno), ["eid", "bed"]));
    for i = 1:numel(fis), geno.(fis(i)) = vertcat(geno.(fis(i)){:}); end
end

if opts.verbose
    fprintf('\n')
end

% remove duplicates
if all(opts.minorAllele == "") || opts.merge
    return
else
    opts.minorAllele = string(opts.minorAllele);
end

for j = 1:numel(chrom)
    dups = duplicates(geno.("chr" + chrom(j)).snp);
    if ~isempty(dups)
        for i = 1:numel(dups)
            idx = find(geno.("chr" + chrom(j)).snp == dups(i));
            fix = ismember(geno.("chr" + chrom(j)).a2(idx), ...
                opts.minorAllele(snpList == dups(i)));
            idx = idx(~fix);
            geno.("chr" + chrom(j)).snp(idx) = [];
            geno.("chr" + chrom(j)).a1(idx) = []; 
            geno.("chr" + chrom(j)).a2(idx) = [];
            geno.("chr" + chrom(j)).chr(idx) = [];
            geno.("chr" + chrom(i)).pos(idx) = [];
            geno.("chr" + chrom(i)).I_A(idx) = [];
            geno.("chr" + chrom(i)).swap(idx) = [];
            geno.("chr" + chrom(j)).bed(:, idx) = [];
            if ~isempty(opts.tag)
                geno.("chr" + chrom(j)).tag(idx) = [];
            end
        end
    end
end

end % END

%% subfunctions ===========================================================
function [calls, bed] = getgeno(snpChunk, inopts, opts)
inopts.datatype = opts.datatype;
% calls = getrsid(snpChunk, inopts);
if inopts.dosage
    gpType = "GP";
else
    gpType = "GT";
end
if opts.parallel
    numWorkers = inopts.numWorkers; % parallel pool cannot be started from a worker
else
    numWorkers = 0;
    opts.tall = false;
end

if isempty(opts.eid) && ~opts.sexchr
    calls = bgenreader(inopts.bgenfile, "varInfo",  snpChunk, ...
            "datatype", inopts.datatype, "gpType", gpType, ...
            "verbose", inopts.verbose, "numWorkers", numWorkers, ...
            "tall", opts.tall);
else % sample subsetting
    calls = bgenreader(inopts.bgenfile, "varInfo",  snpChunk, ...
        "datatype", inopts.datatype, "gpType", gpType, ...
        "verbose", inopts.verbose, "numWorkers", numWorkers, ...
        "N", opts.N, "samples", opts.samples, "tall", opts.tall);
end

if opts.tall && ~istall(calls.bed)
    calls.bed = tall(calls.bed); % happens with chunk reading
end

rmIdx = find(calls.I_A < opts.infoscore); % remove those below this INFO score
calls.I_A(rmIdx) = []; calls.SWAP_FLAG(rmIdx) = [];
fnames = fieldnames(calls.bim);
for i = 1:numel(fnames)
    calls.bim.(fnames{i})(rmIdx) = [];
end
calls.bed(:, [2.*rmIdx-1; 2.*rmIdx]) = [];

if opts.datatype == "single"
    bed = zeros(numel(calls.fam), numel(calls.I_A), 'single');  
else
    bed = zeros(numel(calls.fam), numel(calls.I_A));   
end

% 'dosage': Estimated Alternate Allele Dosage or expected genotype counts : [P(0/1)+2*P(1/1)]
if opts.dosage
    PPed = 1:numel(calls.I_A);
    if opts.model == "add"
        bed = calls.bed(:, PPed.*2-1)+ 2.*calls.bed(:, PPed.*2);
    elseif opts.model == "dom"
        bed = calls.bed(:, PPed.*2-1)+ calls.bed(:, PPed.*2);
    elseif opts.model == "rec"
        bed = calls.bed(:, PPed.*2);
    end

else % hard calls
    bed = calls.bed(:, 1:2:end) + calls.bed(:, 2:2:end); % does not remove missing calls
    if opts.model == "dom"
        bed(bed > 1) = 1;
    elseif opts.model == "rec"
        bed(bed == 1) = 0;
        bed(bed > 0) = 1;
    end
end

calls = rmfield(calls, 'bed');

end