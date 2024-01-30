function [infile, toolFlag, index] = readGWASfile(infile, opts)

% readGWASfile reads/filters GWAS summary statistics files from different
% tools (now supports PLINK2/SAIGE/BOLT-LMM/REGENIE) and outputs a table
% with column names alike. Input file can be delimited or MATLAB file. 
% 
% Oveis Jamialahmadi, October 2021. Sahlgrenska Academy, GU.
% 
% @30/12/2021: some bugs were fixed.
% @11/01/2022: a bug was fixed.
% @12/01/2022: 'n' flag was added to read sample size column. 
% @24/01/2022: 'index' and 'win' options were added. Only applicable when
%               infile is not a MAT file
% @03/02/2022: a bug was fixed with 'light' flag and mat infile.
% @02/07/2022: some bugs were fixed.
% @14MARCH2023: 'tall' flag was added. In case of 'parallel' being true,
%               'tall' flag returns the tall table before gathering. This
%               is useful for cases when there are other downstream
%               filtering to be applied to the table.
%               'legacy' flag was added (default is true). If false, uses
%               a newer set of column names (annotate -> snp, afreq ->
%               a1freq).

arguments
    infile {mustBeFile}
    opts.light (1,1) logical = true % excludes variants with P < opts.
    opts.p (1,1) double {mustBeInRange(opts.p, 0, 1)} = 0.05 % P value cutoff for 'light' option
    opts.parallel (1,1) logical = false % use tall/gather for reading 'infile' summary stat file
    opts.full (1,1) logical = false % if set to true reads also AF, beta and se columns
    opts.n (1,1) logical = false % read also column N (if available)
    opts.index = "" % index variants for which only variants around them will be kept (using opts.win)
    opts.win (1,1) double = 500 % in kbp. Window around the index variants 
    opts.workers (1,1) = nan
    opts.tall (1,1) logical = false % to return as a tall table (true). Works only if 'parallel' is true.
    opts.legacy (1,1) logical = true % use old column names, or use 'snp', and 'a1freq' for column names (false).
    opts.maxIndexSize (1,1) double = 100 % max number of index variants to be filtered in parallel mode. If 'index' size exceeds this value, filtering is done after gathering the tall ds.
end

% check parallel pool
if opts.parallel && isempty(gcp('nocreate'))
    if isnan(opts.workers)
        % number of physical cores usually gives the best efficiency for
        % regression analyses. Although this has not been thoroughly
        % tested.
        maxNumCompThreads('automatic'); % reset (only if was changed in current session)
        opts.workers = maxNumCompThreads; 
    end
    parc = parcluster('Processes');
    parc.NumWorkers = opts.workers;
    parpool(parc); 
end

% check if input is a mat file
if endsWith(infile, '.mat')
    opts.mat = true;
else
    opts.mat = false;
end

if opts.mat
    infile = load(infile);
    sname = fieldnames(infile);
    infile = infile.(sname{1});
    header = string(infile.Properties.VariableNames);
else
    
    infileChar = infile;
    [~, ~, infileExt] = fileparts(infile);
    infile = tabularTextDatastore(infile, 'TextType', 'string',...
        'FileExtensions', infileExt, 'VariableNamingRule', 'preserve');
    if numel(infile.SelectedVariableNames) < 3
        infile.Delimiter = '\t';
    end

    fid = fopen(infileChar, 'r');
    header = fgetl(fid); header = split(string(header));
    fclose(fid);
end

if any(contains(header, "LOG10P")) || any(header == "ID")
    toolFlag = "REGENIE";
elseif any(header == "A1FREQ")
    toolFlag = "BOLT";
elseif any(header == "AF_Allele2")
    toolFlag = "SAIGE";
elseif any(header == "A1_FREQ")
    toolFlag = "PLINK";
elseif opts.mat
    toolFlag = "MAT";
else
    error('readGWASfile::unknownFormat:if file has a custom format, save it as a MAT file first!')
end

if ~opts.mat
    
    if toolFlag == "REGENIE"

        extra_col = strcmpi(infile.VariableNames, 'extra');
        if any(extra_col)
            extra_col = find(extra_col);
            infile.TextscanFormats(extra_col(1)) = {'%q'};
        end
        infile.TreatAsMissing = 'NA';

        try % old REGENIE version
            infile.SelectedVariableNames = {'CHROM', 'GENPOS', 'ID', 'LOG10P.Y1'}; % only keep first phenotype
        catch
            try
                infile.SelectedVariableNames = {'CHROM', 'GENPOS', 'ID', 'LOG10P'};
            catch
                infile.SelectedVariableNames = {'CHROM', 'GENPOS', 'ID', 'P'};
            end
        end
    elseif toolFlag == "BOLT"
        try
            infile.SelectedVariableNames = {'SNP', 'CHR', 'BP', 'P_BOLT_LMM'};
        catch
            infile.SelectedVariableNames = {'SNP', 'CHR', 'BP', 'P_BOLT_LMM_INF'};
        end

    elseif toolFlag == "SAIGE"
        varNames = {'CHR', 'POS', 'rsid', 'p_value'};
        if ~isempty(setdiff(infile.VariableNames, varNames))
            % to see if already existing variable names are a subset of varNames
            infile.SelectedVariableNames = varNames;
        end

    elseif toolFlag == "PLINK" % PLINK  
        infile.SelectedVariableNames = {'x_CHROM', 'POS', 'ID', 'P'};
    end

    if any(contains(lower(header), 'gene'))
        annotateCols = header(contains(lower(header), 'gene'));
        annotateCols = annotateCols(1); % pick one!
        infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols).'];
    end

    if any(ismember(lower(header), 'clump'))
        annotateCols = header(ismember(lower(header), 'clump'));
        infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols).'];
    end

    if opts.n && any(ismember(lower(header), 'n'))
        Ncol = header(ismember(lower(header), 'n'));
        infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(Ncol(1))];
    end
    
    if opts.full % read extra columns 
        if any(contains(lower(header), 'beta'))
            annotateCols = header(contains(lower(header), 'beta'));
            infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols(1)).'];
        end
        if any(contains(lower(header), 'se'))
            annotateCols = header(ismember(lower(header), ["se", "standard_error", "standard error"]));
            infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols(1)).'];
        end
        if any(contains(lower(header), ["af_allele2", "a1freq", "a1_freq", "maf", "a2freq"]))
            annotateCols = header(contains(lower(header), ["af_allele2", "a1freq", "a1_freq", "maf", "a2freq"]));
            infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols(1)).'];
        end
        if any(ismember(lower(header), ["allele0", "allele1", "allele2", "a1", "a2"]))
            annotateCols = header(ismember(lower(header), ["allele0", "allele1", "allele2", "a1", "a2"]));
            infile.SelectedVariableNames = [infile.SelectedVariableNames, cellstr(annotateCols).'];
        end
    end
end

% unify column names
if opts.mat
    colNames = lower(string(infile.Properties.VariableNames));
else
    colNames = string(lower(infile.SelectedVariableNames));
    [~, selIdx] = ismember(infile.SelectedVariableNames, infile.VariableNames);
end

colNames(contains(colNames, ["p_bolt_","p_value", "pval", "p_adj"])) = "p";
colNames(contains(colNames, "log10p")) = "log10p";
colNames(ismember(colNames, "bp")) = "pos";
colNames(contains(colNames, "genpos")) = "pos";
colNames(contains(colNames, "chr")) = "chr";
if opts.full
    colNames(ismember(colNames, ["beta", "b", "beta_y1", "beta.y1"])) = "beta";
    if opts.legacy, freqcol = "afreq"; else, freqcol = "a1freq"; end
    colNames(ismember(colNames, ["af_allele2", "a1freq", "a1_freq", "maf", "a2freq"])) = freqcol;
    colNames(startsWith(colNames, "se")) = "se";
end

if any(ismember(colNames, ["gene", "nearest_genes"]))
    if any(colNames == "annotate")
        fprintf('WARNING::readGWASfile:annotate column already exists, will overwrite it!\n')
        colNames(colNames == "annotate") = "annotate1";
    end
    colNames(ismember(colNames, ["gene", "nearest_genes"])) = {'annotate'};
else
    snpcolidx = ismember(colNames, ["snp", "varid", "rsid", "rsids", "id"]);
    if any(colNames == "annotate") && any(snpcolidx)
        fprintf('WARNING::readGWASfile:annotate column already exists, will overwrite it!\n')
        colNames(colNames == "annotate") = "annotate1";
    end

    if sum(snpcolidx) > 1, snpcolidx = find(snpcolidx, 1); end
    if opts.legacy
        colNames(snpcolidx) = {'annotate'};
    else
        colNames(snpcolidx) = {'snp'};
    end
end

% check index variants ----------------------------------------------------
if ~istable(opts.index) && numel(string(opts.index)) < 2 && isfile(opts.index)
    try % matfile
        opts.index = load(opts.index);
        fi = fieldnames(opts.index);
        opts.index = opts.index.(fi{1});
    catch % otherwise
        opts.index = bfilereader(opts.index, 'parallel', opts.parallel, 'verbose', 'off');
    end
end

if istable(opts.index)
    indexcols = string(lower(opts.index.Properties.VariableNames));
    snpcol = find(ismember(indexcols, ["snp", "variant", "annotate", "id", "variantid"]));
    indexcols(snpcol(1)) = "snp";
    poscol = find(indexcols == "bp" | indexcols == "genpos" | startsWith(indexcols, "pos"));
    indexcols(poscol(1)) = "pos";
    chrcol = find(ismember(indexcols, ["chr", "chrom", "chromosome"]));
    indexcols(chrcol(1)) = "chr";
    opts.index.Properties.VariableNames = indexcols;
    opts.index = opts.index(:, ["chr", "snp", "pos"]);
else
    opts.index = table(opts.index, 'VariableNames', {'snp'}); 
end
opts.index(opts.index.snp == "", :) = [];
index = opts.index;

% read summary stat file --------------------------------------------------
if opts.mat
    infile.Properties.VariableNames = colNames;
    if any(contains(infile.Properties.VariableNames, 'log10p')) % REGENIE
        if opts.light
            infile(infile.log10p < -log10(opts.p), :) = []; 
        end
        infile.p = 10.^-infile.log10p;
        infile.log10p = [];
    else
        infile(infile.p > opts.p, :) = [];
    end
else
    infile.VariableNames(selIdx) = cellstr(colNames);
    infile.SelectedFormats(ismember(infile.SelectedVariableNames, {'id', 'annotate', 'snp'})) = {'%q'};
 
    if opts.parallel
        chrcol = strcmp(infile.SelectedVariableNames, "chr");
        infile.SelectedFormats{chrcol} = '%q'; % if X/Y have not been numbered as 23/24
        infile = tall(infile);
        if opts.light && any(contains(infile.Properties.VariableNames, 'log10p')) % REGENIE
            infile(infile.log10p < -log10(opts.p), :) = []; 
        elseif opts.light
            infile(infile.p > opts.p, :) = [];
        end

        if any(contains(infile.Properties.VariableNames, 'log10p')) % REGENIE
            infile.p = 10.^-infile.log10p;
            infile.log10p = [];
        end
        
        % a large index size slows down the filtering process; so, it's
        % usually faster to read whole dataset and then perform the
        % filtering.
        if height(opts.index) < opts.maxIndexSize
            if ~isempty(opts.index) % filer using index variants
                indexIdx = fliterIndexVariants(infile, opts);
                infile(~indexIdx, :) = [];
            end
        end

        matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
        if ~opts.tall
            infile = gather(infile);
        end
        matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
        
        % skip other steps since they would 'gather' the tall table
        if opts.tall, return; end

        % convert X/Y to 23/24
        checkChr = double(infile.chr);
        if any(isnan(checkChr))
            infile.chr(infile.chr == "X") = "23";
            infile.chr(infile.chr == "Y") = "24";
            infile.chr = double(infile.chr);
        else
            infile.chr = checkChr;
        end

        if height(opts.index) >= opts.maxIndexSize
            if ~isempty(opts.index) % filer using index variants
                indexIdx = fliterIndexVariants(infile, opts);
                infile(~indexIdx, :) = [];
            end
        end

        
    else
        infile = readall(infile);
        if opts.light && any(contains(infile.Properties.VariableNames, 'log10p')) % REGENIE
            infile(infile.log10p < -log10(opts.p), :) = []; 
        elseif opts.light
            infile(infile.p > opts.p, :) = [];
        end
        if any(contains(infile.Properties.VariableNames, 'log10p')) % REGENIE
            infile.p = 10.^-infile.log10p;
            infile.log10p = [];
        end

        if ~isempty(opts.index) % filer using index variants
            indexIdx = fliterIndexVariants(infile, opts);
            infile(~indexIdx, :) = [];
        end
    end
end

end % END

%% subfunction ============================================================
function idx = fliterIndexVariants(infile, opts)
opts.win = opts.win.*1e3; % kbp -> bp
if any(ismember(opts.index.Properties.VariableNames, 'chr')) && ...
        isstring(opts.index.chr) && ~any(isnan(double(opts.index.chr)))
    opts.index.chr = double(opts.index.chr);
end

if (istall(infile) && isaUnderlying(infile.chr, 'string')) || isUnderlyingType(infile.chr, "string")
    opts.index.chr = string(opts.index.chr);
end

if any(ismember(opts.index.Properties.VariableNames, {'pos'})) % chr must be present in this case
    for i = 1:numel(opts.index.snp)
        if i == 1
            idx = infile.pos >= (opts.index.pos(i) - opts.win) &...
                infile.pos <= (opts.index.pos(i) + opts.win) & ...
                infile.chr == opts.index.chr(i);
        else
            idx = idx | (infile.pos >= (opts.index.pos(i) - opts.win) &...
                infile.pos <= (opts.index.pos(i) + opts.win) & ...
                infile.chr == opts.index.chr(i));
        end
    end
else
    if opts.legacy, idname = "annotate"; else, idname = "snp"; end
    [~, varidx] = ismember(infile.(idname), opts.index.snp);
    varidx = find(varidx);
    for i = 1:numel(opts.index.snp)
        if i == 1
            idx = infile.pos >= (infile.pos(varidx(i)) - opts.win) &...
                infile.pos <= (infile.pos(varidx(i)) + opts.win) & ...
                infile.chr == infile.chr(varidx(i));
        else
            idx = idx | (infile.pos >= (infile.pos(varidx(i)) - opts.win) &...
                infile.pos <= (infile.pos(varidx(i)) + opts.win) & ...
                infile.chr == infile.chr(varidx(i)));
        end
    end
end

end