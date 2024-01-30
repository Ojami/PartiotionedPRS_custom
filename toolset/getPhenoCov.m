function [phenotypes, covariates] = getPhenoCov(opts)
% getPhenoCov processes locally saved trait structures (mat files), and
% returns two tables for phenotype(s) and covariate(s) (+ genomic PCA and
% array batch for UKBB). This function is basically implements the methods
% already integrated in gwasrunner function.
% 
% Oveis Jamialahmadi, Sahlgrenska Akademy, July 2021.
% 
% @15/10/2021: a bug was fixed with 'verbose' flag.
% @15/10/2021: 'defaultcovar' flag was added (default: true) to include
%               default covariates for UKBB analysis: first 10 genomic PC
%               and array batch. 15 October 2021.
% @14/04/2022: covariates with exclusion criteria ('exeid' field in pheno
%              struct) are now supported. 
% @13APR2023: 'legacy' flag was added to use either file names or tags for
%              variable names in the returned table.
% @14MAY2023: 'getbulkgeno' arguments were added to fetch genetic (imputed)
%             data.
% @14MAY2023: 'covarTab' was added to include additional covariates (as a
%             table).

arguments
    opts.trait {mustBeText} = "" % if left empty, UKBtrait_toolbox will fetch traits
    opts.qc (1,1) double {mustBeMember(opts.qc, 0:6)} = 0 % see getQCEID function for more details
    opts.qcinc (:, 1) double {mustBeNumeric} = [] % custom eids to be included; in this case intersect of qc and qcinc will be used for analyses.
    opts.phenopath {mustBeFolder} = "D:UKB_PHENO"
    opts.transform {mustBeMember(opts.transform, ["none", "log", "log10", "irnt"])} = "irnt" % transformation function for continuous variables.
    opts.covar {mustBeText} = ["Age", "Sex", "BMI"]; % list of covariates, it automatically adds PC1-10 and array batch.
    opts.covarTab {mustBeA(opts.covarTab, "table")} % additional covariates to be added (a table with covariates + 'eid' column).
    opts.wes (1, 1) logical = false % for WES analysis (array batch is removed from covariates)
    opts.removenan (1,1) logical = false % remove samples with missing phenotype (only effective is number of query traits is 1).
    opts.verbose (1,1) logical = true 
    opts.defaultcovar (1,1) logical = true % to automatically add PC1-10 and array batch to covariates.
    opts.legacy (1,1) logical = true % true: use 'tag' field for variable names; false: use file names for variable names.

    % genetic data
    opts.snp
    opts.chr
    opts.dosage (1,1) logical = true % dosage or hard genotype calls
    opts.model {mustBeMember(opts.model, ["add", "dom", "rec"])} = "add"
    opts.genohome {mustBeFolder} = "D:\Imputed" % must only contains unique BGEN files for each chr
end

% get QCed samples --------------------------------------------------------
qc_eid = getQCEID(opts.qc, opts.verbose);
if opts.verbose
    fprintf('\nsample N = %d\n', numel(qc_eid))
end
if ~isempty(opts.qcinc) % custom samples (eid)
    qc_eid = intersect(opts.qcinc, qc_eid);
    if opts.verbose
        fprintf('sample N after applying custom samples = %d\n', numel(qc_eid))
    end
    opts.qcinc = [];
end

% get covariates ----------------------------------------------------------
tfile = fullfile(opts.phenopath, opts.covar);
get_covars = struct;

if ~opts.legacy % use file names as variable names
    [~, covar_names_raw] = fileparts(tfile);
end

if opts.covar ~= ""
    covar_names = strings(numel(tfile), 1); % original covariate names
    for i = 1:numel(tfile)
        t = load(tfile{i});
        covar_names(i) = string(t.UKB_STRUCT_ALL.tag);
        if opts.legacy
            get_covars.(matlab.lang.makeValidName(t.UKB_STRUCT_ALL.tag)) = t.UKB_STRUCT_ALL; 
        else
            get_covars.(covar_names_raw(i)) = t.UKB_STRUCT_ALL;
        end
    end

    covariates = getcovar(get_covars, qc_eid, opts);
    get_covars = fieldnames(get_covars);
    extracov = "cov" + (numel(get_covars)+1:size(covariates, 2)).';
    get_covars = [get_covars; extracov];
    covar_names = [covar_names; extracov];

    %@14MAY2023: additional covariates
    if isfield(opts, "covarTab")
        ctab = rmmissing(opts.covarTab); 
        opts = rmfield(opts, "covarTab");
        [f1, f2] = ismember(qc_eid, ctab.eid); f2(f2<1) = [];
        cols = setdiff(colnames(ctab), "eid", "stable");
        covariates2 = nan(numel(qc_eid), numel(cols));
        covariates2(f1, :) = ctab{f2, cols};
        covariates = [covariates, covariates2];
        clear covariates2 ctab
        get_covars = [get_covars; cols'];
        covar_names = [covar_names; cols'];
    end

    if opts.verbose
        fprintf('input covariates: %s\n\n', strjoin(get_covars, ','))
    end
    
    if opts.removenan
        f_rem = any(isnan(covariates), 2);
        if any(f_rem)
            if opts.verbose
                fprintf('%d missing samples were removed from covar_matrix!\n', sum(f_rem))
            end
            covariates(f_rem, :) = [];
            qc_eid(f_rem, :) = [];
        end
    end
    
    covariates = array2table([qc_eid, covariates], 'VariableNames', ...
        ["eid"; string(get_covars)]);
    covariates.Properties.VariableDescriptions = ["eid"; covar_names];
        
else
    covariates = []; % no covariates are queried 
end

% get phenotypes ----------------------------------------------------------
% check traits
if opts.trait == ""
    opts.trait = getPheno(true, false, opts.phenopath);
else
    opts.trait = string(opts.trait);
end

[trait, trait_desc] = deal(strings(numel(opts.trait), 1));
bt = false(numel(opts.trait), 1);
phenotypes = nan(numel(qc_eid), numel(opts.trait));
notfoundIdx = false(numel(opts.trait), 1);
for n = 1:numel(opts.trait)
    
    if isstruct(opts.trait)
        D = opts.trait(n);
        opts.trait(n).eid = "";
        opts.trait(n).rawUKB = "";
    else
        dv = fullfile(opts.phenopath, opts.trait(n));
        if ~endsWith(dv, '.mat'); dv = dv + ".mat"; end
        if ~isfile(dv)
            if opts.verbose
                fprintf("%s not found!\n", opts.trait(n))
            end
            notfoundIdx(n) = true;
            continue
        end
        D = load(dv); 
        D = D.UKB_STRUCT_ALL;
    end
    
    trait_desc(n ,1) = string(D.tag); % for variable description
    if opts.legacy || isstruct(opts.trait(n))
        trait(n, 1) = string(D.tag);
    else
        [~, trait(n, 1), ext] = fileparts(opts.trait(n));
        trait(n, 1) = trait(n, 1) + ext;
    end

    if opts.verbose
        fprintf('%d of %d-processing trait %s\n', n, numel(opts.trait), trait(n))
    end
    
    f_exeid = [];
    if ~isfield(D, 'exeid')
        D.exeid = [];
    else
        f_exeid = ismember(qc_eid, D.exeid);
        if opts.verbose
            fprintf('\t%d samples were removed due to exeid constraint.\n',...
                sum(f_exeid))
        end
    end
    
    if D.numericFlag || strcmp(D.tag, 'Sex')
        % for sex 0 is female and 1 is male. Note that get_covar also codes
        % Sex in the same way. So in descriptive stat, n(%) shows men %
        [f1, f2] = ismember(qc_eid, D.eid); f2(f2 < 1) = [];
        DV_des = double(D.rawUKB(f2));
        if ~strcmpi(D.tag, 'sex')
            if opts.verbose
                fprintf('\ttransformation function: %s\n', opts.transform)
            end
            if ~strcmp(opts.transform, 'none')
                DV_des = feval(opts.transform, DV_des);
            end
        end
        
        phenotypes(f1, n) = DV_des;
    else
        phenotypes(ismember(qc_eid, D.eid), n) = 1;
        phenotypes(~ismember(qc_eid, D.eid), n) = 0;
    end
    
    bt(n, 1) = ~D.numericFlag;
    phenotypes(f_exeid, n) = NaN; 
    if opts.verbose; fprintf('\n'); end
end

phenotypes(:, notfoundIdx) = [];
trait(trait == "") = [];

try
    phenotypes = array2table([qc_eid, phenotypes], 'VariableNames', ...
        ["eid"; trait]);
catch
    phenotypes = array2table([qc_eid, phenotypes], 'VariableNames', ...
        ["eid"; matlab.lang.makeValidName(trait)]);
end

% @14MAY2023: get and merge genotype data
if all(isfield(opts, ["snp", "chr"]))
    parflag = false;
    if ~isempty(gcp("nocreate")), parflag = true; end
    geno = getbulkgeno(opts.snp, opts.chr, home=opts.genohome, ...
        model=opts.model, parallel=parflag, infoscore=2e-6, ...
        verbose=opts.verbose, merge=true, dosage=opts.dosage);
    [f1, f2] = ismember(phenotypes.eid, geno.eid); f2(f2<1) = [];
    phenotypes{:, geno.snp} = nan;
    for k = 1:numel(geno.snp)
        phenotypes.(geno.snp(k))(f1) = geno.bed(f2, k);
    end
    
    trait_desc = [trait_desc; join([geno.chr, geno.pos, geno.a1, geno.a2], ":")];
    clear geno
end

phenotypes.Properties.VariableDescriptions = ["eid"; trait_desc];
phenotypes.Properties.UserData = bt; % binary traits
if ~isempty(covariates) && all(phenotypes.eid ~= covariates.eid)
    error('getPhenoCov: QC eids are not consistent between pheno/cov!')
end

if opts.removenan && numel(trait) == 1 % remove samples with missing phenotype
    fnan = isnan(phenotypes.(2));
    phenotypes(fnan, :) = [];
    if ~isempty(covariates); covariates(fnan, :) = []; end
end

end % END 

%% subfunctions
function covar_matrix = getcovar(get_covars, eid, opts)
f_names = fieldnames(get_covars);
covar_matrix = nan(numel(eid), numel(f_names));

exeid = cell(numel(f_names), 1); % exclusion criteria for covariates
for i = 1:numel(f_names)
    if isfield(get_covars.(f_names{i}), 'exeid')
        exeid{i} = get_covars.(f_names{i}).exeid;
    end
    [f1, f2] = ismember(eid, get_covars.(f_names{i}).eid);
    f2(f2<1) = [];
    if get_covars.(f_names{i}).numericFlag
        covar_matrix(f1, i) = double(get_covars.(f_names{i}).rawUKB(f2));
    else
        if ~strcmpi(get_covars.(f_names{i}).tag, 'sex') % ignore categories
            covar_matrix(:, i) = zeros(numel(eid), 1);
            covar_matrix(f1, i) = 1; 
            continue
        end
        u_char = unique(get_covars.(f_names{i}).rawUKB(f2));
        if numel(u_char) < 2
            error('covariates cannot contain < 2 categories!')
        end
        covar_temp = zeros(numel(f2), 1);
        for j = 1:numel(u_char)
            covar_temp(get_covars.(f_names{i}).rawUKB(f2) == u_char(j)) = j-1;
        end
        covar_matrix(f1, i) = covar_temp;
    end
end

exeid(cellfun(@isempty, exeid)) = [];
if ~isempty(exeid)
    exeid = unique(vertcat(exeid{:}));
    exeidx = ismember(eid, exeid);
    if any(exeidx)
        covar_matrix(exeidx, :) = nan;
        if opts.verbose
            fprintf('%d samples will be remove due to exclusion criteria.\n', sum(exeidx))
        end
    end
end

if ~opts.defaultcovar % don't need to include additional covariates
    return
end
% add extra covariates
ext_covar = load(fullfile(opts.phenopath,'UKB_CONFOUNDERS_NUMERIC.mat'));
ext_covar = ext_covar.UKB_CONFOUNDERS;


[f1, f2] = ismember(eid, ext_covar.pcaeid); f2(f2<1) = [];
covar_matrix(f1, numel(f_names)+1:numel(f_names)+10) = ext_covar.pcadata(f2, :);

if ~opts.wes % ignore this for WES
    [f1, f2] = ismember(eid, ext_covar.binarrayeid); f2(f2<1) = [];
    covar_matrix(f1, size(covar_matrix, 2) + 1) = ext_covar.binarraydata(f2);
end

if ~isempty(exeid) % not needed though 
    covar_matrix(exeidx, :) = nan;
end

end % END