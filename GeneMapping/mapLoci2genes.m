function mapLoci2genes

% NOTE: should be run after call_regenie_gene_based
% @05MAY2023: gene-based results should be run over novel loci identified
% via different adiposity adjustments. These loci usually point to
% different nearby genes w/wo supporting eQTL/CI/MAGAM gene-based evidence.
% Therefore, this function unifies the list of genes returned from FUMA
% gene-mapping over different adiposity adjustments. This list can be
% further narrowed down to include only specific criteria, such as
% eQTL/MAGMA/positionally mapped genes.
%
% NOTE: we keep FUMA mapped genes in an adipo-specific manner. This means
% from each set of mapped genes by FUMA, only mapped genes for loci in each
% adipo-specific adjustments are kept. For instance, if there are 4 SNPs in
% the final table with strongest p-value (by final table we mean after
% clumping and keeping the strongest variant per each loci along with its
% adipo adjustment) that have been adjusted for BMI, then only mapped genes
% for those 4 variants will be kept from the FUMA results (i.e. FUMA PDFF
% or cT1 adjusted for BMI). This better reflects the adipo adjusted
% approach of this project, and further reduces the computational time
% (because otherwise there will be overlapping genes for the same locus
% from all FUMA files).

%@30AUG2023: for novel loci, we keep only replicated loci on external
%cohorts and discard the rest.

% @01SEP2023: we add fine-mapped genes (nearest gene to variant with the
% highest PIP at each locus), colocalized genes and genes with highest V2G
% score from Open Targets Genetics to the list of FUMA genes.

% @04SEP2023: for all independent loci (overrides @30AUG2023). To be used
% for tables/figures before replication using all loci identified in UKBB.

mergeFlag = true; % merge PDFF and cT1 genes.
intersectOnly = false; % to test only intersection of genes over adjustments in 'adiAdj'

if intersectOnly
    res_tag = "intersect";
else
    res_tag = "union";
end

if mergeFlag
    res_tag = res_tag + "_merged";
end

wd = fullfile(fileparts(pwd), "FUMA", "results", "publish");
fuma_folders = getfilenames(wd, "", "dir", true);
gene_folders = getfilenames(pwd, "", "dir", true);

% only these adjustmetns (from penalized phenotypic predictive models).
% adiAdj is used also to pick the overlapping genes. For instance, if gene
% A appears in both BMI and WFM, the summary stat in BMI is picked, since
% it is the first index in 'adiAdj' vector. This order is based on the
% sample size (N BMI > WFM > VAV) and adjustment (NA is the last)
adiAdj = ["BMI", "WFMm2", "VAV", "NA"]; 
fuma_folders(~endsWith(fuma_folders, adiAdj)) = [];
gene_folders(~endsWith(gene_folders, replace(adiAdj, "WFMm2", "WFM_m2"))) = [];

% -------------------------------------------------------------------------
%@04SEP2023: using all identified loci in UKBB (all loci in Table 1 of
%manuscript)
gss.pdff = readtable(fullfile(fileparts(pwd), "tables", "publish", ...
    "PDFF.xlsx"), TextType="string", ...
    VariableNamingRule="preserve", Sheet="publish.minorAllele");

gss.ct1 = readtable(fullfile(fileparts(pwd), "tables", "publish", ...
    "LiverIronCorrectedT1.xlsx"), TextType="string", ...
    VariableNamingRule="preserve", Sheet="publish.minorAllele");

% filter variants in LD (decide based on the strongest association in
% each locus)
gss.pdff = groupfilter(gss.pdff, "Locus ID", @(x) x == min(x), "P");
gss.ct1 = groupfilter(gss.ct1, "Locus ID", @(x) x == min(x), "P");

% add (NA) tag to unadjusted traits, and adipo tag
idx = gss.pdff.Trait.contains("(");
gss.pdff.Trait(~idx) = gss.pdff.Trait(~idx) + " (NA)";
gss.pdff.adipo = extractBetween(gss.pdff.Trait, "(", ")");
idx = gss.ct1.Trait.contains("(");
gss.ct1.Trait(~idx) = gss.ct1.Trait(~idx) + " (NA)";
gss.ct1.adipo = extractBetween(gss.ct1.Trait, "(", ")");

% load fine-mapped genes 
fmtmp = load(fullfile(fileparts(pwd), "Finemapping\publish\PDFF.cT1.finemap.mat")).res;
idx = fmtmp.Trait.contains("PDFF");
fmt = struct;
fmt.pdff = fmtmp(idx, :);
fmt.ct1 = fmtmp(~idx, :);
clear fmtmp

% for each locus keep the genes with largest PIP
fi = ["pdff", "ct1"];
fm = struct;
for j = 1:numel(fi)
    fmtmp = fmt.(fi(j));
    leadVars = unique(fmtmp.Index);
    fmtab = cell(numel(leadVars), 1);
    for k = 1:numel(leadVars)
        idx = fmtmp.Index == leadVars(k);
        tmp = fmtmp(idx, :);
        tmp = tmp(tmp.PIP == max(tmp.PIP), :);
        fmtab{k} = tmp;
    end
    fm.(fi(j)) = vertcat(fmtab{:});
    fm.(fi(j)) = groupsummary(fm.(fi(j)), ["Index", "Locus", "Trait"], ...
        @(x)join(unique(x), "|"), "nearestGene");
    fm.(fi(j)) = renamevars(fm.(fi(j)), "fun1_nearestGene", "nearestGene");
end


% load highest V2G genes for each locus from Open Targets Genetics
v2g.pdff = load(fullfile(fileparts(pwd), "OpenTargetsGenetics\PDFF.V2G.mat")).tab;
v2g.ct1 = load(fullfile(fileparts(pwd), "OpenTargetsGenetics\cT1.V2G.mat")).tab;
v2g.pdff.Symbol(v2g.pdff.Symbol == "-") = missing;
v2g.ct1.Symbol(v2g.ct1.Symbol == "-") = missing;


fi = ["pdff", "ct1"];
for k = 1:numel(fi)
    
    gss.(fi(k)) = renamevars(gss.(fi(k)), "coloc genes", "coloc");

    % add fine-mapped genes (highest PIP at each locus)
    [~, idx] = ismember(gss.(fi(k)).SNP, fm.(fi(k)).Index);
    gss.(fi(k)).finemap = fm.(fi(k)).nearestGene(idx);

    % add V2G
    [~, idx] = ismember(gss.(fi(k)).SNP, v2g.(fi(k)).SNP);
    gss.(fi(k)).v2g = v2g.(fi(k)).Symbol(idx);
end

% -------------------------------------------------------------------------
% load ENSEMBL GTF files
g38 = struct2table(load("Homo_sapiens.GRCh38.107.gtf.gene.mat"));
g37 = struct2table(load("Homo_sapiens.GRCh37.87.gtf.gene.mat")); % to exclude non-novel loci's nearest genes from _novel gene sets

for k = 1:numel(fuma_folders)
    if startsWith(fuma_folders(k).lower, "pdff")
        term = "protondensity";
        pheno = "PDFF (" + extractAfter(fuma_folders(k), "_") + ")";
        phenoFi = "pdff"; 
    elseif startsWith(fuma_folders(k).lower, "ct1")
        term = "liveriron";
        pheno = "Liver iron corrected T1 (" + extractAfter(fuma_folders(k), "_") + ")";
        phenoFi = "ct1"; 
    end
    
    file = fullfile(wd, fuma_folders(k), fuma_folders(k) + ".xlsx");
    genes = readtable(file, Sheet="genes", TextType="string", VariableNamingRule="preserve");
    magma = readtable(file, Sheet="MAGMA.gene-based", TextType="string", VariableNamingRule="preserve");
    genes.ci = genes.ciMap.lower == "yes";
    genes.eqtl = genes.eqtlMapSNPs > 0;
    genes.pos = genes.posMapSNPs > 0;

    % remove chromatin interaction pairs
    s = arrayfun(@(x)split(x, ";"), unique(genes.IndSigSNPs), uni=false);
    s = unique(vertcat(s{:}));
    s = s(s.contains(":") & ~s.contains("_"));
    genes.IndSigSNPs = regexprep(genes.IndSigSNPs, "(|[;])" + s + "(|[;])", "");

    % expand 'genes' table so that each row contains exactly only one locus
    idx = genes.IndSigSNPs.contains(";");
    if any(idx)
        tgn = genes(idx, :);
        genes(idx, :) = [];
        genes2 = cell(height(tgn), 1);
        for j = 1:height(tgn)
            tmpgenes = tgn(j, :);
            tmpindsig = unique(split(tmpgenes.IndSigSNPs, ";"));
            tmpgenes = repmat(tmpgenes, numel(tmpindsig), 1);
            tmpgenes.IndSigSNPs = tmpindsig;
            genes2{j} = tmpgenes;
        end
        genes = [genes; vertcat(genes2{:})];
    end
    genes.snps = genes.IndSigSNPs;
    genes = unique(genes, "rows");

    % sanity check: all gss SNPs should be present in FUMA SNPs. Moreover
    % (see second note in the beginning of this function), only mapped
    % genes of loci for current adipo adjustment will be only kept.
    fuma_adipo = extractAfter(fuma_folders(k), "_");
    subgss = gss.(phenoFi)(gss.(phenoFi).adipo == fuma_adipo, :); % for mapped genes to corresponding loci connection
    checkSNP = subgss.SNP;

    if isempty(checkSNP)
        % i.e. this FUMA adipo adjusted table has no overlap with the final
        % loci (e.g. if in this loop FUMA_NA is being processed, but there
        % is no NA adjusted loci in the gss struct)
        continue
    else
        assert(all(ismember(checkSNP, genes.snps))) % sanity check
        keepIdx = ismember(genes.snps, checkSNP);
        genes(~keepIdx, :) = [];

    end

    % remove empty 'snps': only CI pairs. Note that this doesn't
    % necessarily mean we remove all CI pairs (e.g. rs73221948 for PDFF VAT
    % has only CI evidence but is not empty).
    genes(genes.snps == "", :) = [];


    % IDs/symbols should be based on GRCh38
    genes = map2grch38(genes, g38);
    magma = unique(magma.("Gene ID")); % parse_fuma_jobs function has already mapped them to GRCh38

    % @31AUG2023/@01SEP2023 add other evidence: add nearest gene (using
    % "locus" column in the final clumped tables), V2G, finemapped, and
    % colocalized genes to the "genes" table if FUMA didn't add them to
    % positional mapped gene set due to window constraint (currently: 50
    % kb).
    genes = addOtherEvidence2FUMA(genes, subgss, g37, g38);

    idx_novel = gss.(phenoFi).GC_novel.lower == "true";
    matchingSNPs = gss.(phenoFi).SNP(idx_novel);
    genes.novel = ismember(genes.snps, matchingSNPs);

    % map loci's nearest genes (for previousely known loci) to GRCh38
    gss_loci_g38 = g37.gene_id(ismember(g37.gene_name, gss.(phenoFi).Locus(~idx_novel)));
    if ~isempty(setdiff(gss_loci_g38, g38.gene_id))
        disp("warning: some genes from common variant analysis (GWAS) cannot be mapped from v37->38!")
    end
    
    % @30AUG2023
    % since we don't want to lose the connection between each mapped gene
    % and it's related loci (SNP + nearest gene), we tag gene names by
    % their corrsponding loci.
    [~, idx] = ismember(genes.snps, subgss.SNP);
    genes.locusTag = subgss.SNP(idx) + " " + subgss.Locus(idx);
    genes.Trait = subgss.Trait(idx);

    
    % different set of genes
    gset.pheno(k, 1) = pheno;
    gset.term(k, 1) = term;
    gset.locusTag{k, 1} = genes(:, ["ensg", "locusTag", "Trait", "pos", "eqtl", "ci", "coloc", "finemap", "v2g"]); 
    gset.allevidence{k, 1} = unique(genes.ensg(genes.eqtl | genes.ci | genes.pos | genes.coloc | genes.finemap | genes.v2g)); % pos,ci, eqtl,coloc,finemap or v2g
    gset.magma{k, 1} = magma;

end

% get the union of gene sets for each trait and its adiposity adjustments
all_genes = union(vertcat(gset.allevidence{:}), vertcat(gset.allevidence_magma{:})); % complete set of genes to be tested

% create a table to pairt each locus to it's mapped genes. Some genes maye
% share more than one locus
loctab = vertcat(gset.locusTag{:});
loctab = renamevars(loctab, ["ensg", "locusTag"], ["Gene", "Locus"]);
loctab.magma(:) = false;

% add MAGMA genes: add another column for MAGMA evidence. But "onlyMagma"
% genes don't have any specific locus
onlyMagma.Gene = setdiff(all_genes, loctab.Gene);
loctab(~ismember(loctab.Gene, all_genes), :) = [];
onlyMagma = struct2table(onlyMagma);
om = cell(height(gset.pheno), 1);
for k = 1:numel(gset.pheno)
    idx = ismember(onlyMagma.Gene, gset.magma{k});
    tmp = onlyMagma(idx, :);
    tmp.Trait(:) = gset.pheno(k);
    tmp.magma(:) = true;
    om{k} = tmp;

    idx = ismember(loctab.Gene, gset.magma{k});
    loctab.magma(idx) = true;
end
onlyMagma = vertcat(om{:});
onlyMagma.Locus(:) = "MAGMA";
[onlyMagma.pos(:), onlyMagma.eqtl(:), onlyMagma.ci(:), onlyMagma.coloc(:), onlyMagma.finemap(:), onlyMagma.v2g(:)] = deal("undefined");
loctab = [loctab; onlyMagma];

loctab.Symbol(:) = "";
[f1, f2] = ismember(loctab.Gene, g38.gene_id);
loctab.Symbol(f1) = g38.gene_name(f2(f1));
idx = loctab.Symbol == "";
if any(idx) % no gene name in GRCh38, try with GRCh37
    loctab2 = loctab(idx, :);
    loctab(idx, :) = [];
    [f1, f2] = ismember(loctab2.Gene, g37.gene_id);
    loctab2.Symbol(f1) = g37.gene_name(f2(f1));
    loctab = [loctab; loctab2];
end
loctab = movevars(loctab, "Symbol", After="Gene");
cols = ["pos", "eqtl", "ci", "coloc", "finemap", "v2g"];
for k = 1:numel(cols)
    thisscore = loctab.(cols(k));
    thisscore = replace(thisscore, ["true", "false", "undefined"], ["1", "0", "0"]);
    loctab.(cols(k)) = logical(double(thisscore));
end
loctab = renamevars(loctab, [cols, "magma"], ["Pos-mapped (FUMA)", ...
    "eQTL (FUMA)", "Chromatin interaction (FUMA)", "Colocalization",...
    "Fine-mapped", "V2G", "MAGMA"]);
loctab.Rank = sum(loctab{:, 5:end}, 2); % based on sum of all evidence

% remove duplicates: FUMA was run twice: with UKBB and 1000 G, so there
% maybe duplicates. Whenever there is duplicate, keep the one with highest
% rank
loctab = unique(loctab, "rows");
loctab.id = loctab.Gene + ":" + loctab.Locus + ":" + loctab.Trait;
if ~isempty(duplicates(loctab.id))
    loctab = groupfilter(loctab, "id", @(x) x == max(x), "Rank");
end
loctab.id = [];

% sort based on rank at each locus
loctab.id = loctab.Locus + ":" + loctab.Trait;
uid = unique(loctab.id);
ltab = cell(numel(uid), 1);
for k = 1:numel(uid)
    idx = loctab.id == uid(k);
    tmp = loctab(idx, :);
    loctab(idx, :) = [];
    ltab{k} = sortrows(tmp, "Rank", "descend");
end
loctab = vertcat(ltab{:});
loctab.Trait = replace(loctab.Trait, ["(WFMm2)", "(VAV)"], ["(WFM)", "(VAT)"]);
loctab.id = [];
save("locus2gene.all.mat", "loctab")

idx = loctab.Trait.startsWith("PDFF");
writetable(loctab(idx, :), "locus2gene.PDFF.xlsx")
writetable(loctab(~idx, :), "locus2gene.cT1.xlsx")

end % END

%% subfunctions ===========================================================
function gs = map2grch38(gs, g38)

if any(colnames(gs) == "Gene") % MAGMA
    sy = "Gene";
    id = "Gene ID";
else
    sy = "symbol";
    id = "ensg";
end

[f1, f2] = ismember(gs.(id), g38.gene_id); f2(f2<1) = [];
gs.(sy)(f1) = g38.gene_name(f2);

if any(~f1)
    f1 = find(~f1);
    gs_tmp = gs(f1, :);
    [g1, g2] = ismember(gs_tmp.(sy), g38.gene_name); g2(g2<1) = [];
    gs.(id)(f1(g1)) = g38.gene_id(g2);

    % check genes that cannot be mapped to GRCh38. Use ENSEMBL xref API:
    gs_tmp = gs_tmp(~g1, :);
    f1 = f1(~g1);
    if ~isempty(gs_tmp)
        fprintf("\t%d genes cannot be mapped to GRCh38, trying Ensembl API\n", height(gs_tmp))
        gs_tmp.g38(:) = "";
        for k = 1:height(gs_tmp)
            res = EnsemblREST(gs_tmp.(sy)(k), "xref", "refGenome", "38");
            if ~isempty(res)
                gs_tmp.g38(k) = intersect(res.id, g38.gene_id);
            end
        end

        if any(gs_tmp.g38 == "")
            fprintf("\t\tstill could not map %d genes\n", sum(gs_tmp.g38 == ""))
        end

        [g1, g2] = ismember(gs_tmp.g38, g38.gene_id); g2(g2<1) = [];
        gs.(id)(f1(g1)) = g38.gene_id(g2);
        gs.(sy)(f1(g1)) = g38.gene_name(g2);
    end
end

end

%% ------------------------------------------------------------------------
function fuma_all = addOtherEvidence2FUMA(fuma, subgss, g37, g38)
% appends and adds other genes found from V2G, colocalization, fine-mapping
% and nearest genes to FUMA 'genes' table.

[fuma.coloc(:), fuma.finemap(:), fuma.v2g(:)] = deal(false);
cc = ["Locus", "coloc", "finemap", "v2g"];
fuma_all = fuma; % original fuma table

for k1 = 1:height(subgss)
    idx = fuma_all.snps == subgss.SNP(k1);
    fuma = fuma_all(idx, :);
    fuma_all(idx, :) = [];

    % loop over each evidence column
    for k2 = 1:numel(cc)
        g1 = struct;
        g1.name = subgss.(cc(k2))(k1);
        g1.name = g1.name.split([",", ";", "|"]);
        g1.name(g1.name == "" | g1.name == "-" | ismissing(g1.name)) = [];
        if isrow(g1.name), g1.name = g1.name'; end
        if any(ismissing(g1.name)) || any(ismember(g1.name, ["", "-"]))
            continue
        end

        % now get gene id for gene names in 'g' struct
        g1.ensg = strings(numel(g1.name), 1);
        g1 = struct2table(g1);
        [f1, f2] = ismember(g1.name, g38.gene_name);
        g1.ensg(f1) = g38.gene_id(f2(f1));
        idx = g1.ensg == "";
        if any(idx) % missing gene id in GRCh38, try with GRCh37
            g2 = g1(idx, :);
            g1(idx, :) = [];
            [f1, f2] = ismember(g2.name, g37.gene_name);
            g2.ensg(f1) = g37.gene_id(f2(f1));
            g = [g1; g2];
        else
            g = g1;
        end
        g.snp(:) = subgss.SNP(k1);
        
        % now check if this gene is in FUMA table for this locus
        [f1, f2] = ismember(g.ensg, fuma.ensg);
        f1m = find(f1); f2m = f2(f1);
        f1nm = find(~f1); %f2nm = f2(~f1);

        if cc(k2).lower == "locus"
            col = "pos"; % equivalent to positionally mapped genes
        else
            col = cc(k2);
        end

        % for those already existing in FUMA table add their evidence
        for j = 1:numel(f1m)
            fuma.(col)(f2m(j)) = true;
            % tmpSNPs = fuma.snps{(f2m(j))};
            % tmpSNPs = union(subgss.SNP(k1), tmpSNPs);
            % if isrow(tmpSNPs), tmpSNPs = tmpSNPs'; end
            % fuma.snps{(f2m(j))} = tmpSNPs;
        end

        % for those missing from FUMA table add them to FUMA table
        if ~isempty(f1nm)
            g2 = g(f1nm, :);
            fuma2 = fuma(1:height(g2), :);
            fuma2{:, ["eqtl", "pos", "v2g", "finemap", "coloc", "ci"]} = false;
            fuma2.(col)(:) = true;
            fuma2.ensg = g2.ensg;
            fuma2.symbol = g2.name;
            fuma2.chr(:) = double(subgss.CHR(k1));
            % fuma2.snps(:) = {subgss.SNP(k1)};
            fuma = [fuma; fuma2];
        end

    end

    fuma_all = [fuma_all; fuma]; % add corrected 'fuma' table to 'fuma_all' table

end

end