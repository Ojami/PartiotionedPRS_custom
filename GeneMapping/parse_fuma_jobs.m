function parse_fuma_jobs
% instruction: 
% download the results from FUMA webserver and extract them in folder
% "results"
clc
pth = fullfile(pwd, "results", "publish");
if ~isfolder(pth), mkdir(pth); end

parsetype = "snp2gene"; %"gene2func"; % FUMA functions SNP2GENE or GENE2FUNC

% only protein-coding genes for enrichr
g38 = struct2table(load("path_to_gtf\GTF\Homo_sapiens.GRCh38.107.gtf.gene.mat"));

if parsetype == "snp2gene"
    % parse SNP2GENE ----------------------------------------------------------
    fis = getfilenames(fullfile(pwd, "results", "SNP2GENE"), dir=true, fullpath=true);

    for k = 1:numel(fis)
        [~, name] = fileparts(fis(k)); % file name
        outpath = fullfile(pth, name);
        if ~isfolder(outpath), mkdir(outpath); end
    
        xls = fullfile(outpath, name + ".xlsx");
    
        % read MAGMA gene-based
        [info, magma] = readFUMAtab(fullfile(fis(k), "magma.genes.out"), xls);
    
        % read MAGMA gene-set 
        info = readFUMAtab(fullfile(fis(k), "magma.gsa.out"), xls, info);
    
        % read MAGMA gene-property
        info = readFUMAtab(fullfile(fis(k), ...
            "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"), xls, info);
    
        info = readFUMAtab(fullfile(fis(k), ...
            "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"), xls, info);
    
        writetable(struct2table(info), xls, Sheet="Legend")
    
        % read genes table (positional/eQTL/chromatin interaction mapping)
        genes1 = readtable(fullfile(fis(k), "genes.txt"), ...
            VariableNamingRule="preserve", TextType="string", Delimiter="\t");
        try
            genes2 = readtable(fullfile(fis(k), "1000G", "genes.txt"), ...
                VariableNamingRule="preserve", TextType="string", Delimiter="\t");
            idx = ismember(genes2.ensg, genes1.ensg) & ismember(genes2.IndSigSNPs, genes1.IndSigSNPs);
            genes = [genes1; genes2(~idx, :)];
            skip_enrichr = true;
        catch
            genes = genes1;
            skip_enrichr = true;
        end
        writetable(genes, xls, Sheet="genes")

        if skip_enrichr, continue; end
    
        % perform enrichr based on positional/eQTL/chromatin interaction
        % mapping. This is done because with chromatin interaction mapping,
        % many enriched highly significant terms will belong to chromatin
        % folding processes and not the biology of the underlying phenotype.
        % Here, I use enrichr and skipp output of FUMA GENE2FUNC module.
        % Note: MAGMA_Drugs_and_Diseases ref: https://www.medrxiv.org/content/10.1101/2022.09.06.22279660v1.full-text
        libs = ["GTEx_Tissue_Expression_Down", "GTEx_Tissue_Expression_Up",...
            "GTEx_Tissues_V8_2023", "Descartes_Cell_Types_and_Tissue_2021",...
            "MAGMA_Drugs_and_Diseases", "MSigDB_Hallmark_2020", "Reactome_2022",...
            "KEGG_2021_Human", "GO_Biological_Process_2021"];
        enrichr_path = fullfile(outpath, "enrichr");
        genes.ci = genes.ciMap.lower == "yes";
        genes.eqtl = genes.eqtlMapSNPs > 0;
        genes.pos = genes.posMapSNPs > 0;
        
        % IDs/symbols based on  GRCh38
        magma_38 = map2grch38(magma, g38);
        genes_38 = map2grch38(genes, g38);

        % MAGMA only
        magma_genes = unique(magma.Gene);
        enrichr(magma_genes, db=libs, output=fullfile(enrichr_path, "magma", name))
        gset.magma = unique(magma_38.("Gene ID"));
    
        % intersection of all
        geneset = unique(genes.symbol(genes.ci & genes.eqtl & genes.pos));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "pos_eqtl_ci_intersect", name))
        gset.pos_eqtl_ci_intersect = unique(genes_38.ensg(genes_38.ci & genes_38.eqtl & genes_38.pos));
    
        % union of all
        geneset = unique(genes.symbol(genes.ci | genes.eqtl | genes.pos));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "pos_eqtl_ci_union", name))
        gset.pos_eqtl_ci_union = unique(genes_38.ensg(genes_38.ci | genes_38.eqtl | genes_38.pos));
        
        % pos and (eqtl or ci)
        idx = (genes.pos & genes.eqtl) | (genes.pos & genes.ci);
        geneset = unique(genes.symbol(idx));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "pos_AND_eqtlORci", name))
        gset.pos_AND_eqtlORci = unique(genes_38.ensg(idx));
    
        % pos and eqtl
        geneset = unique(genes_38.symbol(genes_38.eqtl & genes_38.pos));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "posANDeqtl", name))
        gset.posANDeqtl = unique(genes_38.ensg(genes_38.eqtl & genes_38.pos));
        
        % pos or eqtl
        geneset = unique(genes_38.symbol(genes_38.eqtl | genes_38.pos));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "posOReqtl", name))
        gset.posOReqtl = unique(genes_38.ensg(genes_38.eqtl | genes_38.pos));
    
        % union of magma and (pos and eqtl)
        geneset = union(geneset, magma_genes);
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "posANDeqtl_magma", name))
        gset.posANDeqtl_magma = union(gset.posANDeqtl, gset.magma);
        
        % only positional
        geneset = unique(genes.symbol(genes.pos));
        enrichr(geneset, db=libs, output=fullfile(enrichr_path, "pos", name))
        gset.pos = unique(genes_38.ensg(genes_38.pos));
    
        % write gene sets to a file to be uploaded on FUMA GENE2FUNC
        xls = regexprep(xls, ".xlsx$", ".geneset.xlsx");
        fisets = string(fieldnames(gset));
        for j = 1:numel(fisets)
            writematrix(gset.(fisets(j)), xls, Sheet=fisets(j))
        end
    
    end

else

    % parse GENE2FUNC ---------------------------------------------------------
    fis = getfilenames(fullfile(pwd, "results", "GENE2FUNC"), dir=true, fullpath=true);
    for k = 1:numel(fis)
        [~, name] = fileparts(fis(k)); % file name
        outpath = fullfile(pth, name);
        if ~isfolder(outpath), mkdir(outpath); end
    
        xls = fullfile(outpath, name + ".xlsx");
    
        % gtex_v8_ts_DEG.txt
        readG2Ftab(fullfile(fis(k), "gtex_v8_ts_DEG.txt"), xls);
    
        readG2Ftab(fullfile(fis(k), "gtex_v8_ts_general_DEG.txt"), xls);
        
        % gene set ORA
        readG2Ftab(fullfile(fis(k), "GS.txt"), xls);
    end
end

end % END

%% subfunctions ===========================================================
function [info, out] = readFUMAtab(in, xls, info)

if nargin < 3, info = struct; end

out = readtable(in, FileType="text",...
        VariableNamingRule="preserve", TextType="string", ...
        CommentStyle="#");
out = sortrows(out, "P");

if any(colnames(out).lower == "symbol") % gene-based MAGMA
    out = out(:, ["SYMBOL", "GENE", "CHR", "START", "STOP", ...
        "NSNPS", "P"]);
    out = renamevars(out, ["SYMBOL", "GENE", "START", "STOP", "NSNPS"], ...
        ["Gene", "Gene ID", "Start", "Stop", "N SNPs"]);
end

% gene-set description
if endsWith(in.lower, "magma.gsa.out")
    sheet = "MAGMA.gene-set";
    type = "MAGMA gene-set";
    setname = "Gene sets are from Molecular Signatures Database (MsigDB) v7.0";
    out = out(:, ["FULL_NAME", "NGENES", "BETA", "BETA_STD", "SE", "P"]);
    out = renamevars(out, ["FULL_NAME", "NGENES", "BETA", "BETA_STD"], ...
        ["Gene set", "N genes in set", "Beta", "Beta (Std)"]);
elseif all(~cellfun(@isempty, regexp(in.lower, ["gtex", "magma"])))
    type = "MAGMA gene-property";
    if any(colnames(out) == "FULL_NAME")
        fcol = "FULL_NAME";
    else
        fcol = "VARIABLE";
    end
    out = out(:, [fcol, "BETA", "BETA_STD", "SE", "P"]);
    out = renamevars(out, [fcol, "BETA", "BETA_STD"], ...
        ["Tissue", "Beta", "Beta (Std)"]);
    out.Tissue = replace(out.Tissue, "_", " ");

    if contains(in.lower, "general")
        sheet = "MAGMA.gene-property.general";
        setname = "Tissue specificity analysis of " + height(out) + " GTEx (v8) general tissue types";
    else
        sheet = "MAGMA.gene-property";
        setname = "Tissue specificity analysis of " + height(out) + " GTEx (v8) tissue types";
    end
else
    setname = missing;
    sheet = "MAGMA.gene-based";
    type = "MAGMA gene-based";
end

% apply Bonferroni correction
p_cutoff = (0.05/height(out));
if nargin < 3
    info.name(1) = type;
    info.value(1) = "Bonferroni cutoff = " + p_cutoff;
    info.set(1) = setname;
    info.size(1) = height(out);
else
    info.name = [info.name; type];
    info.value = [info.value; "Bonferroni cutoff = " + p_cutoff];
    info.set = [info.set; setname];
    info.size = [info.size; height(out)];
end
out(out.P > p_cutoff, :) = [];

if ~isempty(out)
    writetable(out, xls, Sheet=sheet)
end

% check header for meta data
fid = fopen(in, "r");
fl = string;
fl(1) = string(fgetl(fid));
if startsWith(fl, "#")
    readMore = true;
    ct = 2;
    while readMore
        fl(ct, 1) = string(fgetl(fid));
        if ~fl(ct).startsWith("#")
            fl(ct) = [];
            readMore = false;
        end
        ct = ct + 1;
    end
else
    fl = string;
end
fclose(fid);

if ~any(fl == "")
    fl = extractAfter(fl, "TOTAL_GENES = ");
    fl(ismissing(fl)) = [];
else
    fl = missing;
end

if nargin < 3
    info.gene_size(1) = fl;
else
    info.gene_size = [info.gene_size; fl];
end

end

%% ------------------------------------------------------------------------
function readG2Ftab(in, xls)

out = readtable(in, FileType="text",...
        VariableNamingRule="preserve", TextType="string", ...
        CommentStyle="#");

% for tissue specifcity analyses using DEG sets, adjusted P-value os
% Bonferroni, but for gene set ORA adjustment is FDR BH
out(out.adjP > 0.05, :) = [];

if isempty(out), return; end

if in.endsWith("GS.txt")
    xls = regexprep(xls, ".xlsx$", ".GS.xlsx");
else % tissue specificity of DEGs from GTEx tissue types
    xls = regexprep(xls, ".xlsx$", ".DEG.xlsx");
end

% write each category to a different sheet
if in.endsWith("GS.txt")
    writetable(out, xls, Sheet="All", ...
        WriteMode="overwritesheet")
end
sheets = unique(out.Category);
for k = 1:numel(sheets)
    idx = out.Category == sheets(k);

    thissheet = sheets(k);
    if in.lower.contains("general") % general tissue types
        thissheet = "general." + thissheet;
    end

    tmp = out(idx, :);
    tmp.Category = [];
    writetable(tmp, xls, Sheet=truncateStr(thissheet, 31), ...
        WriteMode="overwritesheet")
end

end

%% ------------------------------------------------------------------------
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

end

end