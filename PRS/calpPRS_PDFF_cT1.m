function calpPRS_PDFF_cT1
% clc

% @19MAY2023: calculate PRS using effect sizes estimated in UKBB using all
% and previously known variants. This serves two purposes: 
% 1- To compare R2 of old + new vs old variants
% 2- To use these PRS for DE analysis in MAFALDA (comparing 10th vs 90th
% percentiles).

% @18AUG2023: we define partitioned PRS by grouping variants based on their
% associatin with cT1/PDFF and circulating TGs and divide them into two
% groups:
%       1- those with concordant association directions between PDFF/cT1 and TG
%       2- those with opposite directions between PDFF/cT1 and TG
%       3- those associate only with PDFF/cT1 but not with TG


if isempty(gcp("nocreate"))
    parpool("Processes", 40);
end

% create a dir for mat files
if ~isfolder("matfiles"), mkdir("matfiles"); end

% create a dir for quantile plots of PRS
if ~isfolder("figs"), mkdir("figs"); end

% for quantile plots
rawPheno = getPhenoCov(trait=["LiverIronCorrectedT1_ct1_", "ProtonDensityFatFraction_PDFF_"], ...
    covar="", defaultcovar=false, transform="none", qc=0, verbose=false);
rawPheno = renamevars(rawPheno, ...
    ["Liver iron corrected T1 (ct1)", "Proton density fat fraction (PDFF)"], ...
    ["ct1", "pdff"]);

ztab = load("PDFF_cT1_novel.mat").tab2;
ztab.Trait = ztab.SNP.extractAfter(",").strip;
ztab.SNP2 = ztab.SNP.extractBefore(",");
ztab.SNP2 = strip(erase(ztab.SNP2, ztab.Locus));

idx.ns = ztab.Pheno == "PDFF" & ztab.p > 0.05; % not associated with PDFF
idx.s = ztab.Pheno == "PDFF" & ztab.p <= 0.05; % associated with PDFF
ztab.beta(ztab.p > 0.05) = nan;
fi = string(fieldnames(idx));
for k = 1:numel(fi)
    csnps = unique(ztab.SNP(idx.(fi(k))));
    if fi(k) == "s"
        trt = "PDFF";
    else
        trt = "Liver iron corrected T1"; % use cT1 for non-associated loci with PDFF
    end

    fidx = ismember(ztab.SNP, csnps) & ismember(ztab.Pheno, ["Triglycerides", trt]);
    tmp = ztab(fidx, :);
    prsgrp.(fi(k)) = struct;
    [prsgrp.(fi(k)).conc, prsgrp.(fi(k)).opp, prsgrp.(fi(k)).ns] = deal(string);
    for j = 1:numel(csnps)
        tmp2 = tmp(tmp.SNP == csnps(j), :);
        idx1 = tmp2.Pheno == "Triglycerides";
        idx2 = tmp2.Pheno == trt;
        if tmp2.beta(idx1)*tmp2.beta(idx2) > 0
            prsgrp.(fi(k)).conc = union(prsgrp.(fi(k)).conc, csnps(j));
        elseif tmp2.beta(idx1)*tmp2.beta(idx2) < 0
            prsgrp.(fi(k)).opp = union(prsgrp.(fi(k)).opp, csnps(j));
        else
            prsgrp.(fi(k)).ns = union(prsgrp.(fi(k)).ns, csnps(j));
        end
    end
    
    % create tables with adipo adjustments to each
    % group (conc,opp,ns)-pheno (cT1, PDFF) and calculate PRS
    sc = string(fieldnames(prsgrp.(fi(k))));
    for j = 1:numel(sc)
        tmpsnps = prsgrp.(fi(k)).(sc(j));
        tmpsnps(tmpsnps == "") = [];
        if isempty(tmpsnps)
            prsgrp.(fi(k)).(sc(j)) = tmpsnps;
            continue
        end
        fidx = tmp.Pheno == trt & ismember(tmp.SNP, tmpsnps);
        prsgrp.(fi(k)).(sc(j)) = tmp(fidx, :);

        % calculate PRS
        calPRS(tmp(fidx, :), rawPheno, trt + "_" + sc(j))
    end
end


end % END

%% subfunctions ===========================================================
function calPRS(tmp, rawPheno, pheno)

pheno = matlab.lang.makeValidName(pheno);

 % fetch GWAS summary stats for each adjustment-snp pair
if isfile(fullfile("matfiles", pheno + ".gwas.mat"))
    gtab = load(fullfile("matfiles", pheno + ".gwas.mat")).gtab;
else
    gtab = fetchGSS(tmp, false);
    save(fullfile("matfiles", pheno + ".gwas.mat"), "gtab")
end

% calculate different set of PRS:
% 1- mixed: using 'tmp' table, where each locus is strongest
%    association over different adiposity adjustments.
% 2- adiposity specific: using 'gtab' cell of tables, where each cell
%    contains a table using different adjustments.

% get the genotype matrix
if isfile(fullfile("matfiles", pheno + ".geno.mat"))
    geno = load(fullfile("matfiles", pheno + ".geno.mat")).geno;
else
    geno = getbulkgeno(gtab{1}.SNP, string(gtab{1}.CHR), merge=true);
    save(fullfile("matfiles", pheno + ".geno.mat"), "geno")
end

% check SNPxSNP interactions (caveat: qc 0 is used!)
if pheno.lower.startsWith("liver")
    outcome = "LiverIronCorrectedT1_ct1_";
else
    outcome = "ProtonDensityFatFraction_PDFF_";
end

for j = 1:numel(gtab)
    tmp = gtab{j};
    tag = unique(tmp.Trait);
    if numel(tag) > 1
        tag = "mixed";
        cvs = "";
    elseif ~contains(tag, "(")
        tag =  "NA";
        cvs = "";
    else
        tag = extractBetween(tag, "(", ")");
        if tag == "BMI"
            cvs = "BMI_2";
        elseif tag == "VAV"
            cvs = "VisceralAdiposeTissueVolume_VAT_";
        elseif tag == "WFM"
            cvs = "WholeBodyFatMass_2";
        end
    end
    
    if isfile(fullfile("matfiles", pheno + "." + tag + ".mat"))
        ret = load(fullfile("matfiles", pheno + "." + tag + ".mat")).ret;
    else
        cvs = union(cvs, ["AgeMRI", "Sex"]); cvs(cvs == "") = [];
        ret = getPRSpheno(geno, outcome, cvs);
        save(fullfile("matfiles", pheno + "." + tag + ".mat"), "ret")
    end

    [~, idx] = ismember(geno.snp, tmp.SNP);
    tmp = tmp(idx, :);

    [~, ridge_idx] = ismember(matlab.lang.makeValidName(geno.snp), preg.beta.feature);
    preg.beta = preg.beta(ridge_idx, :);

    % note: A2 is minor allele in geno, while A1 is minor allele in 'tmp'
    assert(all(geno.a1 == tmp.A2))
    assert(all(geno.a2 == tmp.A1))

    prstmp = struct;
    prstmp.eid = geno.eid;
    prstmp.full = geno.bed*tmp.BETA; % previously known + new loci
    prstmp.full_ridge = geno.bed*preg.beta.full;% previously known + new loci using Ridge betas

    prstmp = struct2table(prstmp);
    prstmp = rmmissing(prstmp);

    % calculate incremental R2: R2 w/ PRS - R2 w/o PRS
    qfis = "qc";
    for j1 = 1:numel(qfis)
        retTmp = ret.(qfis(j1));
        [f1, f2] = ismember(retTmp.eid, prstmp.eid); f2(f2<1) = [];
        retTmp(~f1, :) = []; 
        retTmp.eid = [];
        mdl1 = fitlm(retTmp);
        schemes = setdiff(colnames(prstmp), "eid");

        for j2 = 1:numel(schemes)
            mdl2 = fitlm([prstmp(f2, schemes(j2)), retTmp]);
            mdl3 = fitlm([prstmp(f2, schemes(j2)), retTmp(:, end)]);
            r2.(qfis(j1)).(schemes(j2)) = mdl3.Rsquared.Adjusted;
            r2inc.(qfis(j1)).(schemes(j2)) = mdl2.Rsquared.Adjusted - mdl1.Rsquared.Adjusted;
            pval.(qfis(j1)).(schemes(j2)) = mdl3.Coefficients(2, :);
            pval.(qfis(j1)).(schemes(j2)).N(:) = mdl3.NumObservations;
        end
    end

    % for quantile plots
    tmpPheno = rawPheno;
    if pheno.lower.startsWith("pdff")
        tmpPheno = tmpPheno(:, ["eid", "pdff"]);
        tmpPheno = rmmissing(tmpPheno);
        tmpPheno = renamevars(tmpPheno, "pdff", "pheno");
        tmpPheno.mark = tmpPheno.pheno > 5.5; % steatosis
        ylabo = "PDFF (%)";
    else
        tmpPheno.mark = tmpPheno.pdff > 5.5 & tmpPheno.ct1 > 800; % NASH
        tmpPheno = tmpPheno(:, ["eid", "ct1", "mark"]);
        tmpPheno = rmmissing(tmpPheno);
        tmpPheno = renamevars(tmpPheno, "ct1", "pheno");
        ylabo = "Liver iron corrected T1 (ms)";
    end
    [f1, f2] = ismember(tmpPheno.eid, prstmp.eid); f2(f2<1) = [];
    tmpPheno(~f1, :) = [];

    for j2 = 1:numel(schemes)
        score = prstmp.(schemes(j2))(f2);
        bins = [min(score), quantile(score, 9), max(score)];
        vizo = struct;
        vizo.idx = discretize(score, bins);
        vizo = struct2table(vizo);
        vizo(:, ["pheno", "mark"]) = tmpPheno(:, ["pheno", "mark"]);

        transFunc = ["none", "irnt", "log"];
        
        for j3 = 1:numel(transFunc)
            if ~isfolder(fullfile("figs", transFunc(j3)))
                mkdir(fullfile("figs", transFunc(j3)))
            end

            fig_name = pheno + "." + tag + "." + schemes(j2) + "." + transFunc(j3);
            fig_name = fullfile("figs", transFunc(j3), fig_name);
            if isfile(fig_name + ".fig"), continue; end

            close gcf force
            viz = vizo;
            ylab = ylabo;
            if transFunc(j3) == "irnt"
                viz.pheno = irnt(viz.pheno);
                ylab = regexprep(ylab, "\(.*?\)", "(IRNT)");
            elseif transFunc(j3) == "log"
                viz.pheno = log(viz.pheno);
                ylab = regexprep(ylab, "\(.*?\)", "(log)");
            end

            fax = figure(Position=[1, 41, 1280, 720], Units="pixels");
            viz.mark = (viz.mark + 1).*4;
            viz.idx = categorical(viz.idx);
            ax = swarmchart(viz, "idx", "pheno", "filled", "ColorVariable", "mark");
            ax.SizeData = 15;
            ax.MarkerEdgeColor = "black";
            ax.MarkerEdgeAlpha = 0.3;
            ax.MarkerFaceAlpha = 0.5;
            
            hold(ax.Parent, "on")
            h = boxchart(viz.idx, viz.pheno);
            h.MarkerStyle = "none";
            h.WhiskerLineColor = "r";
            colormap(ax.Parent, "parula")

            ax.Parent.FontName = "Garamond";
            ax.Parent.FontSize = 14;
            ax.Parent.LineWidth = 1;
            xlabel(ax.Parent, "Quantiles of PRS")
            ylabel(ax.Parent, ylab)
            
            exportgraphics(fax, fig_name + ".jpg", "Resolution", 600)
            savefig(fax, fig_name + ".fig")
        end
    end
    
    prs.(pheno).(tag).r2 = r2;
    prs.(pheno).(tag).stat = pval;
    prs.(pheno).(tag).r2inc = r2inc;
    prs.(pheno).(tag).prs = prstmp;
    prs.(pheno).(tag).beta = tmp;
    prs.(pheno).(tag).ridge = preg;
end
save("PRS.UKBB." + pheno + ".mat", "prs")

% print R2 and incremental R2 to output tables
fi = string(fieldnames(prs.(pheno)));
for j1 = 1:numel(qfis)

    [r2inctab, r2tab, stat] = deal(cell(numel(fi), 1));
    for j2 = 1:numel(fi)
        r2tab{j2} =  prs.(pheno).(fi(j2)).r2.(qfis(j1));
        r2tab{j2}.Adjustment = fi(j2);

        r2inctab{j2} =  prs.(pheno).(fi(j2)).r2inc.(qfis(j1));
        r2inctab{j2}.Adjustment = fi(j2);

        stat{j2} = prs.(pheno).(fi(j2)).stat.(qfis(j1));
        stat{j2}.Adjustment = fi(j2);
    end
    r2tab = struct2table(vertcat(r2tab{:}));
    r2inctab = struct2table(vertcat(r2inctab{:}));

    newstat = cell(numel(stat), 1);
    for m = 1:numel(stat)
        newstat{m} = splitvars(struct2table(stat{m}));
    end
    newstat = vertcat(newstat{:});

    writetable(newstat, pheno + ".goodnesOfFit.xlsx", Sheet=qfis(j1) + ".stat")
    writetable(r2tab, pheno + ".goodnesOfFit.xlsx", Sheet=qfis(j1) + ".R2")
    writetable(r2inctab, pheno + ".goodnesOfFit.xlsx", Sheet=qfis(j1) + ".incremental.R2")
end

end % END

%% subfunctions ===========================================================
function ctab = fetchGSS(snps, revFlag)

if nargin < 2
    revFlag = false; % to fetch associations with PDFF for cT1 loci and vice versa
end

% GWAS files
gfis = getfilenames(fileparts2(pwd, 3), dir=true, fullpath=true);

if any(colnames(snps) == "Pheno")
    ut = unique(snps.Pheno); 
    ut = ut + " (" + ["BMI", "WFMm2", "VAV", "NA"] + ")";
    fiFlag = true;
else
    ut = unique(snps.Trait);
    fiFlag = false;
end

ctab = cell(numel(ut), 1);
for k = 1:numel(ut)
    % fprintf("fetching summary stats: %d of %d\n", k, numel(ut))

    ut_snps = snps;
    if fiFlag
        ut_snps = ut_snps(:, ["SNP2", "CHR", "POS", "A1", "A2", "Locus", "variant"]);
        ut_snps = renamevars(ut_snps, "SNP2", "SNP");
    end
    
    adipo = extractBetween(ut(k), "(", ")");

    if fiFlag
        if adipo == "WFMm2", adipo = "WFM_m2"; end
    else
        if ut(k).contains("(")
            adipo = extractBetween(ut(k), "(", ")");
            if adipo == "WFMm2", adipo = "WFM_m2"; end
        else
            adipo = "NA";
        end
    end

    if ut(k).lower.startsWith("liver")
        if revFlag
            patt = "protondensity";  
            pheno = "PDFF";
        else
            patt = "liveriron";
            pheno = "Liver iron corrected T1";
        end
    else
        if revFlag
            patt = "liveriron";
            pheno = "Liver iron corrected T1";
        else
            patt = "protondensity";  
            pheno = "PDFF";
        end
    end

    idx = gfis.lower.contains(patt) & gfis.endsWith(adipo);
    ssfile = getfilenames(gfis(idx), "txt", fullpath=true).txt;
    ssfile(~endsWith(ssfile, ".step2.txt")) = [];

    gss = readGWASfile(ssfile, index=ut_snps, win=0, full=true,...
        light=false, parallel=true, n=true, legacy=false);
    if adipo ~= "NA", pheno = pheno + " (" + adipo + ")"; end
    gss.Trait(:) = pheno;
    gss = movevars(gss, "Trait", Before=1);
    gss.Properties.VariableNames = createGWASheader(colnames(gss));
    
    gss = renamevars(gss, ["BP", "A1FREQ"], ["POS", "A1Freq"]); % to be consistent with the tables in the publish folder
    gss.N = [];

    % align the variants' effect allele (A1)
    [~, idx] = ismember(ut_snps.SNP, gss.SNP); 
    gss = gss(idx, :);
    idx = ut_snps.A1 ~= gss.A1;
    if any(idx)
        a1 = gss.A1(idx);
        gss.A1(idx) = gss.A2(idx);
        gss.A2(idx) = a1;
        gss.A1Freq(idx) = 1 - gss.A1Freq(idx);
        gss.BETA(idx) = -gss.BETA(idx);
    end

    assert(all(ut_snps.A1 == gss.A1))
    assert(all(ut_snps.A2 == gss.A2))

    % add locus names
    gss(:, ["Locus", "variant"]) = ut_snps(:, ["Locus", "variant"]);  

    ctab{k} = gss;
end

ctab(cellfun(@isempty, ctab)) = [];

if ~fiFlag, return; end

% create aipo adjustment specific summary stat table
mtab = ctab{1};
assert(all(mtab.SNP == snps.SNP2))
adipo = ["BMI", "WFM", "VAT", "NA"];
for k = 1:height(mtab)
    idx = adipo == snps.adipo(k);
    if ~any(idx), error("wrong adipo!"); end
    mtab.BETA(k) = ctab{idx}.BETA(k);
    mtab.SE(k) = ctab{idx}.SE(k);
    mtab.A1Freq(k) = ctab{idx}.A1Freq(k);
    mtab.Trait(k) = ctab{idx}.Trait(k);
    mtab.P(k) = ctab{idx}.P(k);
end
ctab{end+1} = mtab;

end

%% ========================================================================
function ret = getPRSpheno(geno, outcome, cvs)
[c1,c2] = getPhenoCov("covar", cvs, "legacy", false, "trait", outcome);
data = join(c1, c2);
data = rmmissing(data);
data.age2 = data.AgeMRI.^2;
data.ageXSex = data.AgeMRI.*data.Sex;
data.age2XSex = data.age2.*data.Sex;
cols = setdiff(colnames(data)', ["eid", outcome]);

[f1, f2] = ismember(data.eid, geno.eid); f2(f2<1) = [];
data{f1, geno.snp} = nan;
data{f1, geno.snp} = geno.bed(f2, :);
data = rmmissing(data);
data.Properties.VariableNames = matlab.lang.makeValidName(colnames(data));

% for inceremental R2 (goodness of fit)
ret.qc = data(:, ["eid", cols', outcome]);
data.eid = [];

end