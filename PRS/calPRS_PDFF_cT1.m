function calPRS_PDFF_cT1
clc
% @19MAY2023: calculate PRS using effect sizes estimated in UKBB using all
% and previously known variants. This serves two purposes: 
% 1- To compare R2 of old + new vs old variants

if isempty(gcp("nocreate"))
    parpool("Processes", 40);
end

% get list of final independent index variants in UKBB for PDFF and cT1
% (same as Table 1 in the manuscript)
pth = fullfile(fileparts(pwd), "tables", "publish");
xls = getfilenames(pth, "xlsx", "fullpath", true).xlsx;
xls(contains(xls, "~$")) = [];

% create a dir for mat files
if ~isfolder("matfiles"), mkdir("matfiles"); end

% create a dir for quantile plots of PRS
if ~isfolder("figs"), mkdir("figs"); end

% for quantile plots
rawPheno = getPhenoCov(trait=["LiverIronCorrectedT1_ct1_", "ProtonDensityFatFraction_PDFF_"], ...
    covar="", defaultcovar=false, transform="none", qc=0);
rawPheno = renamevars(rawPheno, ...
    ["Liver iron corrected T1 (ct1)", "Proton density fat fraction (PDFF)"], ...
    ["ct1", "pdff"]);

if ~isfile("PRS.UKBB.mat")
    for k = 1:numel(xls)
        tmp = readtable(xls(k), TextType="string", ...
            VariableNamingRule="preserve", Sheet="publish.minorAllele");
        [~, pheno] = fileparts(xls(k));
    
        % filter variants in LD (decide based on the strongest association in
        % each locus)
        tmp = groupfilter(tmp, "Locus ID", @(x) x == min(x), "P");
    
        % fetch GWAS summary stats for each adjustment-snp pair
        if isfile(fullfile("matfiles", pheno + ".gwas.mat"))
            gtab = load(fullfile("matfiles", pheno + ".gwas.mat")).gtab;
        else
            gtab = fetchGSS(tmp);
            save(fullfile("matfiles", pheno + ".gwas.mat"), "gtab")
        end
    
        % calculate different set of PRS:
        % 1- mixed: using 'tmp' table, where each locus is strongest
        %    association over different adiposity adjustments.
        % 2- adiposity specific: using 'gtab' cell of tables, where each cell
        %    contains a table using different adjustments. This is done for
        %    the sake of consistency among different adjustments.
        gtab{end+1} = tmp;
    
        % get the genotype matrix
        if isfile(fullfile("matfiles", pheno + ".geno.mat"))
            geno = load(fullfile("matfiles", pheno + ".geno.mat")).geno;
        else
            geno = getbulkgeno(tmp.SNP, string(tmp.CHR), merge=true);
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
            
            if isfile(fullfile("matfiles", pheno + "." + tag + ".ridge.mat"))
                preg = load(fullfile("matfiles", pheno + "." + tag + ".ridge.mat")).preg;
                ret = load(fullfile("matfiles", pheno + "." + tag + ".ridge.mat")).ret;
            else
                cvs = union(cvs, ["AgeMRI", "Sex"]); cvs(cvs == "") = [];
                [intEffs, preg, ret] = snpXsnpInteraction(geno, outcome, cvs);
                save(fullfile("matfiles", pheno + "." + tag + ".snpXsnp.mat"), "intEffs")
                save(fullfile("matfiles", pheno + "." + tag + ".ridge.mat"), "preg", "ret")
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
        
            tmp_novel = tmp(tmp.GC_novel.lower == "true", :);
            [~, idx] = ismember(tmp_novel.SNP, geno.snp);
            prstmp.novel = geno.bed(:, idx)*tmp_novel.BETA; % novel loci
            prstmp.novel_ridge = geno.bed(:, idx)*preg.beta.full(idx);% previously known + new loci using Ridge betas
            
            prstmp.old = prstmp.full - prstmp.novel;
            prstmp.old_ridge = prstmp.full_ridge - prstmp.novel_ridge;
            prstmp = struct2table(prstmp);
            prstmp = rmmissing(prstmp);

            % calculate incremental R2: R2 w/ PRS - R2 w/o PRS
            qfis = ["qc", "qcmax"];
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
                    pval.(qfis(j1)).(schemes(j2)) = mdl3.Coefficients(2, :);
                    pval.(qfis(j1)).(schemes(j2)).N(:) = mdl3.NumObservations;
                    r2inc.(qfis(j1)).(schemes(j2)) = mdl2.Rsquared.Adjusted - mdl1.Rsquared.Adjusted;
                    
                end
            end

            % for quantile plots
            tmpPheno = rawPheno;
            if pheno.lower == "pdff"
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
        
        % print R2 and incremental R2 to output tables
        fi = string(fieldnames(prs.(pheno)));
        for j1 = 1:numel(qfis)

            [r2inctab, r2tab, stat] = deal(cell(numel(fi), 1));
            for j2 = 1:numel(fi)
                r2tab{j2} = prs.(pheno).(fi(j2)).r2.(qfis(j1));
                r2tab{j2}.Adjustment = fi(j2);
    
                r2inctab{j2} = prs.(pheno).(fi(j2)).r2inc.(qfis(j1));
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

    end
    
    save("PRS.UKBB.mat", "prs")
end

end % END

%% subfunctions ===========================================================
function ctab = fetchGSS(snps)

% GWAS files
gfis = getfilenames(fileparts(fileparts(pwd)), dir=true, fullpath=true);

ut = unique(snps.Trait);
ctab = cell(numel(ut), 1);
for k = 1:numel(ut)
    % fprintf("fetching summary stats: %d of %d\n", k, numel(ut))

    ut_snps = snps;

    if ut(k).contains("(")
        adipo = extractBetween(ut(k), "(", ")");
        if adipo == "WFMm2", adipo = "WFM_m2"; end
    else
        adipo = "NA";
    end

    if ut(k).lower.startsWith("liver")
        patt = "liveriron";
        pheno = "Liver iron corrected T1";
    else
        patt = "protondensity";  
        pheno = "PDFF";
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
    gss(:, ["Locus", "Locus ID", "GC_novel", "variant"]) = ut_snps(:, ["Locus", "Locus ID", "GC_novel", "variant"]);  

    ctab{k} = gss;
end

ctab(cellfun(@isempty, ctab)) = [];

end

%% ========================================================================
function [interactionEffs, preg, ret] = snpXsnpInteraction(geno, outcome, cvs)
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
slist = setdiff(colnames(data), union(["eid", outcome], cols));

% ridge effects on a hold-out CV
eid = getQCEID(1, false); % maximum unrelated Europeans
idx = ismember(data.eid, eid);
[ridge_beta, ridge_stat] = trainCVmodels(data(idx, [slist, outcome]), ...
    "covariates", data(idx, cols), ...
    "kfold", 2, "targetTransform", "none", "tool", "ridge", ...
    "standardize", false); % standardize:false to compare with marginal effects from REGENIE
ridge_beta = ridge_beta{1};
ridge_beta(1, :) = [];
ridge_stat.aic(:, ["outcome", "test_N", "train_N"]) = ridge_stat.best(:, ["outcome", "test_N", "train_N"]);
ridge_stat = [ridge_stat.best; ridge_stat.aic];
preg.beta = ridge_beta;
preg.stat = ridge_stat;

% for inceremental R2 (goodness of fit)
ret.qc = data(:, ["eid", cols', outcome]);
ret.qcmax = data(idx, ["eid", cols', outcome]);
data.eid = [];

interactionEffs = ({}); ct = 1;
for j1 = 1:numel(slist)-1
    fprintf("snp1: %d of %d\n", j1, numel(slist)-1)
    for j2 = j1+1:numel(slist)
        fprintf("\tsnp2: %d of %d\n", j2, numel(slist))
        mdl = fitlm(data, outcome + "~" + slist(j1) + "*" + slist(j2) + "+" + join(cols, "+"));
        idx = ismember(mdl.Coefficients.Row, [slist(j1) + ":" + slist(j2), slist(j2) + ":" + slist(j1)]);
        interactionEffs{ct, 1} = mdl.Coefficients(idx, :);
        ct = ct + 1;
    end
end

interactionEffs = vertcat(interactionEffs{:});

end