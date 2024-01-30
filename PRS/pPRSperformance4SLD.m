function pPRSperformance4SLD
clc
% @21AUG2023: checks the performance and association of the partitioned PDFF
% and cT1 PRS in predicting CLD, cirrhosis and cardiometabolic traits
% after excluding those with available PDFF or cT1.

% disease end-points
traits.ukbb.base = ["HCC", "ChronicLiverDisease", "Cirrhosis",...
    "CardiovascularDisease_CVD_", "Diabetes", "Hypertension",...
    "AutoimmuneLiverDiseases", "AllCancer_excludingBlood_", ...
    "ChronicKidneyFailure", "HeartFailure_HF_", ...
    "Hyperthyroidism", "Dyslipidemia"];
traits.ukbb.mri = "ProtonDensityFatFraction_PDFF_"; % sanity check

% load partitioned PRS files
prs_files = getfilenames(fileparts(pwd), "mat", "fullpath", true).mat;
prs_files(~prs_files.contains(filesep + "PRS.")) = [];

% step 1: test in UKBB ----------------------------------------------------
ukb_prs_files = prs_files(prs_files.contains(filesep + "PRS.UKBB."));
for k = 1:numel(ukb_prs_files)
    prs = load(ukb_prs_files(k)).prs;
    checkPRSperformance(prs, traits, "ukbb")
end

end % END

%% subfcuntions ===========================================================
function checkPRSperformance(prs, traits, pop)
prs_name = string(fieldnames(prs));
traits = traits.(pop);

if pop == "ukbb"
    traits = union(traits.mri, traits.base, "stable");
end

if isfile(upper(pop) + "." + prs_name(1) + "_continuous.xlsx"), return; end


if pop == "ukbb"
    qc = 0; % European & at most 10 putative third-degree relatives
    qc_eid = getQCEID(qc, false);
    
    % get covariates
    covars.base = ["Age", "Sex", "BMI", "Cholesterol"];%, "Hypertension", "Diabetes" "Alcohol_consumption_g_day" % sensitivity analysis
    covars.mri = ["AgeMRI", "Sex", "BMI_2", "Cholesterol"];%, "Hypertension", "Diabetes" % sensitivity analysis
    fi = string(fieldnames(covars));
    for k = 1:numel(fi)
        ctab.(fi(k)) = getPhenoCov(trait=covars.(fi(k)), covar="", defaultcovar=false, ...
            transform="none", qc=qc, verbose=false, legacy=false);
        ctab.(fi(k)) = rmmissing(ctab.(fi(k)));
        if fi(k) == "mri"
            ctab.(fi(k)) = renamevars(ctab.(fi(k)), "AgeMRI", "Age");
        end
        ctab.(fi(k)).Age_2 = ctab.(fi(k)).Age.^2;
        ctab.(fi(k)).Sex_Age = ctab.(fi(k)).Sex.*ctab.(fi(k)).Age;
        ctab.(fi(k)).Sex_Age_2 = ctab.(fi(k)).Sex.*ctab.(fi(k)).Age_2;
    end
end

tt = tic;
res = struct; ct = 1;
for k1 = 1:numel(prs_name)
    
    if pop == "ukbb"
        % to check the predictive performance of PRS, we test them in samples
        % withouth them (independent internal samples): only for UKBB
        if prs_name(k1).startsWith("PDFF")
            rm_eid = load("ProtonDensityFatFraction_PDFF_.mat").UKB_STRUCT_ALL.eid;
        else
            rm_eid = load("LiverIronCorrectedT1_ct1_.mat").UKB_STRUCT_ALL.eid;
        end
    end
    
    prs2 = prs.(prs_name(k1));
    % different flavors of PRS based on adiposity adjustment
    prs_adj = string(fieldnames(prs2));
    for k2 = 1:numel(prs_adj)
        
        prs3 = prs2.(prs_adj(k2)).prs;
        prs_type = ["full", "novel", "old"]; 
        for k3 = 1:numel(prs_type)
            
            try
                prs4 = prs3(:, ["eid", prs_type(k3)]);
            catch
                continue
            end
            
            for k4 = 1:numel(traits)

                sz = numel(prs_name)*numel(prs_adj)*numel(prs_type)*numel(traits);
                fprintf("%d of %d\n", ct, sz)
                
                if pop == "ukbb"    
                    outcome = load(traits(k4)).UKB_STRUCT_ALL;
                end


                data = struct;
                if pop == "ukbb" && any(traits(k4) == mri_triats) % do not test for out of sample subset in case of MRI traits
                    data.eid = qc_eid;
                else
                    data.eid = setdiff(qc_eid, rm_eid); % remove samples used for calculating PRS weights
                end
                

                [data.(prs_type(k3)), data.(traits(k4))] = deal(nan(numel(data.eid), 1));
                [f1, f2] = ismember(data.eid, prs4.eid);
                data.(prs_type(k3))(f1) = prs4.(prs_type(k3))(f2(f1));

                if outcome.numericFlag % cont traits are dichotomized based on median for AUROC
                    [idx, idx2] = ismember(data.eid, outcome.eid);
                    data.(traits(k4))(idx) = double(outcome.rawUKB(idx2(idx)));

                else % binary trait
                    idx = ismember(data.eid, outcome.eid);
                    data.(traits(k4))(idx) = 1;
                    data.(traits(k4))(~idx) = 0;
                end
    
                if isfield(outcome, "exeid") % exclusion criteria
                    idx = ismember(data.eid, outcome.exeid);
                    data.(traits(k4))(idx) = nan;
                end

                data = struct2table(data);
                data = rmmissing(data);

                % remove outcome from covariates
                col1 = colnames(data); col1(col1 == "eid") = [];
                
                if pop == "ukbb"
                    if any(traits(k4) == mri_triats)
                        col2 = colnames(ctab.mri);
                        data = outerjoin(data, ctab.mri(:, setdiff(col2, col1)), "Keys", "eid", "MergeKeys", true);
                    else
                        col2 = colnames(ctab.base);
                        data = outerjoin(data, ctab.base(:, setdiff(col2, col1)), "Keys", "eid", "MergeKeys", true);
                    end
                end
                data.eid = [];
                data = movevars(data, prs_type(k3), Before=1);
                data = movevars(data, traits(k4), After=width(data));
                
                % fit univariate and multiple regression models for
                % continuous PRS

                if outcome.numericFlag
                    data.(traits(k4)) = irnt(data.(traits(k4)));
                    mdl1 = fitlm(data(:, [prs_type(k3), traits(k4)]));
                    mdl2 = fitlm(data);
                    [stmp, dtmp] = getGLMstat(mdl1, mdl2, data(:, [prs_type(k3), traits(k4)]));

                    med = median(data.(traits(k4)), "omitmissing");
                    data.(traits(k4))(data.(traits(k4)) > med) = 1;
                    data.(traits(k4))(data.(traits(k4)) <= med) = 0;
                    
                    mdl1 = fitglm(data(:, [prs_type(k3), traits(k4)]), ...
                        Link="logit", Distribution="binomial");
                    mdl2 = fitglm(data, Link="logit", Distribution="binomial");

                else
                    mdl1 = fitglm(data(:, [prs_type(k3), traits(k4)]), ...
                        Link="logit", Distribution="binomial");
                    mdl2 = fitglm(data, Link="logit", Distribution="binomial");
    
                    [stmp, dtmp] = getGLMstat(mdl1, mdl2, data(:, [prs_type(k3), traits(k4)]));
                end
                
                stmp.PRS(:) = prs_type(k3);
                [stmp.outcome(:), dtmp.outcome(:)] = deal(string(outcome.tag));
                [stmp.adj(:), dtmp.adj(:)] = deal(prs_adj(k2));
                [stmp.score(:), dtmp.score(:)] = deal(prs_name(k1));
                
                % plot ROC
                name = prs_name(k1) + "." + prs_adj(k2) + "." + prs_type(k3) + "." + traits(k4);
                auc_res = getROC(mdl1, mdl2, data(:, [prs_type(k3), traits(k4)]), name, pop);
                stmp.AUC = [auc_res.AUC; auc_res.AUC_adj];

                res.stat{ct, 1} = stmp;
                res.desc{ct, 1} = dtmp;

                ct = ct + 1;

            end % over different diseases

        end % prs type: full, novel or old set of variants
        
    end % prs adipo adjustments (e.g. NA, BMI, ...)

end % prs name

toc(tt)

% write the results
stats = vertcat(res.stat{:, 1});
desc = vertcat(res.desc{:, 1});
stats.AUC = double(stats.AUC);
sc = unique(stats.score);
for k = 1:numel(sc)
    stats_in = stats(stats.score == sc(k), :);
    stats_in.score = [];
    writetable(stats_in, upper(pop) + "." + sc(k) + "_continuous.xlsx", Sheet="stat");

    dec_in = desc(desc.score == sc(k), :);
    dec_in.score = [];
    writetable(dec_in, upper(pop) + "." + sc(k) + "_continuous.xlsx", Sheet="desc");
end


end % END

%% subfunctions ===========================================================
function auc_res = getROC(mdl1, mdl2, data, name, pop)

name = split(name, ".");
wd = fullfile(pwd, "ROCfigs." + upper(pop), name(1), name(2));
if ~isfolder(wd), mkdir(wd); end
name = fullfile(wd, name(3) + "." + name(4));

nboot = 1000;
score = colnames(data, index=1);
data.(2) = categorical(data.(2));
roc1 = rocmetrics(data.(2), mdl1.Fitted.Probability, 1, ...
    NumBootstraps=nboot, ...
    BootstrapOptions=statset(UseParallel=true), ...
    BootstrapType="corrected percentile");

roc2 = rocmetrics(data.(2), mdl2.Fitted.Probability, 1, ...
    NumBootstraps=nboot, ...
    BootstrapOptions=statset(UseParallel=true), ...
    BootstrapType="corrected percentile");

auc_res = struct;
if nboot == 1000
    auc_res.AUC(1) = compose("%.3f", roc1.AUC');
    auc_res.AUC_adj(1) = compose("%.3f", roc2.AUC');
    showCI = false;
else
    auc_res.AUC(1) = compose("%.3f (%.3f-%.3f)", roc1.AUC');
    auc_res.AUC_adj(1) = compose("%.3f (%.3f-%.3f)", roc2.AUC');
    showCI = true;
end

if isfile(name + ".fig"), return, end
close all force
ax1 = plot(roc1, ShowModelOperatingPoint=false, ShowConfidenceIntervals=showCI);
ax1.DisplayName = regexprep(ax1.DisplayName, "^1", score);

hold on
ax2 = plot(roc2, ShowModelOperatingPoint=false, ShowConfidenceIntervals=showCI);
ax2.DisplayName = regexprep(ax2.DisplayName, "^1", score + " (adjusted)");
ax2.Parent.Legend.Location = "southeast";
ax2.Parent.FontName = "Garamond";
ax2.Parent.FontSize = 12;
ax1.Parent.Title.String = [];

savefig(gcf, name + ".fig")
exportgraphics(gcf, name + ".png", "Resolution", 300)

end
%% ------------------------------------------------------------------------

function [tab, res] = getGLMstat(mdl1, mdl2, data)

mdls = {mdl1, mdl2};
tab = struct;
for k = 1:numel(mdls)
    ci = mdls{k}.coefCI;
    ci = ci(2, :);
    ci = "[" + round(ci(1), 2, "decimals") + ", " + round(ci(2), 2, "decimals") + "]";
    tab.BETA(k, 1) = mdls{k}.Coefficients.Estimate(2);
    tab.SE(k, 1) = mdls{k}.Coefficients.SE(2);
    tab.P(k, 1) = mdls{k}.Coefficients.pValue(2);
    tab.CI(k, 1) = ci;
    tab.N(k, 1) = mdls{k}.NumObservations;
end
tab = struct2table(tab);
tab.univariate(:) = false;
tab.univariate(1) = true;

res = innerDescriptive(data, colnames(data, index=width(data)));

end

%% ------------------------------------------------------------------------
function res = innerDescriptive(pheno, outcome)

cols = setdiff(colnames(pheno), ["eid", outcome]);
pheno(ismissing(pheno.(outcome)), :) = [];
idx = pheno.(outcome) == 0;
 
for k = 1:numel(cols)
    tmp = pheno.(cols(k));
    idx_tmp = idx;
    idx_missing = ismissing(tmp);
    idx_tmp(idx_missing) = []; tmp(idx_missing) = [];
    
    res.PRS(k, 1) = cols(k);
    res.N_overall(k, 1) = height(pheno);
    res.N_missing(k, 1) = sum(idx_missing);
    res.N_n(k, 1) = compose("%d (%.2f %%)", sum(idx_tmp), 100*sum(idx_tmp)/numel(idx_tmp));
    res.N_c(k, 1) = compose("%d (%.2f %%)", sum(~idx_tmp), 100*sum(~idx_tmp)/numel(idx_tmp));
    
    if numel(unique(tmp)) > 2
        res.meansd_n(k, 1) = compose("%.2f±%.2f", mean(tmp(idx_tmp)), std(tmp(idx_tmp)));
        res.meansd_c(k, 1) = compose("%.2f±%.2f", mean(tmp(~idx_tmp)), std(tmp(~idx_tmp)));
    
        res.medianiqr_n(k, 1) = compose("%.2f (%.2f)", median(tmp(idx_tmp)), iqr(tmp(idx_tmp)));
        res.medianiqr_c(k, 1) = compose("%.2f (%.2f)", median(tmp(~idx_tmp)), iqr(tmp(~idx_tmp)));

    else % binary predictor
        tmp_n = tmp(idx_tmp);
        res.N_n_n(k, 1) = compose("%d (%.2f %%)", sum(tmp_n == 0), 100*sum(tmp_n)/numel(tmp_n));
        res.N_n_c(k, 1) = compose("%d (%.2f %%)", sum(tmp_n), 100*sum(tmp_n)/numel(tmp_n));

        tmp_n = tmp(~idx_tmp);
        res.N_c_n(k, 1) = compose("%d (%.2f %%)", sum(tmp_n == 0), 100*sum(tmp_n)/numel(tmp_n));
        res.N_c_c(k, 1) = compose("%d (%.2f %%)", sum(tmp_n), 100*sum(tmp_n)/numel(tmp_n));
    end
    
end

res = struct2table(res);

end
