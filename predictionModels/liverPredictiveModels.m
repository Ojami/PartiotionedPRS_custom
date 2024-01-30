function liverPredictiveModels
% regularized linear regression models with adiposity indices as predictors
clc

traits = ["ProtonDensityFatFraction_PDFF_", "LiverIronCorrectedT1_ct1_"];
traits_tag = ["PDFF", "Liver iron corrected T1"];
preds = ["BMI_2", "ImpedanceOfWholeBody_MRI_", "WFM_kgm2_2",...
    "waist_to_hip_ratio_MRI", "VisceralAdiposeTissueVolume_VAT_"];
preds_tag = ["BMI", "Impedence of whole body", "WFM (kg/m2)",...
    "Waist to hip ratio", "Visceral adipose tissue volume (VAT)"];

out_tag = ".VAV";

covars = ["AgeMRI", "Sex"]; % standard covariate effects to be removed
covars_tag = ["Age", "Sex"];

data = getPhenoCov(trait=unique([traits, preds, covars]), covar="", ...
     defaultcovar=false, transform="none", qc=0, verbose=false);
data = renamevars(data, ["Age at MRI", "BMI_2", ...
    "Impedance of whole body (MRI)", "Liver iron corrected T1 (ct1)", ...
    "Proton density fat fraction (PDFF)", "Waist to hip ratio (MRI)"], ...
    ["Age", "BMI", "Impedence of whole body", "Liver iron corrected T1", ...
    "PDFF", "Waist to hip ratio"]);

% correlation matrix
raw = data;
raw(:, ["eid", covars_tag]) = [];
cols = colnames(raw);
cols(cols.startsWith("Impedence")) = "IWB";
cols(cols.startsWith("Waist to hip ratio")) = "WHR";
cols(cols.startsWith("Visceral adipose tissue")) = "VAT";
ct = 1;
for k1 = 1:numel(cols)-1
    for k2 = k1:numel(cols)

        if k1 == k2
            continue
        else
            V1 = raw{:, k1};
            V2 = raw{:, k2};
            idx = ismissing(V1) | ismissing(V2);
            [rg, pval] = corr((V1(~idx)), (V2(~idx)), "type", "Spearman");
        end
        if pval < realmin, pval = realmin; end
        rmat.p1(ct, 1) = cols(k1);
        rmat.p2(ct, 1) = cols(k2);
        rmat.rg(ct, 1) = rg;
        rmat.p(ct, 1) = pval;
        ct = ct + 1;
    end
end
rmat = struct2table(rmat);
gcorrplot(rmat, method="corr", order="hclust", ...
    out="correlation_pearson" + out_tag, type="lower", ...
    height=8, addCoef=true, pdf=true)


outvis = cell(numel(traits_tag), 1);
for k = 1:numel(traits_tag)
    tab = data(:, [preds_tag, traits_tag(k)]);
    Z = data(:, covars_tag);
    Z.age2 = Z.Age.^2;
    Z.ageSex = Z.Age.*Z.Sex;
    Z.age2Sex = Z.age2.*Z.Sex;

    idx = any(ismissing(tab), 2) | any(ismissing(Z), 2);
    tab(idx, :) = [];
    Z(idx, :) = [];

    [beta, stat] = trainCVmodels(tab, kfold=2, covariates=Z,...
            targetTransform="irnt", ignore_aic=false, ...
            tool=["elnet", "lasso", "ridge"], shap=false, standardize=true); ...
    
    beta = beta{1};
    file = traits_tag(k) + out_tag + ".xlsx";

    writetable(beta, file, Sheet="beta." + stat.best.method(1))
    writetable([stat.mse; stat.best], file, Sheet="stat." + stat.best.method(1))
    
    beta(1, :) = [];
    beta.Pheno(:) = traits_tag(k);
    beta.beta = beta.full;
    outvis{k, 1} = beta;

end

outvis = vertcat(outvis{:});
output = "standardized_beta_" + stat.best.method(1) + out_tag;
effPlotter(outvis, output=output, ...
    resolution=300, format="jpg", ...
    groupCol="Pheno", yCol="feature", markersize=100, ybold=false, ...
    savefig=true, save=true, fontsize=16, ...
    addTitle=true, sortEffect=false, fontname="Garamond", merge=true, ...
    ignoreCI=true, legendOrientation="vertical", ...
    legendCols=1, ciLineWidth=0.6, color=["BlueViolet", "OrangeRed"])

end % END