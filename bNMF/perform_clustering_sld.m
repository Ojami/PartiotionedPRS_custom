function perform_clustering_sld
% @18MARCH2024: clusters previously found + novel loci based on their
% associations with continiuous cardiometabolic and liver traits to
% identify different clusters. This is a parallele attempt to see how
% hypothesis-free unsupervised clustering contrass to our first
% hypothesis-deriven discordant/concordant approach.


% traits for clustering: at baseline after exclusion of individuals with
% PDFF
traits = ["ALT", "GlycatedHaemoglobin_HbA1c_", ...
    "LDL", "Triglycerides", "AST", "Glucose", "Creatinine", ...
    "SystolicBloodPressure_AutomatedReading", "CystatinC"];

% previously known + replicated novel loci, that associate with PDFF
% (nominally). As of March 25th 2024, 6 novel replicated + 29 previously
% known compose 35 loci, with 3 overlapping (GPAM, PNPLA3 and TM6SF2),
% resulting in a final 32 loci (initial step in calculating PDFF-TG pPRS).
% Out of these 32, 6 don't associate with PDFF in UKBB. So, a total of
% 32-6=26 loci were considered for bNMF clustering.
pth = fullfile(fileparts2(pwd, 2), "PRS\culsteredPRS");
pfi = getfilenames(pth, "mat").mat;
pfi(~pfi.startsWith("PRS.UKBB.PDFF")) = [];
snps = cell(numel(pfi), 1);
for k = 1:numel(snps)
    fi = extractBetween(pfi(k), "PRS.UKBB.", ".mat");
    tmp = load(fullfile(pth, pfi(k))).prs.(fi).mixed.beta;
    snps{k} = tmp(:, ["Trait", "CHR", "POS", "SNP", "Locus", "A1", "A2", "A1Freq", "BETA"]);
end
snps = vertcat(snps{:});
idx = ~contains(snps.Trait, "(");
snps.Trait(idx) = snps.Trait(idx) + " (NA)";
snps.adipo = extractBetween(snps.Trait, "(", ")");
snps.id = snps.SNP.replace("_", ":") + " " + snps.Locus;

% align variants with negative BETA to positive BETA (increasing PDFF)
idx = snps.BETA < 0;
snps.BETA(idx) = -snps.BETA(idx);
snps.A1Freq(idx) = 1 - snps.A1Freq(idx);
a1 = snps.A1(idx);
snps.A1(idx) = snps.A2(idx);
snps.A2(idx) = a1;

% fetch association z-scores
od = fullfile(pwd, "raw");
if ~isfile(fullfile(od, "merged_phewide.mat"))
    fetchPhewideGss(traits, snps, od);
end
gss = load(fullfile(od, "merged_phewide.mat")).tab;

% perform bNMF clustering
fi = string(fieldnames(gss));
wd = fullfile(pwd, "results");
k_vec = 7;
n_reps = 1e3;
pval_cutoff = 0.05;

if ~isfile(fullfile(wd, "res_full.mat"))
    res = struct;
    for k = 1:numel(fi)
        tab = gss.(fi(k));
        tmp = array2table(tab.Z);
        tmp.Properties.VariableNames = tab.Pheno;
        tmp.Row = tab.SNP;

        tmp.Properties.UserData = tab.N;
        wd2 = fullfile(wd, fi(k) + "_full");
        
        res.(fi(k)).("k" + k_vec) = run_bNMF(tmp, out=fullfile(wd2, "maxK_" + k_vec), K=k_vec, K0=k_vec, n_reps=n_reps, pval_cutoff=pval_cutoff);
    end
    save(fullfile(wd, "res_full.mat"), "res")
end


end % END