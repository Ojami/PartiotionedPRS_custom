function plotpPRSperformance
% clc
% 28AUG2023: plot output of pPRSperformance4SLD for UKBB 

files = getfilenames(pwd, "xlsx").xlsx;
drawClustedFLD(files(files.startsWith("UKBB.")), "UKBB")

end % END

%% subfunctions ===========================================================
function drawClustedFLD(files, term)

wd = fullfile(pwd, "effPlots", term);
if ~isfolder(wd), mkdir(wd); end
viz = struct;

for k = 1:numel(files)
    name = files(k).extractBetween(".", ".");
    name = split(name, "_");
    trait = name(1);
    if trait.lower ~= "pdff"
        % trait = "Liver iron corrected T1";
        continue % @29AUG2023: skip cT1 scores
    end
    
    if name(2).lower == "conc"
        sc = "Concordant";
    elseif name(2).lower == "opp"
        sc = "Discordant";
    else
        sc = "Neutral";
    end
    
    tab = readtable(files(k), TextType="string", VariableNamingRule="preserve");
    tab(tab.univariate | tab.PRS ~= "full", :) = [];
    adi = unique(tab.adj);
    tab.outcome = tab.outcome.eraseBetween(",", textBoundary("end"));
    tab.outcome = erase(tab.outcome, ",");
    
    for j = 1:numel(adi)
        tab2 = tab(tab.adj == adi(j), :);
        tab2.trait(:) = trait;
        tab2.sc(:) = sc;
        viz.(adi(j)){k, 1} = tab2;
    end

end

fis = string(fieldnames(viz));
for k = 1:numel(fis)
    tab = viz.(fis(k));
    tab(cellfun(@isempty, tab)) = [];
    tab = vertcat(tab{:});
    tab.outcome = replace(tab.outcome, "_", " ");

    % @13NOV2023: remove neutral PRS
    tab(tab.sc == "Neutral", :) = [];

     yOrder = ["Chronic liver disease", "Cirrhosis", ...
        "HCC", "Autoimmune liver diseases", "All cancer (excluding blood)", ...
        "Cardiovascular disease (CVD)", "Heart Failure (HF)", "Diabetes",...
        "Hypertension", "Chronic kidney failure"];
    
     tab.BETA = exp(tab.BETA);
     tab.CI = [];
     effPlotter(tab, "binary", true(height(tab), 1), "save", true, ...
        "groupCol", "sc", "yCol", "outcome", "hide", true, "squeeze", false,...
        "output", fullfile(wd, term + "." + fis(k)), "fontname", "Garamond", ...
        "fontsize", 21, "sortEffect", false, "yOrder", yOrder, ...
        "legLocation", "southeast", "legTitle", "PDFF-circulating TGs PRS", ...
        "colormap", "parula", "resolution", 600, "markersize", 100)
end

end % END