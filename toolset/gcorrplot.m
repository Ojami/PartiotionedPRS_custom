function gcorrplot(res, opts)
% correlation plot for bivariate genetic correlation in 'res', calculated
% from ldsc.
% 
% a wrapper for:
%   - https://github.com/mkanai/ldsc-corrplot-rg/blob/master/plot_corrplot_rg.R
%   - https://github.com/mkanai/ldsc-network-plot/blob/master/plot_network_rg.R
% 
% Oveis Jamialahmadi, University of Gothenburg, February 2023.

arguments
    res {mustBeA(res, "table")} % table from LDSC. Example: "gcorrplot_example.mat"
    opts.out {mustBeTextScalar} = "gcorrPlot" % name of output image
    opts.pdf (1,1) logical = false
    opts.res (1,1) double = 300 % image resolution
    opts.width (1,1) double = 8 % image width in inch
    opts.height (1,1) double = 7 % image height in inch
    opts.fontname {mustBeTextScalar} = "Garamond"

    opts.map (1,1) {mustBeMember(opts.map, ["jet", "turbo", "parula", "hot", "winter", "copper", "bone"])} = "turbo"
    opts.method {mustBeMember(opts.method, ["corr", "network"])} = "network"

    % if method == "corr"
    opts.type {mustBeMember(opts.type, ["upper", "lower", "full"])} = "full"
    opts.order {mustBeMember(opts.order, ["original", "AOE", "FPC", "hclust", "alphabet"])} = "original"
    opts.addCoef (1,1) logical = false % to add correlation coefficient
    opts.adj (1,1) logical = true % to adjust for BH
    opts.title {mustBeTextScalar} = ""

end

if opts.addCoef
    opts.addCoef = "'black'";
    opts.sig = "n";
else
    opts.addCoef = "NULL";
    opts.sig = "pch";
end

fontname = opts.fontname;
opts = rmfield(opts, 'fontname');

wd = fileparts(opts.out);
if isempty(wd) || wd == ""
    wd = pwd;
    opts.out = fullfile(wd, opts.out);
end

opts.out = regexprep(opts.out, ".png$", "");
opts.out = opts.out + "_" + opts.method + ".png";

opts.out = replace(opts.out, filesep, "/"); % R path

% filtering and adjusted P-values
res(isnan(res.p), :) = [];
if opts.adj
    res.q = mafdr(res.p, "BHFDR", true);
else
    res.q = res.p;
end
cols = colnames(res);

% category columns for each trait pair
catcols = ["p1", "p2"] + "_category";
for k = 1:numel(catcols)
    if ~any(cols == catcols(k))
        res.(catcols(k)) = res.(regexprep(catcols(k), "_category$", ""));
        res.(catcols(k)) = erase(res.(catcols(k)), "_" + alphanumericsPattern + textBoundary("end"));
    end
end

% create color table
phenonames = union(res.p1, res.p2);
traits.TRAIT = phenonames;
mfis = ["p1", "p2"];
traits.CATEGORY =  traits.TRAIT;
for k = 1:numel(mfis)
    [f1, f2] = ismember(traits.TRAIT, res.(mfis(k))); f2(f2<1) = [];
    traits.CATEGORY(f1) = res.(mfis(k) + "_category")(f2);
end

catHex = rgb2hex(feval(opts.map, numel(unique(traits.CATEGORY))));
[~, idx] = ismember(traits.CATEGORY, unique(traits.CATEGORY, "stable"));
traits.COLOR = catHex(idx);
traits = struct2table(traits);

writetable(res, "rg_gcorrplot.txt", Delimiter="\t")
writetable(traits, "traits_gcorrplot.txt", Delimiter="\t")

r(1) = "rga = read.table('rg_gcorrplot.txt', T, sep='\t', as.is = T)";
r(2) = "trait_all = read.table('traits_gcorrplot.txt', T, sep = '\t', as.is = T, quote = '', fileEncoding='utf-8', comment.char='')";

if opts.method == "network"
    r(3) = "rg_sig = rga %>% filter(q < 0.05)";
    r(4) = "rel_rg = data.frame(from = rg_sig$p1, to = rg_sig$p2, " + ...
        "directed = F, weight = abs(rg_sig$rg), q = rg_sig$q, type = 'rg')";
    r(5) = "g = graph.data.frame(rel_rg, directed=F)";
    r(7) = "########### layout ############";
    r(8) = "l = layout_with_fr(g)";
    r(9) = "l = norm_coords(l)";
    r(10) = "V(g)$x = l[,1]";
    r(11) = "V(g)$y = l[,2]";

    r(13) = "########### vertex visual ###########";
    r(14) = "V(g)$color = trait_all$COLOR[match(V(g)$name, trait_all$TRAIT)]";
    r(15) = "V(g)$shape = 'circle'";
    r(16) = "V(g)$label = V(g)$name";
    r(17) = "V(g)$label.color = 'black'";
    r(18) = "V(g)$label.font = 2";
    r(19) = "V(g)$label.family = 'sans'";
    r(20) = "V(g)$label.cex = 0.7";
    r(21) = "V(g)$label.degree = -pi/2";
    r(22) = "V(g)$label.dist = 1";

    r(23) = "# duplicate lines";
    r(24) = "tmp = rga";
    r(25) = "tmp$p1 = rga$p2";
    r(26) = "tmp$p2 = rga$p1";
    r(27) = "tmp$p1_category = rga$p2_category";
    r(28) = "tmp$p2_category = rga$p1_category";
    r(29) = "rg_dup = rbind(rga, tmp)";
    r(30) = "rg_sig_n = rg_dup %>% filter(q < 0.05) %>% group_by(p1) %>% count()";
    r(31) = "n_v_cls = diff(range(rg_sig_n$n))+1";
    r(32) = "v_cls_rg = cut2(rg_sig_n$n[match(V(g)$name, rg_sig_n$p1)], g = n_v_cls)";
    r(33) = "v_size_rg = seq(6, 20, length.out = n_v_cls)[match(v_cls_rg, levels(v_cls_rg))]";
    r(34) = "V(g)$size = v_size_rg";

    r(36) = "########### edge visual###########";
    r(37) = ['rg_cols = colorRampPalette(c("#67001F", "#B2182B",', ...
        ' "#D6604D", "#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0", ', ...
        '"#92C5DE", "#4393C3", "#2166AC", "#053061"))(200)'];
    r(38) = "edge_cls_rg = cut2(rg_sig$rg, cuts = seq(-1, 1, length.out = 201))";
    r(39) = "edge_color_rg = rg_cols[match(edge_cls_rg, levels(edge_cls_rg))]";
    r(40) = "edge_width = -log10(rg_sig$q)";    
    r(41) = "edge_width[edge_width > 10] = 10";
    r(42) = "E(g)$color = edge_color_rg";
    r(43) = "E(g)$width = edge_width";
    r(44) = "E(g)$curved = .2";
    r(45) = "E(g)$arrow.size = .5";
    r(46) = "trait_all = trait_all %>% distinct(CATEGORY, .keep_all = T)";
    
    r(49) = "png('" + opts.out + "', width = " + opts.width + ...
        ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
    
    r(51) = "plot(g, rescale =F)"; % layout=layout_in_circle(g)
    r(52) = "legend('topright', bty='n',legend=trait_all$CATEGORY, " + ...
        "fill=trait_all$COLOR, inset=c(0,1), xpd=T, cex=0.8, title='Trait')";
    r(53) = "dev.off()";
    
    if opts.pdf
        r(55) = "grDevices::cairo_pdf('" + regexprep(opts.out, ".png$", ".pdf") + "', width = " + opts.width + ...
            ", height = " + opts.height + ")";
        r(56) = "plot(g, rescale =F)";
        r(57) = "legend('topright', bty='n',legend=trait_all$CATEGORY, " + ...
            "fill=trait_all$COLOR, inset=c(0,1), xpd=T, cex=0.8, title='Trait')";
        r(58) = "dev.off()";
    end

    rtmp(1) = "pkg = c('dplyr', 'Hmisc', 'igraph', 'RColorBrewer', 'stringr')";
    rtmp(2) = "lapply(pkg, require, character.only = TRUE)";
    rtmp(3) = "set.seed(123)";
    r = [r(1:2), rtmp, r(3:end)];

elseif opts.method == "corr"
    r(3) = "tmp = rga";
    r(4) = "tmp$p1 = rga$p2";
    r(5) = "tmp$p2 = rga$p1";
    r(6) = "tmp$p1_category = rga$p2_category";
    r(7) = "tmp$p2_category = rga$p1_category";
    r(8) = "rg = rbind(rga, tmp)";
    r(9) = "x2 = data.table::dcast(rg, p1 ~ p2, value.var = 'rg')";
    r(10) = "mat2 = as.matrix(x2[, 2:ncol(x2)])";
    r(11) = "rownames(mat2) = x2$p1";
    r(12) = "mat2[mat2 > 1] = 1";
    r(13) = "mat2[mat2 < -1] = -1";
    r(14) = "mat2[is.na(mat2)] = 0";
    r(15) = "x2 = data.table::dcast(rg, p1 ~ p2, value.var = 'q')";
    r(16) = "qmat2 = as.matrix(x2[, 2:ncol(x2)])";
    r(17) = "rownames(qmat2) = x2$p1";
    r(18) = "qmat2[is.na(qmat2)] = 1";
    r(19) = "if (nrow(mat2) == ncol(mat2)) {";
    r(20) = "  diag(mat2) = 1";
    r(21) = "  diag(qmat2) = -1}";
    r(22) = "png('" + opts.out + "', width = " + opts.width + ...
        ", height = " + opts.height + ", res = " + opts.res + ", units = 'in')";
    r(23) = "corrplot::corrplot(mat2, method = 'psquare', order = '" + opts.order + "'," + ...
        "p.mat = qmat2, sig.level = 0.05, sig = '" + opts.sig + "', pch = '*'," + ...
        "pch.cex = 1.5, full_col=FALSE, na.label = 'square', " + ...
        "na.label.col = 'grey30', diag = F, tl.col = 'black', type='" +...
        opts.type + "', addCoef.col = " + opts.addCoef + ",title='" + opts.title + "')";
    r(24) = "dev.off()";
    
    if opts.pdf
        r(25) = "grDevices::cairo_pdf('" + regexprep(opts.out, ".png$", ".pdf") + "', width = " + opts.width + ...
            ", height = " + opts.height + ")";
        r(26) = "corrplot::corrplot(mat2, method = 'psquare', order = '" + opts.order + "'," + ...
            "p.mat = qmat2, sig.level = 0.05, sig = '" + opts.sig + "', pch = '*'," + ...
            "pch.cex = 1.5, full_col=FALSE, na.label = 'square', " + ...
            "na.label.col = 'grey30', diag = F, tl.col = 'black', type='"...
            + opts.type + "', addCoef.col = " + opts.addCoef + ",title='" + opts.title + "')";
        r(27) = "dev.off()";
    end
end

r = addRfont(r, fontFamily=fontname);
MATLAB2Rconnector("call_gcorrplot.r", code=r);

arrayfun(@delete, ["rg_gcorrplot", "traits_gcorrplot"] + ".txt")

end % END
