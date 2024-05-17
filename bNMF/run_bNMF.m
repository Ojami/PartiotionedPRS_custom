function res = run_bNMF(data, opts)
%@21MARCH2024: a wrapper for bNMF clustering for snp-traits summary
% statistics. The input must be a table with columns as traits and rows as
% independent snps. vmatrix sample size 'N_mat' with the same dimension of
% input table should be provided in tab.Properties.UserData.
% see:
% https://github.com/gwas-partitioning/bnmf-clustering/tree/master
% https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002654#sec008
% https://www.nature.com/articles/s41591-024-02865-3
% 
% Oveis Jamialahmadi, University of Gothenburg, Sweden.
% 
% Note: based on https://arxiv.org/pdf/1111.6085.pdf, K and K0 are the
% same.

arguments
    data {mustBeA(data, "table")}

    opts.run_bNMF {mustBeFile} = fullfile(fileparts(which("run_bNMF.m")), "run_bNMF.R")
    opts.prep_z {mustBeFile} = fullfile(fileparts(which("run_bNMF.m")), "prep_z_matrix.R")

    % prep_z inputs
    % opts.N_mat
    opts.corr_cutoff (1,1) double = 0.85
    opts.pval_cutoff (1,1) double = nan % 0.05 % nan: applies Bonferroni adjustment and set z-score of those below cut-off to zero

    % run_bNMF inputs
    opts.n_reps (1,1) double = 25 % a series of times to generate results and evaluate cluster stability
    opts.random_seed (1,1) double = 123 % random seed
    opts.tol (1,1) double = 1e-7 % Tolerance for convergence of fitting procedure
    opts.K (1,1) double = 5 % Number of clusters to be initialized (algorithm may drive some to zero)
    opts.K0 (1,1) double = 5 % Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
    
    % non-native arguments: output path
    opts.out {mustBeTextScalar} = pwd

    % for formatting results after line 40
    opts.width (1,1) double = 7
    opts.height (1,1) double = 8
end

if isempty(opts.out) || opts.out == ""
	opts.out = pwd;
end
if ~isfolder(opts.out), mkdir(opts.out); end
opts.out = regexprep(opts.out, filesep + "$", "");
opts.rout = replace(opts.out, filesep, "/"); % R path

file = fullfile(opts.out, getRandomName("bNMF_R", 5)); % random tmp file name

% get row and col names: missing values will be set to 0 by prep_z
rownames = data.Row;
if isempty(rownames)
    rownames = "y" + (1:height(data));
end
rnames = cellstr(rownames);
cnames = cellstr(colnames(data));
n_mat = data.Properties.UserData;
data = data{:, :};
assert(all(size(data) == size(n_mat)))
save(file + ".mat", "data", "rnames", "cnames", "n_mat")

% R script is defined in string 'r'; modify it as appropriate.
r = string;
r(1) = "setwd('" + opts.rout + "')";
r(2) = "source('" + opts.run_bNMF.replace("\", "/") + "')";
r(3) = "source('" + opts.prep_z.replace("\", "/") + "')";
r(4) = "tmp = R.matlab::readMat('" + replace(file + ".mat", filesep, "/") + "')";
r(5) = "data = tmp$data";
r(6) = "rownames(data) = unlist(tmp$rnames)";
r(7) = "colnames(data) = unlist(tmp$cnames)";
r(8) = "n_mat = tmp$n.mat";
r(9) = "rownames(n_mat) = unlist(tmp$rnames)";
r(10) = "colnames(n_mat) = unlist(tmp$cnames)";


% prep Z
if isnan(opts.pval_cutoff), opts.pval_cutoff = "NULL"; end
r(12) = "prep_z_output = prep_z_matrix(z_mat = data, N_mat = n_mat,corr_cutoff = " +...
    opts.corr_cutoff + ", pval_cutoff=" + opts.pval_cutoff + ")";
r(13) = "final_zscore_matrix <- prep_z_output$final_z_mat"; % The scaled, non-negative z-score matrix
r(14) = "df_traits_filtered <- prep_z_output$df_traits"; % Results from the trait filtering
r(15) = "readr::write_csv(x = df_traits_filtered,file ='df_traits.csv')";

% run bNMF algorithm
r(17) = "bnmf_reps = run_bNMF(final_zscore_matrix, n_reps=" + opts.n_reps + ...
    ",tolerance = " + opts.tol + ", random_seed = " + opts.random_seed + ...
    ",K=" + opts.K + ", K0=" + opts.K0 + ")";
r(18) = "summarize_bNMF(bnmf_reps, dir_save='" + opts.rout + "')";

% format results: most probable K and optimal cut-off point
r(25) = "pckg = c(""dplyr"", ""data.table"", ""ggplot2"", ""ggrepel"", ""magrittr"", ""strex"")";
r(26) = "lapply(pckg, function(x) suppressMessages(require(x, character.only = TRUE)))";
r(29) = "df_run_summary <- fread('run_summary.txt', stringsAsFactors = F, data.table = F)";
r(30) = "k_counts <- df_run_summary %>% dplyr::count(K) %>% mutate(perc = n/nrow(df_run_summary))";
r(31) = "k <- k_counts$K[which.max(k_counts$n)]";
r(32) = "w <- fread(sprintf(""L2EU.W.mat.%i.txt"", k), stringsAsFactors = FALSE, data.table = F)";
r(33) = "n_variants <- nrow(w)";
r(34) = 'h <- fread(sprintf("L2EU.H.mat.%i.txt",k), stringsAsFactors = FALSE, data.table = F) %>% rename_at(vars(starts_with("X2hr")), function(x){gsub("X","",x) })';
r(35) = "orig_traits <- colnames(h)";
r(36) = "traits_short <- str_before_last(names(h), ""[_]"")";
r(37) = "df_traits <- data.frame(Trait=unique(traits_short))";
r(38) = 'trait_results <- fread("df_traits.csv")';
r(39) = "write.table(k_counts, file='" + opts.rout + "/k_counts.txt', sep = '||', row.names = F, col.names = T, quote = F)";

% Determining the optimal K from the bNMF run summary
r(40) = "run_sum_max <- df_run_summary %>% group_by(K) %>% summarise(convergence = min(evid, na.rm=TRUE)) %>% mutate(K=as.factor(K))";
r(42) = "pdf('" + opts.rout + "/conv_score_by_k.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")";   
r(43) = 'ggplot(data=run_sum_max, aes(x=K, y=convergence, fill=K)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Spectral") + labs(y="Convergence") + theme_minimal() + theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), plot.title = element_text(size = 24, face = "bold"))+ ggtitle("Convergence Scores by K")';
r(44) = "dev.off()";
r(46) = 'df <- k_counts %>% mutate(value=n,group=paste("K",K,sep="=")) %>% arrange(desc(value))';
r(47) = "df2 <- df %>% mutate(csum = rev(cumsum(rev(value))), pos = value/2 + lead(csum, 1), pos = if_else(is.na(pos), value/2, pos))";

r(49) = "pdf('" + opts.rout + "/dist_of_k_converge.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")"; 
r(50) = "ggplot(df, aes(x = '', y = value, fill = forcats::fct_inorder(group))) " + ...
    "+ geom_col(width = 1, color = 1) + coord_polar(theta = 'y') + " + ...
    "scale_fill_brewer(palette = 'Spectral')+ " + ...
    "geom_label_repel(data = df2, aes(y = pos, label = sprintf('K=%i\n(%i%%)',K,round(perc*100))), "+ ...
    " size = 8, nudge_x = 1, show.legend = FALSE) + theme_void()+ " + ...
    "scale_y_continuous(breaks = df2$pos, labels = df$group) + " + ...
    "theme(plot.title = element_text(size = 20, face = 'bold'), " + ...
    "legend.position = 'none') + ggtitle(sprintf('Distribution of K convergences\n%i iterations', sum(df2$n)))";
r(51) = "dev.off()";

% Determining the optimal cutoff for clustering weights (parameters: N, M)
% 1) Fit a 1st line to the top N% of weights
% 2) Fit a 2nd line to the M% of tail weight
% 3) Using the remaining weights from top N% to last M%, check if they have shorter distance to 1st or 2nd line
% 4) The first weight that has a shorter distance to 2nd line (defined by long tail) is selected as the cutoff.
r(54) = 'w_mat <- data.frame(t(w)) %>% set_colnames(.["variant", ]) %>%subset(!rownames(.) %in% "variant")';
r(56) = "dist_point_line <- function(a, slope, intercept) {";
r(57) = "  b = c(1, intercept+slope)";
r(58) = "  c = c(0, intercept)";
r(59) = "  v1 <- b - c";
r(60) = "  v2 <- a - b";
r(61) = "  m <- cbind(v1,v2)";
r(62) = "  return(abs(det(m))/sqrt(sum(v1*v1)))}";

r(64) = "w_list <- as.list(w_mat)";
r(65) = "weight = unlist(w_list)";

r(66) = "weights = as.data.frame(weight)";
r(67) = "weights$weight = as.numeric(as.character(weights$weight))";
r(68) = "weights = as.data.frame(weights[order(weights$weight,decreasing = T),])";
r(69) = 'names(weights)[1] = "weight"';

r(70) = "numWeights <- dim(weights)[1]";
r(71) = "weights$x = c(1:numWeights)";
r(72) = "weights$x_scale = weights$x/numWeights";
r(73) = "weights$w_scale = weights$weight/max(weights$weight)";

r(74) = "cutoff_test <- quantile(weights$weight, prob=1-1/100)";
r(75) = "n = ifelse(sum(weights$weight > cutoff_test)>3, 1, 5)";
r(76) = "m = 80";
r(77) = "top1 = as.data.frame(weights[weights$weight > quantile(weights$weight, prob=1-n/100),])";
r(78) = "low80 = as.data.frame(weights[weights$weight <= quantile(weights$weight, prob=m/100),])";
r(79) = 'weights$group = ifelse(weights$x %in% top1$x,paste0("top",n), ifelse(weights$x %in% low80$x,paste0("bottom",m),"other"))';

r(80) = "top1line <- lm(w_scale ~ x_scale, data=top1)";
r(81) = "bottom80line <- lm(w_scale ~ x_scale, data=low80)";

r(82) = "weights$dist_top1 = rep(0,numWeights)";
r(83) = "for(i in 1:numWeights){";
r(84) = "  x = weights$x_scale[i]";
r(85) = "  y = weights$w_scale[i]";
r(86) = "  weights$dist_top1[i] = dist_point_line(c(x,y),top1line$coefficients[2],top1line$coefficients[1])";
r(87) = "}";

r(88) = "weights$dist_low80 = rep(0,numWeights)";
r(89) = "for(i in 1:numWeights){";
r(90) = "  x = weights$x_scale[i]";
r(91) = "  y = weights$w_scale[i]";
r(92) = "  weights$dist_low80[i] = dist_point_line(c(x,y),bottom80line$coefficients[2],bottom80line$coefficients[1])";
r(93) = "}";

r(94) = "weights$diff = weights$dist_top1 - weights$dist_low80";
r(95) = "cut = weights[weights$diff > 0,]";
r(96) = "cutoff = cut$weight[1]";
r(97) = "cutoff_num = cut$x[1]";

r(98) = "highlight <- weights %>% filter(x == cutoff_num)";

% plot optimal cut-off plot
r(100) = "pdf('" + opts.rout + "/optimal_weight_cutoff.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")"; 
r(101) = "ggplot(weights, aes(x=x_scale, y=w_scale, color=group)) + " + ...
    "geom_point() + geom_abline(intercept = top1line$coefficients[1], " + ...
    "slope = top1line$coefficients[2], color='blue') + " + ...
    "geom_abline(intercept = bottom80line$coefficients[1], " + ...
    "slope = bottom80line$coefficients[2], color='red') + " + ...
    "geom_point(data=highlight, aes(x=x_scale,y=w_scale), " + ...
    "color='red',size=5) + xlab('Scaled count') + ylab('Scaled cluster weight')";
r(102) = "dev.off()";

r(104) = "cutoff_f = format(cutoff, digits = 5)";
r(105) = "cutoff_num_f = format(100*cutoff_num/numWeights, digits=2)";
r(106) = "write(cutoff, '" + opts.rout + "/optimal_cutoff.txt')";

r(107) = "w <- w %>% mutate(gene=variant)";
r(108) = 'cluster_names <- names(w)[names(w) %like% "X"]';

 
% variant * clusters
r(109) = "suppressPackageStartupMessages(library(ComplexHeatmap))";
r(110) = "library(magrittr)";
r(111) = "suppressPackageStartupMessages(library(circlize))";
% assign column label colors based on which cluster they have their maximum weight
r(112) = "palette1_named = setNames(object = scales::hue_pal()(length(cluster_names)), nm = cluster_names)";
r(113) = "ix=sapply(w_mat, which.max)";
r(114) = "label_colors=palette1_named[ix]";
r(115) = 'hm_colors = colorRamp2(c(0, max(weights$weight)), c("white", "black"))';

r(119) = "pdf('" + opts.rout + "/lociXclusters.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")"; 
r(120) = "w %>% set_rownames(make.unique(.[,'gene'])) %>%   dplyr::select(starts_with('X')) %>% t() %>% as.matrix() %>% Heatmap(name = 'Weight', column_title='Loci x Clusters', border_gp = gpar(col = 'black'), show_column_dend = F,column_names_gp = grid::gpar(fontsize = 6),col=hm_colors)";
r(121) = "dev.off()";

% repeat for phenotypes
r(123) = "ix=sapply(h, which.max)";
r(124) = "label_colors=palette1_named[ix]";
r(126) = "pdf('" + opts.rout + "/traitXclusters.pdf', width = " + opts.width + ...
        ", height = " + opts.height + ")"; 
r(127) = "h %>% set_rownames(cluster_names) %>% as.matrix() %>% Heatmap(name='Weight', column_title='Phenotypes x Clusters', border_gp = gpar(col = 'black'), show_column_dend = F, column_names_gp = grid::gpar(fontsize = 8), col=hm_colors)";
r(128) = "dev.off()";

% Cluster Circle Plots
% Only includes variants and phenotypes with weights above cutoff
r(130) = "library(readr)";
r(131) = "loci = w[cluster_names]";
% # colnames(loci) <- c(cluster_names, "name")
r(133) = "loci$name <- make.unique(w$gene)";
r(134) = 'loci$group <- as.factor("Variant")';

r(136) = "trait = data.frame(t(h))";
% # colnames(trait) <- (c("name",cluster_names))
% # trait[cluster_names] <- sapply(trait[cluster_names],as.numeric)
r(139) = 'trait$name = gsub("_neg|_pos","",rownames(trait))';
r(140) = 'trait$group <- as.factor(ifelse(grepl("_pos",rownames(trait)), "Pos_Trait", ifelse(grepl("_neg",rownames(trait)), "Neg_Trait","Misc. Trait")))';

r(142) = "for (i in 1:length(cluster_names)) { #1:length(cluster_names)";
r(143) = "  my_clust=cluster_names[i]";
r(144) = '  df_loci <- loci[loci[my_clust]>=cutoff, c(my_clust, "name","group")]';
r(145) = '  df_trait <- trait[trait[my_clust]>=cutoff, c(my_clust, "name","group")]';
r(146) = '  data <- rbind(df_loci, df_trait) %>% set_colnames(c("value","name","group"))';
r(147) = "  if (nrow(data)>40){";
    % # data = data %>% group_by(group) %>% top_n(10,wt=value) %>%
    % #   arrange(group, desc(value))
r(150) = "    data = data %>% top_n(40,wt=value) %>% arrange(group, desc(value))";
r(151) = "  } else{";
r(152) = "      data = data %>% group_by(group) %>% arrange(group, desc(value))";
r(153) = "      }";
r(154) = '  data$value = ifelse(data$group=="Neg_Trait", data$value*-1.5, data$value*1.5)';%

  % # set number of empty bars
r(156) = "  empty_bar <- 1";
  % # add lines to the initial dataset
r(157) = "  to_add <- data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )";
r(158) = "  colnames(to_add) <- colnames(data)";
r(159) = "  to_add$group <- rep(levels(data$group), each=empty_bar)";
r(160) = "  data <- rbind(data, to_add)";
r(161) = "  data <- data %>% arrange(group)";
r(162) = "  data$id <- seq(1, nrow(data))";
r(163) = "  ymax=(max(data$value,na.rm=T))";
r(164) = "  ymin=(min(data$value,na.rm=T))";

  % # get the name and the y position of each label
r(165) = "  label_data <- data";
r(166) = "  number_of_bar <- nrow(label_data)";
r(167) = "  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)";
r(168) = "  label_data$hjust <- ifelse( angle < -90, 1, 0)";
r(169) = "  label_data$angle <- ifelse(angle < -90, angle+180, angle)";
  
  % # prepare a data frame for base lines
r(171) = "  base_data <- data %>% group_by(group) %>% dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% rowwise() %>% mutate(title=mean(c(start, end)))";
  
  % # make the plot
r(172) = "  p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) " + ...
    " + geom_bar(stat='identity', alpha=0.7,width=0.8) + " + ...
    "ylim(-20, ymax+75) + theme_minimal() + theme(legend.position = 'none', " + ...
    "plot.title = element_text(hjust = 0.5,vjust = -10), " + ...
    "axis.text = element_blank(), axis.title = element_blank(), " + ...
    "panel.grid = element_blank(), plot.margin = unit(rep(-1,4), 'cm')) " + ...
    "+ scale_fill_manual('Legend', values =c('Variant' = 'green3', " + ...
    "'Neg_Trait' = 'blue3', 'Pos_Trait' = 'red3'))+ coord_polar(clip = 'off') " + ...
    "+ geom_text(data=label_data, aes(x=id, y=ifelse(value>0, value+5, 7)," + ...
    "label=name, hjust=hjust), color='black', fontface='bold', alpha=0.7, size=4, angle= label_data$angle, inherit.aes = FALSE )";
  
    % #y=ifelse(value>=0,value+10,10)
    % # Add base line information
    % # geom_segment(data=base_data, aes(x = start-0.5, y = 0, xend = end+0.5, yend = 0), colour = "darkslategrey", alpha=0.8, size=0.2 , inherit.aes = FALSE )  +
    % # geom_text(data=base_data, aes(x = title, y = -22, label=group), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE) 
  
r(173) = "  ggsave(sprintf('/circlePlot_k%i_cluster%i.png', k, i), plot=p, path = getwd(), width = 7, height = 7, dpi = 300, units = 'in', device='png')";
r(174) = "  print(p+ggtitle(sprintf('Cluster %s',i)))";
r(175) = "}";


MATLAB2Rconnector(file + ".r", code=r, delr=true, log=true);
delete(file + ".mat")

% read optimal cut-off
if ~isfile(fullfile(opts.out, "optimal_cutoff.txt"))
    res = struct;
    return % check log file: ran into error probably because W and H are non-conformable so cannot perform W %*% H in summarize_bNMF func
end
optimal_cutoff = readmatrix(fullfile(opts.out, "optimal_cutoff.txt"));

% find most probable K
optimal_k = readtable(fullfile(opts.out, "k_counts.txt"), "TextType", "string",...
    "Delimiter", "||", "ReadVariableNames", true, "NumHeaderLines", 0);
optimal_k = sortrows(optimal_k, "perc", "descend");

% read variant weights for most probable K
v_tab = readtable(fullfile(opts.out, "L2EU.W.mat." + optimal_k.K(1) + ".txt"), ...
    "TextType", "string", "ReadVariableNames", true, "NumHeaderLines", 0);
v_tab = movevars(v_tab, "variant", Before=1);
idx = v_tab{:, 2:end} < optimal_cutoff;
v_tab_clust = v_tab;
for k = 2:width(v_tab)
    v_tab_clust{idx(:, k-1), k} = 0;
end

% read trait weights for most probable K
t_tab = readtable(fullfile(opts.out, "L2EU.H.mat." + optimal_k.K(1) + ".txt"), ...
    "TextType", "string", "ReadVariableNames", true, "NumHeaderLines", 0);
t_tab = rows2vars(t_tab);
t_tab.Properties.VariableNames = ["Pheno", "X" + (1:(width(t_tab)-1))];
t_tab.Pheno = string(t_tab.Pheno);

res = struct;
res.optimal_cutoff = optimal_cutoff;
res.optimal_k = optimal_k;
res.v_tab = v_tab;
res.v_clust = v_tab_clust;
res.t_tab = t_tab;



end % END