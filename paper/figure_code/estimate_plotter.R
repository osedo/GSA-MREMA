
estimate_plots <- function(sd_in, sd_out, N, limit, thresh1, y_top){

fc <- 2^rhalfnorm(100, sd2theta(sd_in))
simulation <- simulation_dist(s = 100, p = 100, N = N, beta = 0.01, gamma = 1, theta = 1, foldchange = fc, inset_fc = sd_out, upreg = 0.5)
# filter out genes not in list of gene sets and lowly expressed genes
igenes <- intersect(rownames(assay(simulation$sum_exp_raw_count)), unique(unlist(simulation$raw.gs)))
if(!length(igenes)) stop("Expression dataset (se)", " and ", "gene sets (gs) have no gene IDs in common")
se <- assay(simulation$sum_exp_raw_count)[igenes,]
group <- colData(simulation$sum_exp_raw_count)$GROUP
keep <- filterByExpr(se, group = group)
se <- se[keep,]
gene_sets <- lapply(simulation$raw.gs, function(s) s[s %in% rownames(se)]) 

### use Deseq2 to get group effect estimate and it's standard error
dea_raw_count<- DESeqDataSetFromMatrix(countData = se, colData = colData(simulation$sum_exp_raw_count), design = ~ GROUP)
dds <- DESeq(dea_raw_count)
res <- results(dds)
postdata <- tibble("Ensembl" = rownames(res), "effect" = res[,2], "variance" = res[,3]^2)
postdata_sim <- postdata[complete.cases(postdata),]

lower_thresh <- mrema(postdata_sim, gene_sets, set_number = 100, DF =6, threshold = thresh1, params = TRUE)
inset <- lower_thresh[c(1,2,3)] 
g <- seq(-limit,limit, by=.01)
lower_value <- inset$alpha[1]*dnorm(g, inset$mu[1], sqrt(inset$var[1])) + inset$alpha[2]*dnorm(g, inset$mu[2], sqrt(inset$var[2])) + inset$alpha[3]*dnorm(g, inset$mu[3], sqrt(inset$var[3]))

value <- c(lower_value)
method <- c(rep("T1.5_DF7", length(g)))
LFC <- rep(g, 1)
data_estimate <- data.frame("Method" = method, "LFC" = LFC, "Density" = value)

plot = .estimate_lfc_plot(data_estimate, thresh1,  y_top)

return(plot)
}


cbbPalette <-c("#00BA38", "#619CFF")

underlying_lfc_plot <- function(data, y_top){
  ggplot(data = data, aes(LFC, Density, group = Set, color = Set)) +
    geom_line() + xlab("Log Fold Change") + ylab("Density") + 
    geom_area(aes( group = Set, fill = Set), alpha = 0.4, position = "identity") +
    scale_color_manual(values = c("#DC3220", "#005AB5")) +
    scale_fill_manual(values = c("#DC3220", "#005AB5")) +
    geom_segment(aes(x = log2(1.1), y = y_top, xend = log2(1.1), yend = 0), linetype=2, color =  "#1E88E5", size = 1, alpha = 0.5) + 
    geom_segment(aes(x = -log2(1.1), y = y_top, xend = -log2(1.1), yend = 0), linetype=2, color =  "#1E88E5", size = 1, alpha = 0.5) +
    geom_segment(aes(x = log2(1.3), y = y_top, xend = log2(1.3), yend = 0), linetype=2, color = "#FFC107", size = 1, alpha = 0.5) + 
    geom_segment(aes(x = -log2(1.3), y = y_top, xend = -log2(1.3), yend = 0), linetype=2, color = "#FFC107", size = 1, alpha = 0.5) +
    geom_segment(aes(x = log2(1.5), y = y_top, xend = log2(1.5), yend = 0), linetype=2, color = "#004D40", size = 1, alpha = 0.5) + 
    geom_segment(aes(x = -log2(1.5), y = y_top, xend = -log2(1.5), yend = 0), linetype=2, color = "#004D40", size = 1, alpha = 0.5) +
    theme_minimal() +
    
    ggtitle("A) Simulated LFC") +
    theme(plot.title = element_text(hjust = 0, size = 22),
          axis.title=element_text(size=axis_title_size),
          axis.text = element_text(size = axis_text),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = legend_text),
          legend.title = element_text(size = legend_text),
          legend.position = "bottom")}

.estimate_lfc_plot <- function(data, thresh1, y_top){
  ggplot(data = data, aes(LFC, Density, group = Method, color = Method)) +
    geom_line() + xlab("Log Fold Change") + ylab(NULL) + 
    geom_area(aes( group = Method, fill = Method), alpha = 0.4, position = "identity") +
    geom_segment(aes(x = log2(thresh1), y = y_top, xend = log2(thresh1), yend = 0), linetype=2, color = "black", size = .35, alpha = 0.1) + 
    geom_segment(aes(x = -log2(thresh1), y = y_top, xend = -log2(thresh1), yend = 0), linetype=2, color = "black", size = .35, alpha = 0.1) +
    scale_colour_manual(values = "#DC3220") + 
    scale_fill_manual(values = "#DC3220") +
    ylim(c(0,y_top)) +
    theme_minimal() +
    
    theme(plot.title = element_text(hjust = 0.5, size = 22),
          axis.title=element_text(size=axis_title_size),
          axis.text = element_text(size = axis_text),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = legend_text),
          legend.title = element_text(size = legend_text),
          legend.position = "bottom")}