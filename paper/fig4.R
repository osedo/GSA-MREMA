source("../mrema.R")
source("../cancer_analysis/GMM_samples.R")
source("../packages.R")
load("fig4.RData")
library(cowplot)

# mTOR
mtor <- plot_distributions(GMM_side_by_side(combined_postdata, gene_set_combined, set_number = 149, threshold = 1.5, DF = 6, params = TRUE), x_lim = c(-5,5), y_lim = c(0,1.25)) + scale_color_manual(labels = c("Young-female", "All other"), values = c('#DC3220', '#005AB5')) +  scale_fill_manual(labels = c("Young-female", "All other"), values = c('#DC3220', '#005AB5'))
# ras
ras <- plot_distributions(GMM_side_by_side(combined_postdata, gene_set_combined, set_number = 120, threshold = 1.5, DF = 6, params = TRUE), x_lim = c(-5,5), y_lim = c(0,1.25))  + scale_color_manual(labels = c("Mentstruating", "Non-menstruating"), values = c('#DC3220', '#005AB5')) +  scale_fill_manual(labels = c("Mentstruating", "Non-menstruating"), values = c('#DC3220', '#005AB5'))
## neuro 
neuro <- plot_distributions(GMM_side_by_side(combined_postdata, gene_set_combined, set_number = 270, threshold = 1.5, DF = 6, params = TRUE), x_lim = c(-5,5), y_lim = c(0,1.25))  + scale_color_manual(labels = c("Mentstruating", "Non-menstruating"), values = c('#DC3220', '#005AB5')) +  scale_fill_manual(labels = c("Mentstruating", "Non-menstruating"), values = c('#DC3220', '#005AB5'))

#all 
all_genes <- list()
all_genes[["all"]] <- unlist(combined_postdata$genesEnsembl)
all_genes <-  plot_distributions(GMM_side_by_side(combined_postdata, all_genes, set_number = 1, threshold = 1.5, DF = 6,  params = TRUE), x_lim = c(-5,5), y_lim = c(0,1))  + scale_color_manual(labels = c("Young-female", "All other"), values = c('#DC3220', '#005AB5')) +  scale_fill_manual(labels = c("Young-female", "All other"), values = c('#DC3220', '#005AB5'))

upper <- plot_grid( mtor + ggtitle("A) All genes") + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()),
                    mtor + ggtitle("B) mTOR signaling pathway") + theme(legend.position = "none" , axis.title = element_blank(), axis.text = element_blank()), 
                    ras  + ggtitle("C) Ras signaling pathway")  + theme(legend.position = "none"), 
                    neuro + ggtitle("D) Pathways of neurodegeneration")  + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()), 
                    nrow =2 , rel_widths = c(1.11, 1))
l2 <- get_legend(mtor  + theme(legend.position = "bottom", legend.title = element_blank()))
figure3 <- plot_grid(upper, l2, nrow = 2, rel_heights = c(1,0.1))

png(filename = "fig4.png", width = 1000, height = 500)
figure3
dev.off()