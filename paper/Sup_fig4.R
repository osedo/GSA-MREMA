load("Sup_fig4.RData")
library(cowplot)
library(ggplot2)

axis_title_size <- 22
axis_text <- 18
legend_text <- 16
ora_colors <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40")
results_plot <- function(data){
  plot <- ggplot(data = data, aes(x = Inset, y = Power , group = Method, color = Method)) + 
    geom_line(size = 1.5) +
    
    geom_point(size = 1.5) + xlab("Inset standard deviation") +
    ylim(c(0,1)) +
    #scale_colour_manual(values = cbbPalette) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0, size = 22),
          axis.title=element_text(size=axis_title_size),
          legend.key.size = unit(1, 'cm'),
          axis.text = element_text(size = axis_text),
          legend.text = element_text(size = legend_text),
          legend.title = element_text(size = legend_text )) 
  
  plot = plot +  geom_errorbar(aes(ymin = Power-ME, ymax = Power+ME), width=.02,size = 1.5) 
  return(plot)
}

plot_ora_sd1_20 <- results_plot(SD_01_data_N20_ora) + scale_color_manual(values = ora_colors, labels = c("None", "FC 1.1", "FC 1.3", "FC 1.5")) + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank())
plot_ora_sd1_40 <- results_plot(SD_01_data_N40_ora) + scale_color_manual(values = ora_colors, labels = c("None", "FC 1.1", "FC 1.3", "FC 1.5")) + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank())
plot_ora_sd4_20 <- results_plot(SD_04_data_N20_ora) + scale_color_manual(values = ora_colors, labels = c("None", "FC 1.1", "FC 1.3", "FC 1.5"))
plot_ora_sd4_40 <- results_plot(SD_04_data_N40_ora) + scale_color_manual(values = ora_colors, labels = c("None", "FC 1.1", "FC 1.3", "FC 1.5")) + theme(legend.position = "none",axis.title.y = element_blank(), axis.text.y = element_blank())

l2 <- get_legend(plot_ora_sd4_20 + guides(color = guide_legend(title.position = "left", title = "Threshold", ncol = 4, override.aes = list(shape = c(15), linetype = NULL, size = 10)))) 
power_ora <- plot_grid(plot_ora_sd1_20 + ggtitle("A) Outset standard deviation 0.1 - 20 samples"),
                       plot_ora_sd1_40 + ggtitle("B) Outset standard deviation 0.1 - 40 samples"),
                       plot_ora_sd4_20 + ggtitle("C) Outset standard deviation 0.4 - 20 samples") + theme(legend.position = "none"),
                       plot_ora_sd4_40 + ggtitle("D) Outset standard deviation 0.4 - 40 samples"),
                       nrow = 2)
power_ora <- plot_grid(power_ora, l2, nrow = 2, rel_heights = c(0.9, 0.1))


png(filename = "Sup_fig4.png", width = 1500, height = 1500)
power_ora
dev.off()