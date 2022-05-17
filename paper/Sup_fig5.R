
load("Sup_fig5.RData")
library(ggplot2)

axis_title_size <- 22
axis_text <- 18
legend_text <- 16

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

N20 <- results_plot(SD_04_data_N20)
N40 <- results_plot(SD_04_data_N40)

l2 <- get_legend(N20 + guides(color = guide_legend(title.position = "left", title = "Method", ncol = 2, override.aes = list(shape = c(15), linetype = NULL, size = 10)))) 
library(cowplot)

plot <- plot_grid(N20 + ggtitle("20 Samples") + theme(legend.position = "none"),
          N40  + ggtitle("40 Samples") + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()),
          nrow =1, rel_widths = c(1.1, 1))

plot <- plot_grid(plot,
          l2,
          nrow = 2, rel_heights = c(0.95, 0.05))

png(filename = "Sup_fig5.png", width = 1000, height = 500)
plot
dev.off()