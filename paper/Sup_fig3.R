source("figure_code/estimate_plotter.R")
source("../simulations/sim_distribution.R")
source("../packages.R")
source("../mrema.R")
j <- read.csv("../simulations/BRCA_sim_param.csv")
library(ssizeRNA)
library(cowplot)

axis_title_size <- 22
axis_text <- 18
legend_text <- 16

set.seed(6)
### Estimated plots
sd4_estimate_N20 <- estimate_plots(0.6, 0.4, 20, thresh1 = 1.5, y_top = 1.5, limit = 2) 
sd4_estimate_N40 <- estimate_plots(0.6, 0.4, 40, thresh1 = 1.5, y_top = 1.5, limit = 2) 

underlying_lfc_plot <- function(data, y_top, thresh){
  ggplot(data = data, aes(LFC, Density, group = Set, color = Set)) +
    geom_line() + xlab("Log Fold Change") + ylab("Density") + 
    geom_area(aes( group = Set, fill = Set), alpha = 0.4, position = "identity") +
    scale_color_manual(values = c("#DC3220", "#005AB5")) +
    scale_fill_manual(values = c("#DC3220", "#005AB5")) +
    geom_segment(aes(x = log2(thresh), y = y_top, xend = log2(thresh), yend = 0), linetype=2, color =  "black", size = 1, alpha = 0.5) + 
    geom_segment(aes(x = -log2(thresh), y = y_top, xend = -log2(thresh), yend = 0), linetype=2, color =  "black", size = 1, alpha = 0.5) +
    theme_minimal() +
    
    ggtitle("A) Simulated LFC") +
    theme(plot.title = element_text(hjust = 0, size = 22),
          axis.title=element_text(size=axis_title_size),
          axis.text = element_text(size = axis_text),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = legend_text),
          legend.title = element_text(size = legend_text),
          legend.position = "bottom")}


SD_in <- 0.6
SD_out <- 0.4
x <- seq(-2,2, by = 0.01)
y_in <- dnorm(x, 0, SD_in)
y_out <- dnorm(x, 0, SD_out)
group <- c(rep("Inset", length(x)), rep("Outset", length(x)))
y <- c(y_in, y_out)
LFC <- c(x,x)
data <- data.frame("Set" = group, "LFC" = LFC, "Density" = y)
SDo_0.4_SDi_0.6_plot <- underlying_lfc_plot(data,  1, 1.5)

sup_fig3 <- plot_grid(SDo_0.4_SDi_0.6_plot, 
                      sd4_estimate_N20 + ggtitle("B) Estimated underlying LFC distribution - 20 samples") + theme(plot.title = element_text(hjust = 0, size = 22)), 
                      sd4_estimate_N40 + ggtitle("C) Estimated underlying LFC distribution - 40 samples") + theme(plot.title = element_text(hjust = 0, size = 22)), 
                      nrow = 3, align = "v")


png(filename = "Sup_fig3.png", width = 1500, height = 2000)
sup_fig3
dev.off()
