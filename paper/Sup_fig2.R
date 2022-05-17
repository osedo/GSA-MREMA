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


sd1_estimate_N20 <- estimate_plots(0.25, 0.1, 20, thresh1 = 1.1, y_top = 3.5, limit = 1) + scale_fill_manual(name = "Method", labels = c("T1.1_DF7"), values = c("#DC3220")) +  scale_colour_manual(name = "Method", labels = c("T1.1_DF7"), values = c("#DC3220"))
sd1_estimate_N40 <- estimate_plots(0.25, 0.1, 40, thresh1 = 1.1, y_top = 3.5, limit = 1)  + scale_fill_manual(name = "Method", labels = c("T1.1_DF7"), values = c("#DC3220")) +  scale_colour_manual(name = "Method", labels = c("T1.1_DF7"), values = c("#DC3220"))


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

SD_in <- 0.25
SD_out <- 0.1
x <- seq(-1,1, by = 0.01)
y_in <- dnorm(x, 0, SD_in)
y_out <- dnorm(x, 0, SD_out)
group <- c(rep("Inset", length(x)), rep("Outset", length(x)))
y <- c(y_in, y_out)
LFC <- c(x,x)
data <- data.frame("Set" = group, "LFC" = LFC, "Density" = y)
SDo_0.1_SDi_0.25_plot <- underlying_lfc_plot(data,  4 , 1.1)





sup_fig2 <- plot_grid(SDo_0.1_SDi_0.25_plot, 
                      sd1_estimate_N20 + ggtitle("B) Estimated underlying LFC distribution - 20 samples") + theme(plot.title = element_text(hjust = 0, size = 22)), 
                      sd1_estimate_N40 + ggtitle("C) Estimated underlying LFC distribution - 40 samples") + theme(plot.title = element_text(hjust = 0, size = 22)), 
                      nrow = 3, align = "v")

png(filename = "Sup_fig2.png", width = 1500, height = 2000)
sup_fig2
dev.off()