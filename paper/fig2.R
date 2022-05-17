library(ggplot2)
library(cowplot)
library(gt)
library(png)
library(grid)
library(gridExtra)
library(dplyr)
load("fig2.RData")

SD_01_data_N40 <- SD_01_data_N40[-c(16,17,18,19,20),]
SD_04_data_N40 <- SD_04_data_N40[-c(16,17,18,19,20),]

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


ora_colors <- c( "#1E88E5", "#FFC107", "#004D40")


plot_sd4_40 <- results_plot(SD_04_data_N40)  + scale_color_manual(values = ora_colors, labels = c( "FC-1.1", "FC-1.3", "FC-1.5")) + theme(  axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

plot_sd1_40 <- results_plot(SD_01_data_N40)  + scale_color_manual(values = ora_colors, labels = c( "FC-1.1", "FC-1.3", "FC-1.5")) 

l2 <- get_legend(plot_sd4_40 + guides(color = guide_legend(title.position = "left", title = "Threshold", ncol = 4, override.aes = list(shape = c(15), linetype = NULL, size = 10))))

power <- plot_grid(plot_sd1_40  + ggtitle("C) Power - Low perturbation") + theme(legend.position = "none"), 
                   plot_sd4_40  + ggtitle("D) Power - High perturbation") + theme(legend.position = "none"), 
                   nrow = 1, align = c("h", "v"), rel_widths = c(1.12, 1))
power <- plot_grid(power, l2, nrow = 2, rel_heights = c(0.9, 0.1))


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
    
    ggtitle("B) Simulated LFC - High perturbation") +
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
SDo_0.4_SDi_0.6_plot <- underlying_lfc_plot(data,  4)


SD_in <- 0.25
SD_out <- 0.1
x <- seq(-2,2, by = 0.01)
y_in <- dnorm(x, 0, SD_in)
y_out <- dnorm(x, 0, SD_out)
group <- c(rep("Inset", length(x)), rep("Outset", length(x)))
y <- c(y_in, y_out)
LFC <- c(x,x)
data <- data.frame("Set" = group, "LFC" = LFC, "Density" = y)
SDo_0.1_SDi_0.25_plot <- underlying_lfc_plot(data,  4) + ggtitle("A) Simulated LFC - Low perturbation")

l1 <- get_legend(SDo_0.1_SDi_0.25_plot)

simulation <- plot_grid(SDo_0.1_SDi_0.25_plot + theme(legend.position = "none"), 
                        SDo_0.4_SDi_0.6_plot + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()), 
                        nrow = 1, align = c("h", "v") , rel_widths = c(1.12, 1))
simulation <- plot_grid(simulation, l1, nrow = 2, rel_heights = c(0.9, 0.1))

fig1 <- plot_grid(simulation, power, nrow = 2)
png(filename = "fig2.png", width = 1250, height = 1250)
fig1
dev.off()
