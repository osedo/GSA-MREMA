library(ggplot2)
library(cowplot)
library(gt)
library(png)
library(grid)
library(gridExtra)
library(dplyr)
load("fig1.RData")


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


ora_colors <- c(  "#FFC107","#D81B60", "#1E88E5", "#004D40")


plot_sd4_20 <- results_plot(SD_04_data_N20) + scale_color_manual(values = ora_colors, labels = c( "GSEA", "ORA", "SAFE", "1DF-T1.5")) + theme(legend.position = "none")
plot_sd4_40 <- results_plot(SD_04_data_N40)  + scale_color_manual(values = ora_colors, labels = c( "GSEA", "ORA", "SAFE", "1DF-T1.5")) + theme(  axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 


l2 <- get_legend(plot_sd4_40 + guides(color = guide_legend(title.position = "left", title = "Method", ncol = 5, override.aes = list(shape = c(15), linetype = NULL, size = 8))))
# + scale_y_continuous(position = "right")

power <- plot_grid(plot_sd4_20  + ggtitle("C) Power - 20 samples "), 
                   plot_sd4_40  + ggtitle("D) Power - 40 samples") + theme(legend.position = "none"), 
                   nrow = 1, align = c("h", "v"), rel_widths = c(1.12, 1))
power_sd4 <- plot_grid(power, l2, nrow = 2, rel_heights = c(0.9, 0.1))

cbbPalette <-c("#00BA38", "#619CFF")

underlying_lfc_plot <- function(data, y_top){
  ggplot(data = data, aes(LFC, Density, group = Set, color = Set)) +
    geom_line() + xlab("Log Fold Change") + ylab("Density") + 
    geom_area(aes( group = Set, fill = Set), alpha = 0.4, position = "identity") +
    scale_color_manual(values = c("#DC3220", "#005AB5")) +
    scale_fill_manual(values = c("#DC3220", "#005AB5")) +
    #geom_segment(aes(x = log2(1.1), y = y_top, xend = log2(1.1), yend = 0), linetype=2, color =  "#1E88E5", size = 1, alpha = 0.5) + 
    #geom_segment(aes(x = -log2(1.1), y = y_top, xend = -log2(1.1), yend = 0), linetype=2, color =  "#1E88E5", size = 1, alpha = 0.5) +
    #geom_segment(aes(x = log2(1.3), y = y_top, xend = log2(1.3), yend = 0), linetype=2, color = "#FFC107", size = 1, alpha = 0.5) + 
    #geom_segment(aes(x = -log2(1.3), y = y_top, xend = -log2(1.3), yend = 0), linetype=2, color = "#FFC107", size = 1, alpha = 0.5) +
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





SD_in <- 0.6
SD_out <- 0.4
x <- seq(-2,2, by = 0.01)
y_in <- dnorm(x, 0, SD_in)
y_out <- dnorm(x, 0, SD_out)
group <- c(rep("Inset", length(x)), rep("Outset", length(x)))
y <- c(y_in, y_out)
LFC <- c(x,x)
data <- data.frame("Set" = group, "LFC" = LFC, "Density" = y)
SDo_0.4_SDi_0.6_plot <- underlying_lfc_plot(data,  1)

#  Tables

sd <- seq(0.4,.7, by = 0.05)
t1.5_nonDE <-  pnorm(log2(1.5), 0, sd = sd) - pnorm(-log2(1.5), 0, sd = sd)
SD_1_table <- data_frame("Standard deviation" = sd, "Prop. DE at 1.5 threshold" = (1-signif(t1.5_nonDE,2)))

gt_table <- gt(SD_1_table)%>%
  cols_align(align = "center") %>%
  tab_options(row.striping.include_table_body = FALSE) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "#BBBBBB",
      weight = px(1.5),
      style = "solid"
    ),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) 
gt_table <- gt_table %>% tab_style(
  style = list(
    cell_fill(color = '#005AB5', alpha = 0.5)
  ),
  locations = cells_body(
    # not needed if coloring all columns
    rows = c(1)))

gt_table_4 <- gt_table %>% tab_style(
  style = list(
    cell_fill(color = '#DC3220', alpha = 0.5)
  ),
  locations = cells_body(
    # not needed if coloring all columns
    rows = c(5)))


gtsave(gt_table_4, file = "table_SD_4.png")

n <- ggplot() + theme_void()
d <- SDo_0.4_SDi_0.6_plot
img <- readPNG("table_SD_4.png")
g <- rasterGrob(img) 
t <- grid.arrange(g, n, nrow = 2, heights = c(1,.25), top=textGrob("B) DE proportions",gp=gpar(fontsize=22), hjust = 1.5))
top <- grid.arrange(d,t,nrow =1)
figure_1 <- grid.arrange(top, power_sd4, nrow = 2, heights = c(1,1))
ggsave(filename = "fig1.png", figure_1, width = 15, height = 15)
