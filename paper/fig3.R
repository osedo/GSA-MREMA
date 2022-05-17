load("fig3.RData")
library(ggplot2)


axis_title_size <- 26
axis_text <- 20
legend_text <- 18


cbbPalette <- c( "#619CFF", "#00BA38",    "#D39200",  "#F8766D")
subset_plot <- ggplot(data = data, aes(x = samples, y = proportion, group = Method, color = Method)) +
  geom_line(size = 1.5) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = cbbPalette) + 
  xlab("Number of samples") +
  ylab("Proportion of time gene set is significant") +
  theme_minimal() +
  geom_errorbar(aes(ymin = proportion-ME, ymax = proportion+ME), width=.75, size = 1.5) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        axis.title=element_text(size=axis_title_size),
        legend.key.size = unit(1, 'cm'),
        axis.text = element_text(size = axis_text),
        legend.text = element_text(size = legend_text),
        legend.title = element_text(size = legend_text )) 


png(filename = "fig3.png", width = 1000, height = 750)
subset_plot
dev.off()
