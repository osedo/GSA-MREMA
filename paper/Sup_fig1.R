load("Sup_fig1.RData")
library(ggplot2)
fpr_plot <- ggplot(fpr, aes(x = Method, y = FPR, fill = `N-samples`)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap( ~ Outset, ncol = 1) + theme_minimal()

png(filename = "Sup_fig1.png", width = 700, height = 500)
fpr_plot
dev.off()