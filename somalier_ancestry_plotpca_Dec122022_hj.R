library(ggplot2)
library(dplyr)

#setwd("/Users/hojin/Desktop/playground")

dat <- read.table("pca.dat.txt", header = T, sep = "\t")
dat <- dat %>%
  mutate(given_ancestry = replace(given_ancestry, given_ancestry == "", NA)) %>%
  mutate(predicted_ancestry = if_else(is.na(given_ancestry) == TRUE, "control", predicted_ancestry))

dat$predicted_ancestry <- factor(dat$predicted_ancestry, levels = c("control", "AFR", "EAS", "SAS", "EUR", "AMR"))

pdf("/path/to/file", height = 6, width = 6)
dat %>% ggplot(aes(x = PC1, y = PC2, group = predicted_ancestry, size = predicted_ancestry)) +
  geom_point(aes(shape=predicted_ancestry, color=predicted_ancestry, alpha = predicted_ancestry), stroke = 1) +
  ggtitle("PCA plot, PC1 VS PC2") + 
  scale_shape_manual(values=c(3, 1, 1, 1, 1, 1, 1)) +
  scale_size_manual(values=c(1, 2, 2, 2, 2, 2, 2)) +
  scale_alpha_manual(values = c(0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) + 
  scale_color_manual(values=c('#FF0000', '#333333','#3300FF', '#FFCC00', '#33FF00', '#9966CC')) +
  theme_bw() +
  theme(aspect.ratio = 1.5,
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()



