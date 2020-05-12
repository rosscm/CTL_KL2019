library(ggplot2)

ggplot(data = df, aes(x = df$CORRELATE, y = df$Spearman, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values= c("blue", "yellow")) +
  xlab("Immune Signatures") + ylab("Spearman Correlation")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14,  family="sans"))

ggsave("coreVSrandom.pdf", plot = last_plot(), 
       width = 18, height = 8, dpi = 300)