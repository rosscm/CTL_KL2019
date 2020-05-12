# Script to analyze mValidation gene distribution + other functions (tbd)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)

dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/res/%s_out_%s_plots", dataDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

# BAGEL data + query names
dat <- read.delim(sprintf("%s/bftable_all.txt", geneDir))
dat_names <- unique(substr(names(dat)[2:ncol(dat)], 0, nchar(names(dat)[2:ncol(dat)])-4))

# Melt data for use in ggplot2
dat_melt <- melt(dat)
dat_melt$variable <- as.character(dat_melt$variable)
dat_melt$GENE <- as.character(dat_melt$GENE)

df <- dat_melt %>% separate(variable, c("SCREEN", "TP"), "_")

# Remove MC38-QR and B16-QR
df <- filter(df, SCREEN != "MC38.OVA.QR" & SCREEN != "B16.OVA.QR")

plotGuides <- function(gene) {
  #ctrl <- "Intergenic"

  #df_ctrl <- filter(df_all, CLASS == ctrl)
  #mean_ctrl <- ddply(df_ctrl, .(SCREEN), summarise, mean=mean(Tmid_dlogFC))
  #mean_ctrl$SCREEN <- factor(mean_ctrl$SCREEN, levels=dat_names)

  df_plot <- filter(df, GENE == gene)
  df_plot$SCREEN <- factor(df_plot$SCREEN, levels=rev(unique(df$SCREEN)))

  #df_plot$test <- NA
  #for (i in as.character(unique(df_plot$SCREEN))) {
  #  df_plot[which(df_plot$SCREEN == i),]$test <-
  #      df_plot[which(df_plot$SCREEN == i),]$Tmid_dlogFC - mean_ctrl[which(mean_ctrl$SCREEN == i),]$mean
  #}

  #plot_name <- sprintf("%s/%s_%s_guideDistribution_dlog2FC_boxplot.pdf", outDir, gene, ctrl)
  plot_name <- sprintf("%s/%s_guideDistribution_BF_point.pdf", outDir, gene)

  p <- ggplot(df_plot, aes(x=SCREEN, y=value, colour=SCREEN, label=TP)) +
    coord_flip() +
  #  facet_grid(.~TP) +
    geom_point(size=3.5, alpha=0.9) +
    geom_hline(yintercept=5, linetype="dashed", size=0.8) +
    #geom_hline(data=mean_ctrl, aes(yintercept=mean), linetype="dashed", size=1, colour="red") +
    labs(x=NULL, y="Bayes factor score", title=sprintf("Distribution of %s bayes factor score", gene)) +
    theme_bw() +
    geom_text(aes(label=TP), hjust=0.5, vjust=-1) +
    scale_colour_brewer(palette="Dark2") +
    theme(text=element_text(family="sans", size=14),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          plot.title=element_text(hjust=0.5),
          axis.line=element_line(colour="black"),
          #strip.text=element_text(face="bold"),
          legend.position="none",
          legend.title=element_blank())

  pdf(plot_name, height=3.5, width=7, useDingbats=FALSE)
  print(p)
  dev.off()
}

# Plot dropout for specific genes
gene_list <- "Adar"
mapply(plotGuides, gene_list)
