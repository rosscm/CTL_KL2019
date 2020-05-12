# Script to analyze mValidation gene distribution + other functions (tbd)
library(reshape2)
library(xlsx)
library(plyr)
library(dplyr)
library(ggplot2)

dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)

dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/res/%s_out_%s_plots", dataDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

# Foldchange matrices (sheet 2 = guide-level mean; 3 = gene-level mean)
sheet = 2
datF <- list.files(pattern="MUS.*xlsx", path=geneDir, full.names=T)
dat <- lapply(datF, function(x) read.xlsx(x, sheetIndex=sheet))

# Get query names
dat_names <- substr(basename(datF), 8, nchar(basename(datF))-40)

# mVal library gene key (define controls)
gene_list <- read.delim(sprintf("%s/mVal_gene_list.txt", geneDir), h=T, as.is=T)

# Restructure gene key list for easy merging with foldchange matrices
gene_key <- melt(gene_list, measure.vars=colnames(gene_list))
colnames(gene_key) <- c("CLASS", "GENE")

# Remove empty rows
gene_key[which(gene_key$GENE == ""),] <- NA
gene_key <- na.omit(gene_key)

# Fix mismatching gene names / symbols (matching to foldchange matrices)
gene_key$GENE[gene_key$GENE == "Luciferase"] <- "luciferase"
gene_key$GENE[gene_key$GENE == "Intergenic"] <- "intergenic"
gene_key$GENE[gene_key$GENE == "CD274"] <- "Cd274"
gene_key$GENE[gene_key$GENE == "CD47"] <- "Cd47"

# Merge foldchange data with gene keys (eg core, control, intergenic, etc)
dat_key <- lapply(dat, function(x) join(x, gene_key, by="GENE"))

# Define timepoint columns
col_names <- data.frame(tp1=c("T15.*drop.out.", "T15_Control"),
                        tp2=c("T22.*drop.out", "T24_Control"),
                        tp3=c("T11.*drop.out", "T11_Control"),
                        tp4=c("T12.*drop.out", "T12_Control"))

df_list <- list()
# Rename timepoint columns for rbind
for (i in 1:length(dat_key)) {
  df <- dat_key[[i]]
  df$SCREEN <- dat_names[i]

  tp3_drop <- grep(col_names[1,"tp3"], colnames(df))
  tp3_ctrl <- grep(col_names[2,"tp3"], colnames(df))
  tp4_drop <- grep(col_names[1,"tp4"], colnames(df))
  tp4_ctrl <- grep(col_names[2,"tp4"], colnames(df))
  tp1_drop <- grep(col_names[1,"tp1"], colnames(df))
  tp1_ctrl <- grep(col_names[2,"tp1"], colnames(df))
  tp2_drop <- grep(col_names[1,"tp2"], colnames(df))
  tp2_ctrl <- grep(col_names[2,"tp2"], colnames(df))

  if (i == 5) {
    df$Tmid_dlogFC <- df[,tp3_drop] - df[,tp3_ctrl]
    df$Tlate_dlogFC <- df[,tp4_drop] - df[,tp4_ctrl]
    df <- df[,-c(tp3_drop, tp3_ctrl, tp4_drop, tp4_ctrl)]

  } else if (i == 6) {
    df$Tmid_dlogFC <- df[,tp4_drop] - df[,tp4_ctrl]
    df$Tlate_dlogFC <- df[,tp1_drop] - df[,tp1_ctrl]
    df <- df[,-c(tp4_drop, tp4_ctrl, tp1_drop, tp1_ctrl)]

  } else {
    df$Tmid_dlogFC <- df[,tp1_drop] - df[,tp1_ctrl]
    df$Tlate_dlogFC <- df[,tp2_drop] - df[,tp2_ctrl]
    df <- df[,-c(tp1_drop, tp1_ctrl, tp2_drop, tp2_ctrl)]

  }
  df_list[[i]] <- df
}

# Combine 6 screen data
df_all <- do.call("rbind", df_list)

plotGuides <- function(gene) {
  #ctrl <- "Targeting_Controls"
  ctrl <- "Intergenic"

  df_ctrl <- filter(df_all, CLASS == ctrl)
  mean_ctrl <- ddply(df_ctrl, .(SCREEN), summarise, mean=mean(Tmid_dlogFC))
  mean_ctrl$SCREEN <- factor(mean_ctrl$SCREEN, levels=dat_names)

  df_plot <- filter(df_all, GENE == gene)
  df_plot$SCREEN <- factor(df_plot$SCREEN, levels=rev(dat_names))

  df_plot$test <- NA
  for (i in as.character(unique(df_plot$SCREEN))) {
    df_plot[which(df_plot$SCREEN == i),]$test <-
        df_plot[which(df_plot$SCREEN == i),]$Tmid_dlogFC - mean_ctrl[which(mean_ctrl$SCREEN == i),]$mean
  }

  plot_name <- sprintf("%s/%s_%s_guideDistribution_dlog2FC_boxplot.pdf", outDir, gene, ctrl)

  p <- ggplot(df_plot, aes(x=SCREEN, y=test, colour=SCREEN)) +
    coord_flip() +
  #  facet_grid(SCREEN~.) +
    geom_hline(yintercept=0, linetype="dashed", size=0.8) +
    #geom_hline(data=mean_ctrl, aes(yintercept=mean), linetype="dashed", size=1, colour="red") +
    geom_point(size=3.5, alpha=0.9) +
    #geom_boxplot() +
    labs(x=NULL, y=expression(atop(paste("Differential log"["2"], " foldchange"), paste("(normalized to intergenic)"))),
        title=sprintf("Distribution of %s Guide Dropout", gene)) +
    theme_bw() +
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

  pdf(plot_name, height=3, width=6, useDingbats=FALSE)
  print(p)
  dev.off()
}

# Plot dropout for specific genes
gene_list <- c("Fitm2", "Jmjd6", "Atg12")
mapply(plotGuides, gene_list)
