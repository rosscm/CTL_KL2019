# Script to analyze mValidation gene distribution + other functions (tbd)
library(reshape2)
library(xlsx)
library(plyr)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpubr)

dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)

# Write plots to results directory
dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/res/%s_out_%s_plots", dataDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

# Mouse pathway annotation files (interested in autophagy)
pathGO <- sprintf("%s/anno/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt", dataDir)
pathReact <- sprintf("%s/anno/Mouse_Human_Reactome_April_01_2019_symbol.gmt", dataDir)

# Foldchange matrices (sheet 2 = guide-level mean; 3 = gene-level mean)
sheet = 2
datF <- list.files(pattern="MUS.*xlsx", path=geneDir, full.names=T)
dat <- lapply(datF, function(x) read.xlsx(x, sheetIndex=sheet))

# Get query names
dat_names <- substr(basename(datF), 8, nchar(basename(datF))-40)

# mVal library gene key + z score data (defines suppressors / sensitizers)
full_lib <- read.delim(sprintf("%s/mVal_full_lib.txt", geneDir), h=T, as.is=T)
gene_list <- read.delim(sprintf("%s/mVal_gene_list.txt", geneDir), h=T, as.is=T)
z_score <- read.delim(sprintf("%s/Core_MeanNormZ.txt", geneDir), h=T, as.is=T)
z_score <- z_score[,-c(1,9)]
names(z_score)[2:7] <- dat_names[c(1,4,2,3,5,6)]

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

# Merge with z score data (z>0 = supressor; z<0 = sensitizer)
z_melt <- melt(z_score, id.vars="GENE", measure.vars=dat_names)
z_melt$CLASS <- "Non-hit"
z_melt[which(z_melt$value > 0),]$CLASS <- "Suppressor"
z_melt[which(z_melt$value < 0),]$CLASS <- "Sensitizer"

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

plotStuff <- function(anno, pathID, pathName, tp) {

  # annotating genes in pathway
  path <- gmtPathways(anno)
  df_all$PATHWAY <- "Other"

  for (i in seq_along(pathID)) {
    # Get pathway genes
    path2 <- path[grep(pathID[i], names(path))]
    path2 <- as.character(unlist(path2))
    # Indicate pathway genes in data
    df_all[which(df_all$GENE %in% path2),]$PATHWAY <- pathName[i]
    #df_all$PATHWAY <- ifelse(df_all$CLASS == "Targeting_Controls", "Targeting_Controls", df_all$PATHWAY)
    df_all$PATHWAY <- ifelse(df_all$CLASS == "Intergenic", "Intergenic", df_all$PATHWAY)
  }

  # Rename Tmid_dlogFC and Tlate_dlogFC columns
  colnames(df_all)[6:7] <- c("mid", "late")

  # Set factor levels
  df_all$SCREEN <- factor(df_all$SCREEN, levels=dat_names)
  df_plot <- melt(df_all, id.vars=c("GENE", "CLASS", "SCREEN", "PATHWAY"), measure.vars=c("mid", "late"))

  # Filter by time point
  df_plot_tp <- filter(df_plot, variable==tp)

  # Get intergenic control data
  df_ctrl <- filter(df_plot_tp, CLASS == "Intergenic")
  mean_ctrl <- ddply(df_ctrl, .(SCREEN), summarise, mean=mean(value))
  mean_ctrl$SCREEN <- factor(mean_ctrl$SCREEN, levels=dat_names)

  # Normalize gene fold change values to intergenic
  df_plot_tp$test <- NA
  for (i in as.character(unique(df_plot_tp$SCREEN))) {
    df_plot_tp[which(df_plot_tp$SCREEN == i),]$test <-
        df_plot_tp[which(df_plot_tp$SCREEN == i),]$value - mean_ctrl[which(mean_ctrl$SCREEN == i),]$mean
  }

  #df_final <- filter(df_plot_tp, PATHWAY != "Other" & PATHWAY != "Targeting_Controls")
  df_final <- filter(df_plot_tp, PATHWAY != "Other")
  df_final$PATHWAY <- factor(df_final$PATHWAY, levels=c("Autophagy", "NF-KB", "Intergenic"))

  # Plot parameters
  plot_name <- sprintf("%s/mVal_logFcBoxplot_%s_%s_%s", outDir, tp, deparse(substitute(anno)),
                      paste(c(pathName, pathID), collapse="_"))
  title <- sprintf("Depletion of %s at %s timepoint\nsensitized to T-cell killing",
                    paste(pathName, collapse=" and "), tp)
  cols <- c("#66C2A5", "#FC8D62", "#C0BFBF")
  print(table(df_final$SCREEN, df_final$PATHWAY))

  # Style from boxplot_CTLpaper.R
  p <- ggplot(df_final, aes(x=PATHWAY, y=test)) +
          facet_grid(.~SCREEN, scales="free") +
          geom_hline(yintercept=0, linetype="dashed", size=0.5) +
          geom_boxplot(aes(fill=PATHWAY)) +
          labs(y="Differential LFC\n(normalized to intergenic)",
              title=title, x=NULL) +
          theme_bw() +
          scale_fill_manual(values=cols) +
          #scale_fill_brewer(palette="Set2") +
          theme(text=element_text(family="sans", size=14),
                panel.grid=element_blank(),
                panel.border=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line=element_line(colour="black"),
                legend.position="bottom",
                legend.title=element_blank())

  p <- p + stat_compare_means(label="p.signif", method="t.test", ref.group="Intergenic")

  # Extract ggplot data
  pg <- ggplot_build(p)
  pg_p <- pg$data[[3]]
  write.table(pg_p, sprintf("%s_pvals.txt", plot_name), col=T, row=F, quote=F, sep="\t")

  pdf(sprintf("%s.pdf", plot_name), width=8, height=5, useDingbats=FALSE)
  print(p)
  dev.off()
}

plotStuff(anno=pathReact, pathName=c("Autophagy", "NF-KB"), pathID=c("1632852", "975138.1"), tp="mid")
