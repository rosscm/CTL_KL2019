# Script to analyze mValidation gene distribution + other functions (tbd)
library(reshape2)
library(xlsx)
library(plyr)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpubr)
#library(Cairo)
#library(extrafont)
#loadfonts()

dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)
inDir <- sprintf("%s/res/out_drugz/20190418_out_drugz_forKeith_mVal", dataDir)

dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/res/%s_out_%s_boxplots", dataDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

# Mouse pathway annotation files (interested in autophagy)
pathGO <- sprintf("%s/anno/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt", dataDir)
pathReact <- sprintf("%s/anno/Mouse_Human_Reactome_April_01_2019_symbol.gmt", dataDir)
#pathKegg <- sprintf("%s/anno/Mouse_Human_KEGG_December_24_2015_symbol.gmt", dataDir)

# Drugz normZ (using only targeting controls as non-essentials)
datF <- list.files(pattern="MUS.*vOld_allControl.txt", path=inDir, full.names=T)
dat <- lapply(datF, function(x) read.delim(x, h=T, as.is=T))

# Get query names
dat_names <- sapply(strsplit(basename(datF), "_"), function(x) x[2])

# mVal library gene key + z score data (defines suppressors / sensitizers)
full_lib <- read.delim(sprintf("%s/mVal_full_lib.txt", geneDir), h=T, as.is=T)
gene_list <- read.delim(sprintf("%s/mVal_gene_list.txt", geneDir), h=T, as.is=T)
z_score <- read.delim(sprintf("%s/Core_MeanNormZ.txt", geneDir), h=T, as.is=T)
z_score <- z_score[,-c(1,9)]
names(z_score)[2:7] <- unique(dat_names)[c(1,4,2,3,5,6)]

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
z_melt <- melt(z_score, id.vars="GENE", measure.vars=unique(dat_names))
z_melt$CLASS <- "Non-hit"
z_melt[which(z_melt$value > 0),]$CLASS <- "Suppressor"
z_melt[which(z_melt$value < 0),]$CLASS <- "Sensitizer"

# Merge foldchange data with gene keys (eg core, control, intergenic, etc)
dat_key <- lapply(dat, function(x) join(x, gene_key, by="GENE"))

for (i in 1:length(dat_key)) {
  dat_key[[i]]$SCREEN <- NA
  dat_key[[i]]$SCREEN <- dat_names[i]

  dat_key[[i]]$TP <- NA
  if (i %in% seq(1,length(datF),2)) {
    dat_key[[i]]$TP <- "mid"
  } else {
    dat_key[[i]]$TP <- "late"
  }

  # further classifying core genes (suppressor vs. sensitizer)
  dat_key[[i]]$CLASS <- as.character(dat_key[[i]]$CLASS)
  z_melt_dat <- z_melt[z_melt$variable == dat_names[i],]

  for (x in which(dat_key[[i]]$GENE %in% z_melt_dat$GENE)) {
    gene <- as.character(dat_key[[i]][x,]$GENE)
    dat_key[[i]][dat_key[[i]]$GENE==gene,]$CLASS <- z_melt_dat[z_melt_dat$GENE==gene,]$CLASS
  }
}

dat_all <- do.call("rbind", dat_key)
dat_all <- na.omit(dat_all)

plotStuff <- function(anno, pathID, pathName, classes, tp, plot_val) {

  # annotating genes in pathway
  path <- gmtPathways(anno)
  path2 <- path[grep(pathID, names(path))]
  path2 <- as.character(unlist(path2))

  dat_all$PATHWAY <- NA
  dat_all$PATHWAY <- ifelse(dat_all$GENE %in% path2, pathName, "Other")
  dat_all$PATHWAY <- ifelse(dat_all$CLASS == "Targeting_Controls", "Targeting_Controls", dat_all$PATHWAY)

  if (classes == 4) {
    # Keep all columns but re-arrange
    dat_filt <- filter(dat_all, CLASS != "Non-hit")
    dat_filt[dat_filt$CLASS == "Targeting_Controls",]$CLASS <- "Targeting control"
    dat_filt[dat_filt$PATHWAY == "Targeting_Controls",]$PATHWAY <- "Targeting control"
    dat_filt$CLASS[dat_filt$CLASS == "Intergenic"] <- "Other control"
    dat_filt$CLASS[dat_filt$CLASS == "Non.targeting.controls"] <- "Other control"
    dat_filt$CLASS[dat_filt$CLASS == "Others"] <- "Other control"
    dat_filt$CLASS <- factor(dat_filt$CLASS, levels=c("Suppressor", "Sensitizer", "Targeting control", "Other control"))
  } else if (classes == 3) {
    # Only keep suppressor, sensitizer, targeting control columns
    dat_filt <- dat_all[which(dat_all$CLASS %in% c("Suppressor", "Sensitizer", "Targeting_Controls")),]
    dat_filt[dat_filt$CLASS == "Targeting_Controls",]$CLASS <- "Targeting control"
    dat_filt[dat_filt$PATHWAY == "Targeting_Controls",]$PATHWAY <- "Targeting control"
    dat_filt$CLASS <- factor(dat_filt$CLASS, levels=c("Suppressor", "Sensitizer", "Targeting control"))
  }

  dat_filt$SCREEN <- factor(dat_filt$SCREEN, levels=unique(dat_names))

  if (tp == "mid") {
    dat_plot <- filter(dat_filt, TP == "mid")
  } else if (tp == "late") {
    dat_plot <- filter(dat_filt, TP == "late")
  }

  num_class <- length(unique(dat_plot$CLASS))

  if (plot_val == "PATHWAY") {
    dat_plot <- filter(dat_plot, PATHWAY != "Other")
    plot_name <- sprintf("%s/mVal_normzBoxplot_aggregated_%s_%s_%s_%s.pdf", outDir, tp, deparse(substitute(anno)), pathName, pathID)
    cols=c("#1B9E77", "#C0BFBF")
    title <- sprintf("Depletion of %s at %s timepoint sensitized to T-cell killing", pathName, tp)
    print(table(dat_plot$SCREEN, dat_plot$PATHWAY))
  }

  if (plot_val == "CLASS") {
    plot_name <- sprintf("%s/mVal_normzBoxplot_aggregated_%s_%sclass.pdf", outDir, tp, num_class)
    cols=c("#F6EB13", "#6D90CA", "#C0BFBF")
    title <- sprintf("Degree of suppressor and sensitizer T-cell dropout vs. control at %s timepoint", tp)
  }

  p <- ggplot(dat_plot, aes_string(x=plot_val, y="normZ")) +
        facet_grid(.~SCREEN, scales="free") +
        geom_hline(yintercept=0, linetype="dashed", size=0.5) +
        geom_boxplot(aes_string(fill=plot_val)) +
        labs(y="DrugZ normalized z-score", title=title, x="Genes") +
        theme_bw() +
        scale_fill_manual(values=cols) +
        theme(text=element_text(family="sans", size=14),
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line=element_line(colour="black"),
              legend.position="bottom",
              legend.title=element_blank())

  if (plot_val == "CLASS") {
    p <- p + stat_compare_means(label="p.signif", method="t.test", ref.group="Targeting control")
  }
  if (plot_val == "PATHWAY") {
    p <- p + stat_compare_means(label="p.signif", label.x=1.5, method="t.test")
  }
  ggsave(plot_name, p, width=10, height=4, dpi=300)
}

#plotStuff(anno=pathReact, pathID="1632852", classes=3, tp="mid", plot_val="CLASS")
#plotStuff(anno=pathReact, pathID="1632852", classes=3, tp="late", plot_val="CLASS")
plotStuff(anno=pathReact, pathName="NF-KB", pathID="975138.1", classes=4, tp="mid", plot_val="PATHWAY")
plotStuff(anno=pathReact, pathName="NF-KB", pathID="937041", classes=4, tp="mid", plot_val="PATHWAY")
plotStuff(anno=pathReact, pathName="NF-KB", pathID="168927.2", classes=4, tp="mid", plot_val="PATHWAY")
