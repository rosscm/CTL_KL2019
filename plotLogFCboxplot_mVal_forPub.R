# Script to analyze mValidation gene distribution + plot
## Use CRAN package librarian to install / update / attach R packages
pkg <- "librarian"[!("librarian" %in% installed.packages()[,"Package"])]
if(length(pkg)) install.packages(pkg)
library(librarian)
shelf(reshape2, xlsx, plyr, dplyr, fgsea, ggplot2, ggpubr, quiet=TRUE)

## USER INPUT ##
dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)

# Output directory
dt <- format(Sys.Date(), "20%y%m%d")
outDir  <- sprintf("%s/res/%s_out_forKeith_mVal_plots", dataDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

# Paths to mouse pathway annotation files (interested in autophagy)
## http://download.baderlab.org/EM_Genesets/April_01_2019/Mouse/symbol/ (GO + Reactome)
## http://download.baderlab.org/EM_Genesets/December_24_2015/Mouse/symbol/Pathways/ (KEGG; last update)
pathGO <- sprintf("%s/anno/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt", dataDir)
pathReact <- sprintf("%s/anno/Mouse_Human_Reactome_April_01_2019_symbol.gmt", dataDir)
pathKegg <- sprintf("%s/anno/Mouse_Human_KEGG_December_24_2015_symbol.gmt", dataDir)

## START ANALYSIS ##
# Foldchange matrices (sheet 2 = guide-level mean; 3 = gene-level mean)
sheet = 3
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

for (i in 1:length(dat_key)) {
  dat_key[[i]]$SCREEN <- NA
  dat_key[[i]]$SCREEN <- dat_names[i]

  # further classifying core genes (suppressor vs. sensitizer)
  dat_key[[i]]$CLASS <- as.character(dat_key[[i]]$CLASS)
  z_melt_dat <- z_melt[z_melt$variable == dat_names[i],]

  for (x in which(dat_key[[i]]$GENE %in% z_melt_dat$GENE)) {
    gene <- as.character(dat_key[[i]][x,]$GENE)
    dat_key[[i]][dat_key[[i]]$GENE==gene,]$CLASS <- z_melt_dat[z_melt_dat$GENE==gene,]$CLASS
  }

  # Calculating differential log2 fold change (dropout - control)
  ## Defining mid vs late timepoints per screen
  ## Setting up rules manually; can be automated, but will take time
  dat_cols <- data.frame(tp1=c("T15.*drop.out.", "T15_Control"),
                         tp2=c("T22.*drop.out", "T24_Control"),
                         tp3=c("T11.*drop.out", "T11_Control"),
                         tp4=c("T12.*drop.out", "T12_Control"))

  if (dat_names[i] == "MC38-OVA") {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp3"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp3"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp4"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp4"], colnames(dat_key[[i]]))]
  } else if (dat_names[i] == "B16-OVA") {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp4"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp4"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp1"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp1"], colnames(dat_key[[i]]))]
  } else {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp1"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp1"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp2"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp2"], colnames(dat_key[[i]]))]
  }

  # Remove control + drop out columns
  cols_del <- grep(paste("*Control","*drop.out*", sep="|"), colnames(dat_key[[i]]))
  dat_key[[i]] <- dat_key[[i]][,-cols_del]
}

# Combine all screen data
dat_all <- do.call("rbind", dat_key)
dat_all <- na.omit(dat_all)

plotStuff <- function(plot_val, anno, pathID, classes, tp, width, height) {

  # Only keep suppressor, sensitizer, targeting control columns
  dat_filt <- dat_all[which(dat_all$CLASS %in% c("Suppressor", "Sensitizer", "Targeting_Controls")),]
  dat_filt[dat_filt$CLASS == "Targeting_Controls",]$CLASS <- "Targeting control"
  dat_filt$CLASS <- factor(dat_filt$CLASS, levels=c("Suppressor", "Sensitizer", "Targeting control"))

  if (plot_val == "PATHWAY") {
    # Annotating genes in autophagy pathway
    path <- gmtPathways(anno)
    autopath <- path[grep(pathID, names(path))]
    autopath2 <- as.character(unlist(autopath))

    dat_filt$PATHWAY <- NA
    dat_filt$PATHWAY <- ifelse(dat_filt$GENE %in% autopath2, "Autophagy", "Other")
    dat_filt$PATHWAY <- ifelse(dat_filt$CLASS == "Targeting control", "Targeting control", dat_filt$PATHWAY)
  }

  # Rename Tmid_dlogFC and Tlate_dlogFC columns
  colnames(dat_filt)[4:5] <- c("mid", "late")
  dat_filt$SCREEN <- factor(dat_filt$SCREEN, levels=dat_names)

  # Melt data by time point and filter by either 'mid' or 'late'
  dat_plot <- melt(dat_filt, measure.vars=c("mid", "late"))
  dat_plot_tp <- filter(dat_plot, variable==tp)

  # Define boxplot colours and title based on whether plotting suppressor / sensitizer
  # vs control (CLASS) or autophagy vs control (PATHWAY)
  if (plot_val == "CLASS") {
    plot_name <- sprintf("%s/mVal_logfcBoxplot_aggregated_%s", outDir, tp)
    cols=c("#F6EB13", "#6D90CA", "#C0BFBF")
    title <- sprintf("Degree of suppressor and sensitizer T-cell dropout vs. control at %s timepoint", tp)
  }

  if (plot_val == "PATHWAY") {
    dat_plot_tp <- filter(dat_plot_tp, PATHWAY != "Other")
    plot_name <- sprintf("%s/mVal_logfcBoxplot_aggregated_%s_autophagy_%s", outDir, tp, deparse(substitute(anno)))
    cols=c("#1B9E77", "#C0BFBF")
    title <- sprintf("Autophagy validation at %s timepoint", tp)
  }

  # Style from boxplot_CTLpaper.R
  p <- ggplot(dat_plot_tp, aes_string(x=plot_val, y="value")) +
          facet_grid(.~SCREEN, scales="free") +
          geom_hline(yintercept=0, linetype="dashed", size=0.5) +
          geom_boxplot(aes_string(fill=plot_val)) +
          labs(y="Differential LFC", title=title, x=NULL) +
          theme_bw() +
          scale_fill_manual(values=cols) +
          theme(text=element_text(family="sans", size=14),
                panel.grid=element_blank(),
                panel.border=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line=element_line(colour="black"),
                #strip.text=element_text(face="bold"),
                legend.position="bottom",
                legend.title=element_blank())

  # Calculate p value
  if (plot_val == "CLASS") {
    p <- p + stat_compare_means(label="p.signif", method="t.test", ref.group="Targeting control")
  }
  if (plot_val == "PATHWAY") {
    p <- p + stat_compare_means(label="p.signif", label.x=1.5, method="t.test")
  }

  # Extract ggplot data
  pg <- ggplot_build(p)
  pg_p <- pg$data[[3]]
  write.table(pg_p, sprintf("%s_pvals.txt", plot_name), col=T, row=F, quote=F, sep="\t")

  pdf(sprintf("%s.pdf", plot_name), width=10, height=5, useDingbats=FALSE)
  print(p)
  dev.off()
}

# Plot class boxplots (sensitizer / suppressor vs control)
#w=10; h=5
plotStuff(plot_val="CLASS", tp="mid")
plotStuff(plot_val="CLASS", tp="late")

# Plot pathway boxplots (autophagy vs control)
plotStuff(plot_val="PATHWAY", anno=pathGO, pathID="0006914", tp="mid", width=w, height=h) # GO, mid
plotStuff(plot_val="PATHWAY", anno=pathReact, pathID="1632852", tp="mid", width=w, height=h) # Reactome, mid
plotStuff(plot_val="PATHWAY", anno=pathKegg, pathID="HSA04140", tp="mid", width=w, height=h) # KEGG, mid
plotStuff(plot_val="PATHWAY", anno=pathGO, pathID="0006914", tp="late", width=w, height=h) # GO, late
plotStuff(plot_val="PATHWAY", anno=pathReact, pathID="1632852", tp="late", width=w, height=h) # Reactome, late
plotStuff(plot_val="PATHWAY", anno=pathKegg, pathID="HSA04140", tp="late", width=w, height=h) # KEGG, late
