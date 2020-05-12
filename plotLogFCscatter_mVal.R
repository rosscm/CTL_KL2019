# Script to analyze mValidation gene distribution + other functions (tbd)
library(reshape2)
library(xlsx)
library(plyr)
library(ggplot2)

dataDir <- "/Users/catherineross/GIN"
geneDir <- sprintf("%s/data/forKeith", dataDir)

dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/res/%s_out_%s", dataDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

# Foldchange matrices (sheet 2 = guide-level mean; 3 = gene-level mean)
sheet = 2
datF <- list.files(pattern="MUS.*xlsx", path=geneDir, full.names=T)
dat <- lapply(datF, function(x) read.xlsx(x, sheetIndex=sheet))

# Get query names
dat_names <- substr(basename(datF), 8, nchar(basename(datF))-40)

# mVal library gene key
full_lib <- read.delim(sprintf("%s/mVal_full_lib.txt", geneDir), h=T, as.is=T)
gene_list <- read.delim(sprintf("%s/mVal_gene_list.txt", geneDir), h=T, as.is=T)
z_score <- read.delim(sprintf("%s/Core_MeanNormZ.txt", geneDir), h=T, as.is=T)
z_score <- z_score[,-c(1,9)]
names(z_score)[2:7] <- dat_names[c(1,4,2,3,5,6)]

# Restructure list for easy merging with foldchange matrices
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

# Merge with z score data
z_melt <- melt(z_score, id.vars="GENE", measure.vars=dat_names)
z_melt$CLASS <- "Core"
z_melt[which(z_melt$value > 0),]$CLASS <- "Suppressor"
z_melt[which(z_melt$value < 0),]$CLASS <- "Sensitizer"

# Merge foldchange data with gene keys (eg core, control, intergenic, etc)
dat_key <- lapply(dat, function(x) join(x, gene_key, by="GENE"))

plot_list <- list()
plot_name_1 <- sprintf("%s/mVal_geneGuideScatter.pdf", outDir)
plot_name_2 <- sprintf("%s/mVal_geneGuideScatter_long.pdf", outDir)

for (i in 1:length(dat_key)) {
  test <- dat_key[[i]]

  # further classifying core genes (suppressor vs. sensitizer)
  test$CLASS <- as.character(test$CLASS)
  z_melt_test <- z_melt[z_melt$variable == dat_names[i],]

  for (b in which(test$GENE %in% z_melt_test$GENE)) {
    gene <- as.character(test[b,]$GENE)
    test[test$GENE==gene,]$CLASS <- z_melt_test[z_melt_test$GENE==gene,]$CLASS
  }

  dat_cols <- data.frame(tp1=c("T15.*drop.out.", "T15_Control"),
                         tp2=c("T22.*drop.out", "T24_Control"),
                         tp3=c("T11.*drop.out", "T11_Control"),
                         tp4=c("T12.*drop.out", "T12_Control"))

  # Setting up rules manually; can be automated, but will take time
  if (i==5) {
    test$Tmid_dlogFC <- test[,grep(dat_cols[1,"tp3"], colnames(test))] - test[,grep(dat_cols[2,"tp3"], colnames(test))]
    test$Tlate_dlogFC <- test[,grep(dat_cols[1,"tp4"], colnames(test))] - test[,grep(dat_cols[2,"tp4"], colnames(test))]
  } else if (i==6) {
    test$Tmid_dlogFC <- test[,grep(dat_cols[1,"tp4"], colnames(test))] - test[,grep(dat_cols[2,"tp4"], colnames(test))]
    test$Tlate_dlogFC <- test[,grep(dat_cols[1,"tp1"], colnames(test))] - test[,grep(dat_cols[2,"tp1"], colnames(test))]
  } else {
    test$Tmid_dlogFC <- test[,grep(dat_cols[1,"tp1"], colnames(test))] - test[,grep(dat_cols[2,"tp1"], colnames(test))]
    test$Tlate_dlogFC <- test[,grep(dat_cols[1,"tp2"], colnames(test))] - test[,grep(dat_cols[2,"tp2"], colnames(test))]
  }

  test$CLASS <- factor(test$CLASS, levels=c("Core", "Suppressor", "Sensitizer",
                                            "Targeting_Controls",
                                            "Non.targeting.controls",
                                            "Intergenic", "Others"))
  # Rank order genes by mean delta logFC
  gene_unique <- as.character(unique(test$GENE))
  rank_mean <- data.frame(GENE=gene_unique, Tmid_meanDlogFC=NA, Tlate_meanDlogFC=NA)

  for (j in seq_along(gene_unique)) {
    rank_mean[j,"Tmid_meanDlogFC"] <- mean(test[test$GENE==gene_unique[j],]$Tmid_dlogFC)
    rank_mean[j,"Tlate_meanDlogFC"] <- mean(test[test$GENE==gene_unique[j],]$Tlate_dlogFC)
  }

  # Print boxplots
  plot_list2 <- list()
  for (k in c("Tmid_dlogFC", "Tlate_dlogFC")) {

    if (k == "Tmid_dlogFC") {
      rank_mid <- rank_mean[order(rank_mean$Tmid_meanDlogFC, decreasing=F),]
      test$GENE <- factor(test$GENE, levels=rank_mid$GENE)
    }
    if (k == "Tlate_dlogFC") {
      rank_late <- rank_mean[order(rank_mean$Tlate_meanDlogFC, decreasing=F),]
      test$GENE <- factor(test$GENE, levels=rank_late$GENE)
    }

    p <- ggplot(test, aes_string(x="GENE", y=k, colour="CLASS")) +
          geom_point() +
          geom_hline(yintercept=0, linetype="dashed", size=1) +
          labs(x="Gene (ranked by mean dlogFC)", y=sprintf("%s per guide\n", k), title=dat_names[i]) +
          theme_classic() +
          scale_colour_brewer(palette="Dark2") +
          theme(text=element_text(size=18),
                axis.text.x=element_text(hjust=1, angle=45))

    plot_list2[[k]] <- p
  }
  plot_list[[i]] <- plot_list2
}

pdf(plot_name_1, height=6.5, width=20)
for (i in seq_along(plot_list)) { print(plot_list[[i]]) }
dev.off()

pdf(plot_name_2, height=6.5, width=100)
for (i in seq_along(plot_list)) { print(plot_list[[i]]) }
dev.off()
