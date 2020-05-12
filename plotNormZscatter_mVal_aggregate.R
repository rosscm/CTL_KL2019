# Plot drugz normz scatterplots
library(ggplot2)
library(ggrepel)
library(dplyr)

# Run in results directory
setwd("~/GIN/res/20190418_out_drugz_forKeith_mVal")
outDir <- getwd()

# Get data
resF <- list.files(pattern="*_forDrugz_out_.*_vOld_allControl*", path=outDir, full.names=T)
res <- lapply(resF, function(x) read.delim(x, h=T, as.is=T))
dat_names <- sapply(strsplit(basename(resF), "_"), function(x) x[2])

qthres=0.2 # define fdr significance cutoff
dat_all <- list()

for (i in seq_along(resF)) {
  dat <- res[[i]]
  dat <- cbind(ind=1:nrow(dat), dat)

  dat$TP <- NA
  if (i %in% seq(1,length(resF),2)) {
    dat$TP <- "mid"
  } else {
    dat$TP <- "late"
  }

  dat$screen <- NA
  dat$screen <- dat_names[i]

  dat$class <- "Insignificant"
  dat$label <- ""

  sig_supp  <- which(dat$fdr_synth <= qthres)
  sig_synth <- which(dat$fdr_supp <= qthres)

  if (length(sig_supp)==0) {
    dat$class <- "Insignificant"
    dat$label <- ""
  } else {
    dat[sig_supp,]$class <- sprintf("Sensitizer",qthres)
    dat[sig_supp,]$label <- dat[sig_supp,]$GENE
  }
  if (length(sig_synth)==0) {
    dat$class <- "Insignificant"
    dat$label <- ""
  } else {
    dat[sig_synth,]$class <- sprintf("Suppressor",qthres)
    dat[sig_synth,]$label <- dat[sig_synth,]$GENE
  }
  # only plot top 5 sensitizers / suppressors
  dat[6:(nrow(dat)-5),]$label <- ""

  dat_all[[i]] <- dat
}

dat_all <- do.call("rbind", dat_all)
dat_plot <- dat_all
dat_plot$screen <- factor(dat_plot$screen, levels=unique(dat_names))
dat_plot$class <- factor(dat_plot$class, levels=c("Insignificant", "Suppressor", "Sensitizer"))

plotStuff <- function(tp) {
  dat_plot_tp <- filter(dat_plot, TP == tp)

  plot_name <- sprintf("%s/mVal_normzScatter_aggregated_%s_fdr%s.pdf", outDir, tp, qthres*100)
  pdf(plot_name, height=7, width=15.5)
  p <- ggplot(dat_plot_tp, aes(x=ind, y=normZ, colour=class, label=label)) +
          facet_grid(.~screen, scales="free") +
          geom_point(size=1.5, alpha=0.75) +
          geom_text_repel(data=subset(dat_plot_tp, class=="Sensitizer"&label!=""),
                          nudge_x=30, show.legend=F, size=5,
                          direction="y", hjust=0, segment.size=0.2) +
          geom_text_repel(data=subset(dat_plot_tp, class=="Suppressor" & label!=""),
                          nudge_x=-30, show.legend=F, size=5,
                          direction="y", hjust=1, segment.size=0.2) +
          geom_hline(yintercept=0, linetype="dashed", color="black", size=0.7) +
          scale_colour_manual(breaks=c("Suppressor", "Sensitizer"), values=c("lightgrey","#1b9e77","#d95f02")) +
          labs(x="Gene", y="DrugZ normalized z-score") +
          theme_bw() +
          theme(text=element_text(size=13),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                strip.text=element_text(face="bold", size=10),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                legend.position="bottom",
                legend.title=element_blank())
  print(p)
  dev.off()
}

plotStuff(tp="mid")
plotStuff(tp="late")
