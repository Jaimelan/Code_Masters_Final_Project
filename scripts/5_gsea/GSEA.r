# Title: GSEA of differential expressed genes (edger)
# Author: Jaime Llera (modification of an original script by Francisco Garcia)
# Date: 01/06/22

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
R.version.string 

# clean the working space
rm (list = ls ())

# Set the working directory
setwd("~/Results/MA_SexDiff/")


### STEP 1. Differential Expression Analysis for miRNAs
### ===============================================================
### starting
# Library imports:
library (edgeR)
library (mdgsa)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(devtools)
library (mirbaseID) # Cambiar por multiMiR?
library(multiMiR)

# Loading MA Sex Differences genes from a tsv (all.genes.tsv)
mat <- read_tsv(file="~/MA/MA_SexDiff_microRNA/all.genes.tsv")


# MA results
tag <- "SexDiff"   #name of study


# we need gene level p-value and statistic
pvalue    <- as.numeric(mat$pvalue)
statistic <- as.numeric(mat$logFC)
names (pvalue) <- names (statistic) <- mat$gene

# generating a new index using pvalue and statistic from miRNA differential expression
rindex0 <- rindexT <- rindex <- list ()
?pval2index
rindex0[[tag]] <- pval2index (pval = pvalue, sign = statistic)
rindex0[[tag]][1:3]
# number of miRNAs
length(rindex0[[1]])
boxplot(statistic)
boxplot(rindex0[[1]])

# raincloud plot of the data dispersion
df_rindex0 <- as.data.frame(rindex0)
dataset <- "#f4a261"
ggplot(df_rindex0, aes(x = "Dataset", y = SexDiff, fill=dataset)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    color = "#2a9d8f",
    size = 1.3,
    alpha = .05,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

rm(df_rindex0, dataset)


# transfer a ranking index from miRNAs to genes
# We use the get_multimir function to get experimentally validated genes associated
# with all the miRNAs studied

# Function to get a list for each miRNAs targets Gene Symbols
make_targetlist <- function(miRname) {
    
    miRtargets <- (get_multimir(mirna = miRname, org="hsa", table="validated")@data[["target_symbol"]])

    return(miRtargets)
  }


# We get the list of each mirna with its targets, in batches of 500 miRNAs 

# system.time(tar_list_1_500 <- lapply((names(rindex0$SexDiff[1:500])), make_targetlist))
# names(tar_list_1_500) <- names(rindex0$SexDiff[1:500])
# 
# tar_list_501_1000 <- lapply((names(rindex0$SexDiff[501:1000])), make_targetlist)
# names(tar_list_501_1000) <- names(rindex0$SexDiff[501:1000])
# 
# tar_list_1001_1500 <- lapply((names(rindex0$SexDiff[1001:1500])), make_targetlist)
# names(tar_list_1001_1500) <- names(rindex0$SexDiff[1001:1500])
# 
# tar_list_1501_2000 <- lapply((names(rindex0$SexDiff[1501:2000])), make_targetlist)
# names(tar_list_1501_2000) <- names(rindex0$SexDiff[1501:2000])
# 
# tar_list_2000_2498 <- lapply((names(rindex0$SexDiff[2000:2498])), make_targetlist)
# names(tar_list_2000_2498) <- names(rindex0$SexDiff[2000:2498])

# We join the lists together into a tar_list list
tar_list <- c(tar_list_1_500, tar_list_501_1000, tar_list_1001_1500, tar_list_1501_2000, tar_list_2000_2498)
save(tar_list, file = "target_list.RData")

# Loading previously gathered data
load(file="target_list.RData")

# Removing the partial lists
rm(tar_list_1_500, tar_list_501_1000, tar_list_1001_1500, tar_list_1501_2000, tar_list_2000_2498)



# transfer a ranking index from miRNAs to their validated genes
?transferIndex
rindexT[[tag]] <- transferIndex (index = rindex0[[tag]], targets = tar_list, method = "sum") 
rindexT[[tag]][1:3]
boxplot(rindexT[[1]])
# number of genes with transfer
length(rindexT[[1]])
head(rindexT[[1]])


# transform to normal distribution
rindex[[tag]] <- indexTransform (index = rindexT[[tag]], method = "normalize")
rindex[[tag]][1:3]
length(rindex[[1]])
head(rindex[[1]])
boxplot(rindex[[1]])


# exploring different ranking index
t (sapply (rindex0, summary))
t (sapply (rindexT, summary))
t (sapply (rindex, summary))

rindexT[[1]]["PSEN1"]
rindex[[1]]["PSEN1"]

rindexT[[1]]["PSEN2"]
rindex[[1]]["PSEN2"]


### STEP 3. Gene Set Analysis
### ===============================================================

### preparing functional annotation
#  format GO annotation from Ensembl-Biomart
sapply (rindex, class)
sapply (rindex, head)
sapply (rindex, length)
sapply (rindex, summary)
sapply (rindex, function (x) table (x == 0))
# Conservamos los genes ranked
genes <- unique (unlist (lapply (rindex, names)))
length (genes)


# we downloaded the initial annotation from http://www.ensembl.org/biomart/ 
ensgo <- read.table ("mart_export_not_uniq.txt",  header = TRUE, sep = "\t", quote = "", 
                     na.strings = "", colClasses = "character")               
dim (ensgo)
ensgo[1:3,]


# checking and cleaning annotation
table (duplicated (ensgo))
table (ensgo$GO.Term.Accession == "")
table (ensgo$HGNC.symbol       == "")

na.go   <- is.na (ensgo$GO.term.accession)
na.gene <- is.na (ensgo$HGNC.symbol)
table (na.go, na.gene) ## some missing
touse <- !na.go & !na.gene
table (touse)
ensgo <- ensgo[touse,]
dim (ensgo)


# keep just genes in the dataset
touse <- ensgo[,"HGNC.symbol"] %in% genes  ##most of them are present
table (touse)


ensgo <- ensgo[touse,]
ensgo[1:3,]
length (unique (ensgo[,"HGNC.symbol"]))
length (unique (ensgo[,"GO.term.accession"]))

# list format
ensgo <- ensgo[, c("HGNC.symbol", "GO.term.accession")]
ensgo[1:3,]

system.time (annot <- annotMat2list (ensgo))
length (annot)

# propagating ontology
system.time (annot <- propagateGO (annot))
length(annot)



# filtering: better to be done for each dataset
# annot <- annotFilter (annot, minBlockSize = 10, maxBlockSize = 500)
# length (annot)

# split ontologies
annot <- splitOntologies (annot, na.rm = TRUE, verbose = TRUE)
sapply (annot, length)

ontologias <- names (annot)
ontologias

# index
tags <- names (rindex)
tags


### counting GOs by ontology

lons <- list ()
for (onto in ontologias) {
  cat ("\n=============== ", tags, ":", onto, " ===============\n")
  anotacion <- annotFilter (annot[[onto]], index = rindex[[tag]], minBlockSize = 10, maxBlockSize = 300)
  lons[[onto]][tag] <- length (anotacion)
}

as.data.frame (lons)
rowSums (as.data.frame (lons))
sort (unique (rowSums (as.data.frame (lons))))

# We save the session so we can launch the GSEA via SLURM or selecting the number
# of processors to use
save.image(file = "GSEA_Obj.RData")

# Execute the separate script to get the GSE

# ### Gene Set Analysis
# 
# tags <- names (rindex)
# tags
# length(tags)
# 
# for (onto in ontologias) {
#   cat ("\n=============== ", tags, ":", onto, " ===============\n")
#   anotacion <- annotFilter (annot[[onto]], index = rindex[[tag]], minBlockSize = 10,
#                             maxBlockSize = 300)
#   #if (.job$testmode) anotacion <- anotacion[1:3]
#   res <- try (uvGsa (rindex[[tag]], anotacion, fulltable = TRUE))
#   fichero <- paste (tag, "_", onto, ".RData", sep = "")
#   save (list = "res", file = fichero)
# }

load("GSEA_Obj.RData")

### collecting all GSA results

options (width = 170)
corte <- 0.05
ontologias <- c ("bp", "cc", "mf")

# index
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS


# GSA list format
res.gsa.unpa <- list ()
cat ("\n", tag, fill = TRUE, sep = "")
cat ("====================================================================================================", fill = TRUE)
for (onto in ontologias) {
  fichero <- paste (tag, "_", onto, ".RData", sep = "")
  load (fichero)
  res[,"pat"] <- uvPat (res, cutoff = corte)
  res[,"index"] <- pval2index (pval = res[,"pval"], sign = res[,"lor"])
  res[,"Name"] <- getGOnames (res, verbose = FALSE)
  res.gsa.unpa[[onto]][[tag]] <- res
}


res[1:3,]
dim(res)
# are there significant functional terms?
table(res$padj < 0.05)
# what about "transcription coactivator activity"?
res[res$Name == "transcription coactivator activity",]


for (i in 0:length(which(res$padj < 0.05))) {
  cat(rownames(res[which(res$padj < 0.05),])[i])
  cat(" ")
  cat(res$padj[which(res$padj < 0.05)][i])
  cat("\n")
}



# Plots for the top significant terms by ontology

############# BP ############
res_sig_bp_top30 <- res.gsa.unpa[[1]][[1]]
res_sig_bp_top30 <- res_sig_bp_top30[res_sig_bp_top30$pat!=0,]
res_sig_bp_top30 <- head(res_sig_bp_top30[order(res_sig_bp_top30$padj),], n=30)

svg(filename="GO_BP_dotplot.svg", width = 10)
ggplot(res_sig_bp_top30[order(res_sig_bp_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
  geom_point() +
  scale_size_area(max_size = 5) +
  scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 4.53e-21) +
  ggtitle("Diferencias de sexo en la ontología\nProceso biológico") +
  ylab("BP") + xlab("LOR")
dev.off()

############# CC ############
res_sig_cc_top30 <- res.gsa.unpa[[2]][[1]]
res_sig_cc_top30 <- res_sig_cc_top30[res_sig_cc_top30$pat!=0,]
res_sig_cc_top30 <- head(res_sig_cc_top30[order(res_sig_cc_top30$padj),], n=30)

svg(filename="GO_CC_dotplot.svg", width = 10)
ggplot(res_sig_cc_top30[order(res_sig_cc_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
  geom_point() +
  scale_size_area(max_size = 5) +
  scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 1.5e-11) +
  ggtitle("Diferencias de sexo en la ontología\nComponente celular") +
  ylab("CC") + xlab("LOR")
dev.off()


############# MF ############
res_sig_mf_top30 <- res.gsa.unpa[[3]][[1]]
res_sig_mf_top30 <- res_sig_mf_top30[res_sig_mf_top30$pat!=0,]
res_sig_mf_top30 <- head(res_sig_mf_top30[order(res_sig_mf_top30$padj),], n=30)

svg(filename="GO_MF_dotplot.svg", width = 10)
ggplot(res_sig_mf_top30[order(res_sig_mf_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
  geom_point() +
  scale_size_area(max_size = 5) +
  scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 5e-10) +
  ggtitle("Diferencias de sexo en la ontología\nFunción molecular") +
  ylab("MF") + xlab("LOR")
dev.off()


# Plot option 2
ggpubr::ggarrange(
  ggplot(res_sig_bp_top30[order(res_sig_bp_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
    geom_point() +
    scale_size_area(max_size = 5) +
    scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 0) +
    ggtitle("Sex differences in Biological Process Ontology") +
    ylab("BP") + xlab("BP") +
    theme(axis.text.y = element_text(size = 8, angle = 0)),
  
  ggplot(res_sig_cc_top30[order(res_sig_cc_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
    geom_point() +
    scale_size_area(max_size = 5) +
    scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 1.5e-11) +
    ggtitle("Sex differences in Cellular Component Ontology") +
    ylab("CC") + xlab("LOR") +
    theme(axis.text.y = element_text(size = 8, angle = 0)),
  
  ggplot(res_sig_mf_top30[order(res_sig_mf_top30$lor),], aes(x=lor, y=reorder(Name, lor), size=N, color=padj)) + 
    geom_point() +
    scale_size_area(max_size = 5) +
    scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 5e-10) +
    ggtitle("Sex differences in Molecular Function Ontology") +
    ylab("MF") + xlab("LOR") +
    theme(axis.text.y = element_text(size = 8, angle=0)),
  ncol=1 
)



# Plot option 2
# We combine all ontologies in one table
res_sig_bp_top30[,"onto"] <- "bp"
res_sig_cc_top30[,"onto"] <- "cc"
res_sig_mf_top30[,"onto"] <- "mf"


res_sig_byonto <- dplyr::bind_rows(res_sig_bp_top30, res_sig_cc_top30, res_sig_mf_top30)

plt <- ggplot(res_sig_byonto[order(res_sig_byonto$lor),], 
              aes(x=as.factor(pat), y=reorder(Name, lor), size=N, color=lor))

svg(filename="GSEA_dotplot.svg", width = 12, height = 14)
  plt + 
  geom_point() +
  scale_size_area(max_size = 3) +
  scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 0) +
  ggtitle("Sex differences in GO terms") +
  ylab("GSEA") + xlab("LOR") +
  facet_wrap(~ onto, ncol=2, scales = "free_y") +
  scale_x_discrete(labels=c("-1" = "+", "1" = "-")) +
  theme_light() +
  theme(axis.text.y = element_text(size = 6, angle = 0), panel.grid = element_line("transparent"))

dev.off()



# axis.text.y = element_text(size = 6, angle = 0) , panel.background = element_rect(fill='transparent')
ggplot(res_sig_bp_top30[order(res_sig_bp_top30$lor),], aes(x="", y=reorder(Name, lor), size=N, color=lor)) + 
  geom_point() +
  scale_size_area(max_size = 5) +
  scale_colour_gradient2(low="#fde725",mid="#1f9e89", high="#440154", midpoint = 0) +
  ggtitle("Sex differences in Biological Process Ontology") +
  ylab("BP") + xlab("BP") +
  theme(axis.text.y = element_text(size = 8, angle = 0))




### EXIT
warnings ()
sessionInfo ()
q ("no")

