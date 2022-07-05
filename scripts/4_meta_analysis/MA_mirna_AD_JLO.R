
######################################
### metagene.R
### 14/06/2020, fgarcia@cipf.es
### 22/09/2020, jfcatala@cipf.es (adapted for MS output)
### 11/12/2020, jfcatala@cipf.es (RNA-Seq with DESeq2)
### 05/05/2022, jlleraoy@gmail.com  (adapted for microRNA)
######################################


## starting
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
R.version.string


library(argparse)
# Obtenemos por argumentos los datasets:

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-d", "--dir", type="character", default="Datasets", 
                    help="Set the directory containing folders for each GSE studied",
                    metavar="dir name")
parser$add_argument("-c", "--comparison", type="character", default="MF", 
                    help="MA between case and control(CC), female(F), male (M) or sex differences(MF)",
                    metavar="comparison for MA")
parser$add_argument("-t", "--tissue", type="character", default="all", 
                    help="Tissue selected (all, bl for blood or br for brain)",
                    metavar="tissues")
parser$add_argument("-p", "--pvalue", type="double", default=0.05, 
                    help="Pvalue selected as threshold for the analysis, defaults at 0.05",
                    metavar="pval")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


## load libraries
library(Biobase)
library(metafor)
library(tibble)
library(limma)


cat("Aportar un directorio como argumento donde se encuentren en carpetas 
      los estudios analizados, con los datos de expresión diferencial como 
      objetos MArrayLM")


#### CARGA DE DATOS ####



setwd(paste("/clinicfs/userhomes/jllera/",args$d, sep=""))
# setwd(paste("/clinicfs/userhomes/jllera/","Datasets", sep=""))


# Opciones para determinar si el análisis será o no con todos los estudios y si será de diferencias de sexo
sexdiff <- args$c #  -c
mixed <- args$t # Tissue Puede ser all, br (Brain) o bl (Blood) o 120 (menos el GSE120584)
corte <- args$p # Pvalue selected for sig genes

# Ajuste del nombre que tendra el fichero en base a los comandos
if (sexdiff == "CC") {
    comp <- "AD_NC_"
} else if (sexdiff == "MF"){
    comp <- "SexDiff_"
} else if (sexdiff == "F"){
    comp <- "Female_"
} else if (sexdiff == "M"){
    comp <- "Male_"
} else {stop("Specify M(male), F (Female), MF(Female v Male) or Case v Control (CC)")}

if (mixed == "all") {
    tissue <- ""
} else if (mixed == "br") {
    tissue <- "Brain_"
} else if (mixed == "bl") {
    tissue <- "Blood_"
} else if (mixed == "120") {
    tissue <- "Minus120_"
} else {stop("Specify 'br'(brain), 'bl'(blood) or 'all'")}


## load data
GSEnams <- list.files(pattern="GSE[0-9]{5,6}$")
DElist_SexDiff <- DElist_NC_v_AD <- DElist_FNC_v_FAD <- DElist_MNC_v_MAD <- list()

# Loading MArrayLM objects into a list
for (i in 1:length(GSEnams)){
  print(paste("Cargando datos del estudio ", GSEnams[i]))
  
  if (sexdiff == "CC") {
    DElist_NC_v_AD[[GSEnams[i]]] <- readRDS(paste(GSEnams[i],"/",GSEnams[i],"_ED_Case_Ctrl.RDS", sep=""))
  } else if (sexdiff == "MF"){
    DElist_SexDiff[[GSEnams[i]]] <- readRDS(paste(GSEnams[i],"/",GSEnams[i],"_ED_SexDiff.RDS", sep=""))
  } else if (sexdiff == "F") {
    DElist_FNC_v_FAD[[GSEnams[i]]] <- readRDS(paste(GSEnams[i],"/",GSEnams[i],"_ED_FCase_FCtrl.RDS", sep=""))
  } else if (sexdiff == "M") {
    DElist_MNC_v_MAD[[GSEnams[i]]] <- readRDS(paste(GSEnams[i],"/",GSEnams[i],"_ED_MCase_MCtrl.RDS", sep=""))
  }
}


# TODO Hacer que la seleccion del tipo de metaanalisis sea un poco mas elegante
gseBrain <- c("GSE157239","GSE16759","GSE48552")
gseBlood <- c("GSE120584", "GSE46579")
gseMinus120 <- c("GSE157239","GSE16759","GSE48552", "GSE46579")



## directories
wd <- "/clinicfs/userhomes/jllera/MA/" #Output dir

metaanalysis_name <- paste("MA_", comp, tissue, "microRNA", sep="")
metaanalysis_name


wd <- paste0(wd,metaanalysis_name,"")
dir.create(wd, recursive = TRUE)
setwd(wd)

# STEP 0. Pre-processing previous data
# ===============================================================


## Load all ED results in a list

if (sexdiff == "CC"){
  if (mixed == "all") {
  EDs <- DElist_NC_v_AD
  } else if (mixed == "br") {
    EDs <- DElist_NC_v_AD[gseBrain]
  } else if (mixed == "120"){
    EDs <- DElist_NC_v_AD[gseMinus120]
  } else {
    EDs <- DElist_NC_v_AD[gseBlood]
  }
} else if (sexdiff == "MF"){
  if (mixed == "all") {
    EDs <- DElist_SexDiff
  } else if (mixed == "br") {
    EDs <- DElist_SexDiff[gseBrain]
  } else if (mixed == "120"){
      EDs <- DElist_SexDiff[gseMinus120]
    } else {
    EDs <- DElist_SexDiff[gseBlood]
  }
} else if (sexdiff == "F") {
  if (mixed == "all") {
    EDs <- DElist_FNC_v_FAD
  } else if (mixed == "br") {
    EDs <- DElist_FNC_v_FAD[gseBrain]
  } else if (mixed == "120") {
    EDs <- DElist_FNC_v_FAD[gseMinus120]
  } else {
    EDs <- DElist_FNC_v_FAD[gseBlood]
    }
} else if (sexdiff == "M") {
  if (mixed == "all") {
    EDs <- DElist_MNC_v_MAD
  } else if (mixed == "br") {
    EDs <- DElist_MNC_v_MAD[gseBrain]
  } else if (mixed == "120") {
    EDs <- DElist_MNC_v_MAD[gseMinus120]
  } else {
    EDs <- DElist_MNC_v_MAD[gseBlood]
  }
}



rm(list=setdiff(ls(), c("EDs","wd","metaanalysis_name", "mixed", "gseBrain", "gseBlood", "corte"))) 
#Remove unnecessary variables

## load desired comparison
## may be NC vs AD comparison or Female AD vs Male AD comparison

## logFC limma vs log2FC DESeq2


### DESeq2
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
# In the most recent versions of DESeq2 (≥ 1.16.0), the shrinkage of LFC estimates is not performed by default.
# This means that the log2 foldchanges would be the same as those calculated by:
# log2 (normalized_counts_group1 / normalized_counts_group2)

### limma
# https://www.biostars.org/p/100460/#129719
# Limma's "Log(FC)" = mean(log2(Group1)) - mean(log2(Group2))

## Calculate SE
SE_array <- function(fit) {
  if ("lfcSE" %in% colnames(fit)) { #Si los datos son un dataframe de ED con DESeq2
    res <- fit
  } else { #Datos de limma
  #OPTION1: https://support.bioconductor.org/p/70175/ by Gordon Smyth/January Weiner
  #The effect sizes are contained in fit$coefficients
  summary(fit$coefficients)
  head(fit$coefficients)
  #The standard errors can be obtained from 2 sub-options:
  # SE.coef <- sqrt(fit$s2.post) #JANUARY  (Here I have a problem when having several contrasts
  #                                        #because I have the same information for all contrasts)
  # head(SE.coef)
  # summary(SE.coef)
  SE.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled #GORDON
  head(SE.coef)
  summary(SE.coef)

  #OPTION2: https://support.bioconductor.org/p/70175/ by  Steve Lianoglou  (SE PARECE A GORDON)
  allgenes <- topTable(fit, number = "all", confint=TRUE, adjust.method = "fdr")
  dim(allgenes)
  allgenes[, "SE"] <- (allgenes[, "CI.R"] - allgenes[, "CI.L"])/ 3.92
  head(allgenes)

  #final results
  table(rownames(SE.coef) == rownames(fit$coefficients))
  mat <- cbind(fit$coefficients, SE.coef)
  colnames(mat) <- c("coef", "se.coef")
  head(mat)

  int <- intersect(rownames(allgenes), rownames(mat))
  length(int)
  # res <- cbind(allgenes, mat[rownames(allgenes),])
  res <- cbind(allgenes, mat[match(rownames(allgenes),rownames(mat)),])
  head(res)
  dim(res)
  }

  return(res)
}

EDs_sel1 <- lapply(EDs, SE_array)


# STEP 1. Preparing input for meta-analysis: LOR and SE matrix
# ===============================================================

# we search a list including all unique ID genes for all studies
genes <- NULL
for (fi in EDs_sel1){
  genes <- c(genes, rownames(fi))
}

length (genes)
genes <- unique (genes)
length (genes)
genes <- sort (genes)


### generating matrix with all logFC for all studies
mat.logFC <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel1))
rownames (mat.logFC) <- genes
colnames (mat.logFC) <- gsub("_ED", "", names(EDs_sel1))
head (mat.logFC)

for (i in 1:length(EDs_sel1)){
  co <- gsub("_ED", "", names(EDs_sel1[i]))
  res <- EDs_sel1[[i]]

  if ("log2FoldChange" %in% colnames(res)) { # DESeq2
    logFC <- res$log2FoldChange
  } else { # limma
    logFC <- res$logFC
  }

  names (logFC) <- (rownames (res))
  mat.logFC[, co] <- logFC[rownames(mat.logFC)]
}

head (mat.logFC)
tail(mat.logFC)
table (is.na(mat.logFC))
dim (mat.logFC)

# select genes included at least in 2 or more studies
mat.logFC.NA <- is.na(mat.logFC)
head(mat.logFC.NA)
sum.NA <-  apply(mat.logFC.NA, 1, sum)
table(sum.NA)
min.sum.NA <- sum.NA < ncol(mat.logFC) - 1
table(min.sum.NA)

# filter by min.sum.NA
mat.logFC <- mat.logFC[min.sum.NA == T, ]
dim(mat.logFC)

# filter feats that does not appear in both tissues studied
if(mixed == "all"){
  not_all <- apply(mat.logFC, 1, function(x) {sum(is.na(x[gseBlood])) == length(gseBlood) | sum(is.na(x[gseBrain])) == length(gseBrain)}) #No son de los dos tejidos
  mat.logFC <- mat.logFC[!not_all, ]
  dim(mat.logFC)
}


### generating matrix with all SE for all studies
mat.SE <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel1))
rownames (mat.SE) <- genes
colnames (mat.SE) <- names(EDs_sel1)
head (mat.SE)


# (SE FROM GORDON: se.coef)
for (i in 1:length(EDs_sel1)){
  co <- names(EDs_sel1[i])
  res <- EDs_sel1[[i]]

  if ("lfcSE" %in% colnames(res)) { # DESeq2
    SE <- res$lfcSE
  } else { # limma
    SE <- res$se.coef
  }

  names (SE) <- (rownames (res))
  mat.SE[, co] <- SE[rownames(mat.SE)]
}


head (mat.SE)
tail(mat.SE)
table (is.na(mat.SE))
dim (mat.SE)

# filter by min.sum.NA
mat.SE <- mat.SE[min.sum.NA == T, ]
dim(mat.SE)

if(mixed == "all"){
  mat.SE <- mat.SE[!not_all, ]
  dim(mat.SE)
}



# STEP 2. Meta-analysis for genes
# ===============================================================

# suppose between-study variance is non-zero.
# there are different methods to estimate this variance:
# DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
# Now we have logFC and SE  (not VARIANCE), so:
# yi -> logFC   sei -> SE
# result.lor <- rma(yi = mat.logFC[1, ],
#                   sei = mat.SE[1, ],   #pay attention, not vi (varianze)
#                   method = "DL") # DerSimonian-Laird.


# explore the function to do the meta-analysis
#?rma

MA <- lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SE[x, ],
                             method = "DL")})

# MA <- lapply(1:length(rownames(mat.logFC)),
#              function(x){rma(yi = mat.logFC[x, ],
#                              sei = mat.SE[x, ],
#                              method = "FE")})

names (MA) <- rownames(mat.logFC)
class (MA)
length(MA)
head (MA)
MA[[1]]

#result.logFC$pval      #p-value about logFC = 0
#result.logFC$ci.lb     #IC down
#result.logFC$ci.ub     #IC up
#result.logFC$b         #estimation of combined logFC

#data.frame including all detailed results:
result_meta <- as.data.frame(do.call("rbind",
                                     lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub,
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                            })))

colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                           "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")

p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr")
p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
logFDR <- -log10(p.adjust.fdr)
result_meta <-cbind(result_meta, p.adjust.fdr, p.adjust.BY,logFDR)
#result_meta <- cbind(result_meta, p.adjust.fdr, p.adjust.BY)
head(result_meta)


# significant genes
corte # En el caso de SexDiffs no sale significativos con .05, lo triplicamos
table(result_meta[, "pvalue"] < corte)
table(result_meta[, "p.adjust.fdr"] < corte)
table(result_meta[, "p.adjust.BY"] < corte)


# add number of studies where the gene is evaluated
n.studies <-  ncol(mat.logFC) - sum.NA
table(n.studies)
n.studies <- n.studies [rownames(mat.logFC)]
length(n.studies)
result_meta[, "n.studies"]  <- n.studies
head(result_meta)
summary(result_meta$p.adjust.fdr)


sig.genes.df = result_meta[result_meta$p.adjust.fdr < corte,]
dim(sig.genes.df)

write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),] %>% rownames_to_column('gene'), file = "sig.genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = result_meta[order(result_meta$p.adjust.fdr),] %>% rownames_to_column('gene'), file = "all.genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS
# ===============================================================

#add 4 new variables about influence & sensitivity analysis:

for (i in rownames(sig.genes.df)){
  #print(i)
  #define studies for each function (not NA)
  estudios <- colnames(mat.logFC)[!mat.logFC.NA[i,]]

  #influence info 1:
  #number of studies where the sign of the logOR is the same  of the global logOR:
  sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios)))

  #influence info 2: how many studies could be influencers?
  inf <- influence(MA[[i]])
  res <- paste(estudios[inf$is.infl], collapse = ",")
  sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)

  #sensivity analysis
  l1 <-as.data.frame(leave1out(MA[[i]]))
  rownames(l1) <- estudios

  #1. p.value about differences between all estimates from leave one out
  #   and global estimate)
  sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                           mu=as.numeric(MA[[i]]$b))$p.value
  #2. number of  studies where pvalue > 0.05
  # (we hope p-values < 0.05, significant estimates)
  res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")
  sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}


## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:

#1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
table(sig.genes.df$infl.same.sign.logFC)

#2. INFLUENCE STUDIES. How many functions including influence studies?
table(sig.genes.df$infl.nstudies=="non")

#3. SENSITIVITY. In global, are there many functions with differences in the estimate?
table(sig.genes.df$sensi.global < 0.05)

#4. SENSITIVITY.  How many functions including changes in the significance about
# its new estimate  after leave1out?
table(sig.genes.df$sensi.specific == "all.p.values < 0.05")


#save final results:
cat ("ID\t", file = "sig.genes.df.txt")
write.table(sig.genes.df, file = "sig.genes.df.txt", sep ="\t", quote = F,
            append = TRUE, row.names = T)




# STEP 4. Visualization of significant genes
# ===============================================================

#select significant functions to visualize:
sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < corte,]

sig.results
dim(sig.results)

dir.create(paste0(wd,"/","plots"), recursive = TRUE)
setwd(paste0(wd,"/","plots"))

selMethod <- "DL"

for (i in 1:nrow(sig.results)){
  mygenes <- rownames(sig.results)[i]
  res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")

  #FOREST PLOT
  # png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 ,
  #      height = 960, res = 200)
  png (filename = paste(gsub("-","_",mygenes),"_FOREST",".png", sep =""), width = 960 ,
       height = 960, res = 200)
  forest(res,
         slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
         xlab="logFC", cex=0.7,
         mlab=paste(selMethod, "Model for All Studies", sep = " "),
         border = "black", #Color del borde del rombo
         col = "red", #Color del rombo
         main = paste("\n", mygenes, sep=""))
  text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
  dev.off()

  #FUNNEL PLOT
  png (filename = paste(gsub("-","_",mygenes),"_FUNNEL", ".png", sep =""), width = 960 ,
       height = 960, res = 200)
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance",
         back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
  par(mfrow=c(1,1))
  dev.off()

  #INFLUENCE PLOTS
  # That shows various diagnostic measures
  png (filename = paste(gsub("-","_",mygenes), "_INFLUENCE", ".png", sep =""), width = 960 ,
       height = 960, res = 200) ##CAMBIAR
  inf <- influence(res)
  #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
  plot(inf)
  dev.off()

}
print("Generating report")
# STEP 5. Generating report
# ===============================================================
sig.genes.df <- sig.genes.df[order(sig.genes.df$p.adjust.fdr),]
sig.genes.df.5 <- head(sig.genes.df, n=5)

# Summary table
res <- function(df){
  total <- nrow(df)
  if("padj" %in% colnames(df)){
    df$adj.P.Val <- df$padj
  }
  sig <- sum(df$adj.P.Val < 0.05)
  res <- c(total, sig)
  names(res) <- c("genes", "sig.genes")
  return(res)
}

tabla_resumen <- function(x) {
  table <- lapply(x, res)
  # Las uno en una tabla
  table <- dplyr::bind_rows(table, .id = "Study")

    return(table)
}

#Tabla resumen con los genes y los genes significativos (p.val.adj < 0.05)
table_sum <- tabla_resumen(EDs_sel1)
# table_sum$Study <- gsub("_ED","", table_sum$Study)
# if(all(names(apply(!is.na(mat.logFC), 2, sum)) == table_sum$Study)){
#   table_sum <- dplyr::bind_cols(table_sum, apply(!is.na(mat.logFC), 2, sum))
#   table_sum <- table_sum %>% dplyr::relocate(...4, .after = genes)
# }

# Function to create multiple tabs
make.tabs <- function(sig.genes.df){
  res <- NULL
  for(g in rownames(sig.genes.df)){
    file_name <- gsub("-","_", g)
    res <- c(res, '### ', g, ' {-} \n',
             "**Statistics of ", g, " meta-analisys** \n",
             "```{r, fig.align='center'}", '\n',
             "kable(sig.genes.df['",g,"',1:11])", '\n',
             '```', '\n\n',
             "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g, ") \n\n",
             "**Forest plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', wd, '/plots/', file_name, '_FOREST.png")\n',
             '```', '\n',
             "**Funnel plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', wd, '/plots/', file_name, '_FUNNEL.png")\n',
             '```', '\n',
             "**Incluence plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', wd, '/plots/', file_name, '_INFLUENCE.png")\n',
             '```', '\n\n')
  }
  return(res)
}


# Generating volcano plotly
sig_limit <- 0.05
lfc_limit <- 1.5

data <- result_meta
data$gene_name <- rownames(data)
rownames(data) <- NULL
data <- data[,c("gene_name", "logFC", "p.adjust.fdr")]

# add a grouping column; default value is "not significant"
data["Group"] <- "NotSignificant"

# change the grouping for the entries with significance but not a large enough Fold change
data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) < lfc_limit ),"Group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
data[which(data['p.adjust.fdr'] > sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("FoldChange > ",lfc_limit," | < -",lfc_limit)

# change the grouping for the entries with both significance and large enough fold change
data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("Significant & FoldChange > ",lfc_limit," | < -",lfc_limit)

data$texto <- paste0("Gene: ", data$gene_name, "\np.val.fdr: ", data$p.adjust.fdr)

library(ggplot2)
volcano_plot <- ggplot(data, aes(x = logFC, y = logFDR)) +
  geom_point(aes(color = Group, text = texto)) +
  #scale_color_manual(values = c("firebrick1", "firebrick4", "royalblue3", "seagreen3")) +
  scale_color_manual(values = c("#E69F00", "#DB5E00", "#0072B2", "#009E73")) + #https://stackoverflow.com/q/57153428
  xlab("logFC") + ylab("-log10(FDR)")

# Create the Rmd to knit
ns <- nrow(sig.genes.df)
cat(
  '---
title: "Global Meta-analysis"
output:
  html_document:
    toc: false
    toc_float: false
    code_folding: hide
    number_sections: true
    theme: spacelab
---
## ', metaanalysis_name, ' {.tabset .tabset-pills -}

```{r, warning=F, message=F}
library(dplyr)
library(knitr)
library(DT)
library(ggplot2)
library(plotly)
load("', wd, '/', metaanalysis_name, '.RData")
```  \n
### MA results (significant genes) {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "datatable(sig.genes.df[,c(1:4,11,14)], caption='Statistics of ", ns, " significant genes', escape = FALSE, filter = 'top')", '\n',
  '```', '\n\n',
  '### MA results (volcano plot) {-} \n',
  "```{r, fig.align='center', out.width = '99%'}", '\n',
  "ggplotly(volcano_plot)", '\n',
  '```', '\n\n',
  '### MA results (all genes) {-}  \n',
  'The results of the meta-analysis of all the genes can be found in the file ["all.genes.tsv"](./all.genes.tsv).',
  # "```{r, fig.align='center'}", '\n',
  # "datatable(result_meta[order(result_meta$p.adjust.fdr),c(1:4,11,14)], caption='Statistics of all genes', escape = FALSE, filter = 'top')", '\n',
  # '```',
  '\n\n',
  '### Individual results (summary) {-} \n',
  "```{r, fig.align='center'}", '\n',
  "kable(table_sum, col.names = c('Study','All genes','Significant genes'))", '\n',
  '```',
  '\n\n',
  '### sessionInfo {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "date()", "\n",
  "sessionInfo()", '\n',
  '```', '\n\n',
  sep = "",
  file = paste0(wd, "/", metaanalysis_name, ".Rmd"))

cat(
'---
title: "Meta-analysis of genes"
subtitle: "Draft, ', format(Sys.time(), "%d-%m-%Y"),'v6c"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    code_folding: hide
    number_sections: false
    theme: spacelab
---
## ', metaanalysis_name, ' {- .tabset .tabset-fade .tabset-pills}

```{r, warning=F, message=F}
library(dplyr)
library(knitr)
library(DT)
load("', wd, '/', metaanalysis_name, '.RData")
```  \n',
make.tabs(sig.genes.df.5[order(row.names(sig.genes.df.5)),]), "\n\n",
'### sessionInfo {-}  \n',
"```{r, fig.align='center'}", '\n',
"date()", "\n",
"sessionInfo()", '\n',
'```', '\n\n',
  sep = "",
  file = paste0(wd, "/", metaanalysis_name, "_genes.Rmd"))

save(sig.genes.df, result_meta, table_sum, volcano_plot, file = paste0(wd, "/", metaanalysis_name, ".RData"))
# Render the Rmd created into html here
rmarkdown::render(paste0(wd, "/", metaanalysis_name, ".Rmd"))
rmarkdown::render(paste0(wd, "/", metaanalysis_name, "_genes.Rmd"))



rm(list=ls())