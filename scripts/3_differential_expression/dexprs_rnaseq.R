#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: Differential expression for RNA-Seq with limma Voom
##
## Purpose of script: Plots information and basic statistics of a given data matrix (passed as an argument)
##
## Author: Jaime Llera
##
## Date Created: 05-04-2022
##
## Email: jlleraoy@gmail.com
##
## ---------------------------
##

library(edgeR)


#### CARGA DE DATOS ####
# Obtenemos por argumentos los datasets:
args = commandArgs(trailingOnly = T)

# Evaluamos que en los argumentos se aporte una carpeta y una matriz de datos
if (length(args)!=1) {
  stop("Only 1 must be supplied: 1. experiment_id", call.=FALSE)
} else if (length(args)==1) {
  print("Realizando analisis de expresion diferencial...")
}

# elegimos la carpeta de trabajo, cargamos los datos
setwd(args[1])

pDataname <- paste0(args[1], "_pData.RDS", sep="")
pDa <- readRDS(pDataname)

matname <- paste0(args[1], "_mat.RDS", sep="")
datamatrix <- readRDS(matname)

# deshacemos la transformacion logaritmica
datamatrix <- as.matrix(apply(datamatrix, 1:2, function(x) {2^x-1}))


####-------------------ANALISIS DE EXPRESION DIFERENCIAL--------------------####
# Cargamos la matriz de conteos como un objeto DGEList que permita realizar la normalizacion con edgeR
d0 <- DGEList(datamatrix)
d0 <- calcNormFactors(d0, method="TMM")


# Se eliminan genes poco expresados (igual esto hay que hacerlo antes de guardar las matrices de datos
# para que la representacion sea adecuada, en el caso de los RAW data habra que hacerlo, preguntar esto)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) <= cutoff) # Opcion de Jose
# drop <- rowSums(edgeR::cpm(d0$counts)>1) <= cutoff # Opcion de Irene
d <- d0[-drop,] 
dim(d) # number of genes left


casos <- pDa$casos
sex <- pDa$sex
group <- interaction(casos, sex)

# Plot de escalado multidimensional
plotMDS(d, col = as.numeric(group))



# Funcion que realiza el ajuste con limma para realizar el contraste
dif_exp <- function(GSEmat,        # Matriz de datos
                    contraste,     # Contraste
                    grupos = NULL, # Variable donde están los elementos del contraste
                    lotes = NULL,  # Si lo hay, variable para el efecto batch
                    trend = FALSE) {
  
  grupos <- as.factor(grupos) #Los grupos tienen que ser un factor
  
  if (is.null(lotes)) {
    design <- model.matrix(~0+grupos)
    colnames(design) <- levels(grupos)
  } else {
    lotes <- as.factor(lotes)
    design <- model.matrix(~0+grupos+lotes)
    colnames(design) <- c(levels(grupos), levels(lotes)[-2])
    # colnames(design) <- c(levels(grupos), levels(lotes))
  }
  
  cont.matrix <- makeContrasts(contrasts = contraste, levels = design)
  
  set.seed(123)
  
  # Aplicacion del ajuste con limma::voom
  # La matriz se exponencia porque se cuenta con que se le ha aplicado un log2
  
  y <- voom(GSEmat, design, plot = F)
  fit <- lmFit(y, design) # Ajuste a un modelo lineal
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  genes_limma.res <- eBayes(fit2, trend = trend)
  
  return(genes_limma.res)
}

# "Although a number of different multiple testing correction methods exists 
# (for instance see p.adjust documentation in R or permutation-based correction 
# methods), the most preferable approach is controlling FDR as it not only reduces
# false positives, but also minimises false negatives. "
# Jafari, M., & Ansari-Pour, N. (2019). Why, When and How to Adjust Your P Values?. 
# Cell journal, 20(4), 604–607. https://doi.org/10.22074/cellj.2019.5992

ED_Case_Ctrl <- dif_exp(d, "(AD - NC)", casos)
# Guardamos el objeto MArrayLM
saveRDS(ED_Case_Ctrl, file=paste(args[1],"_ED_Case_Ctrl.RDS",sep=""))
tt_c_ctrl <- topTable(ED_Case_Ctrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default


ED_sexdiff <- dif_exp(d, "(AD.female - NC.female) - (AD.male - NC.male)", group)
# Guardamos el objeto MArrayLM
saveRDS(ED_sexdiff, file=paste(args[1],"_ED_SexDiff.RDS",sep=""))
tt_sexdiff <- topTable(ED_sexdiff, number = Inf, adjust.method = "fdr") # Ajuste de BH por default

ED_FCase_FCtrl <- dif_exp(d, "(AD.female - NC.female)", group)
# Guardamos el objeto MArrayLM
saveRDS(ED_FCase_FCtrl, file=paste(args[1],"_ED_FCase_FCtrl.RDS",sep=""))
tt_FCase_FCtrl <- topTable(ED_FCase_FCtrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default

ED_MCase_MCtrl <- dif_exp(d, "(AD.male - NC.male)", group)
# Guardamos el objeto MArrayLM
saveRDS(ED_MCase_MCtrl, file=paste(args[1],"_ED_MCase_MCtrl.RDS",sep=""))
tt_MCase_MCtrl <- topTable(ED_MCase_MCtrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default



# Guardamos la tabla con los resultados del analisis de expresion diferencial
saveRDS(file="tt_case_ctrl.RDS", tt_c_ctrl)
saveRDS(file="tt_sex_case.RDS", tt_sexdiff)
saveRDS(file="tt_FCase_FCtrl.RDS", tt_FCase_FCtrl)
saveRDS(file="tt_MCase_MCtrl.RDS", tt_MCase_MCtrl)

rm(list=ls())
