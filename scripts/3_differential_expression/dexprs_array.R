#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: Differential expression analysis for microarrays
##
## Author: Jaime Llera
##
##
## Email: jlleraoy@gmail.com
##
## ---------------------------
##
##   
##
## ---------------------------
# Analisis de expresion diferencial en microarrays

library(limma)

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
datos <- readRDS(pDataname)

matname <- paste0(args[1], "_mat.RDS", sep="")
matriz <- readRDS(matname)




# Creamos la interaccion entre las variables de Ad y Sexo
datos$interaccion <- interaction(datos$casos, datos$sex)


# Funcion que realiza el ajuste con limma para realizar el contraste
dif_exp <- function(GSEmat,        #Matriz de datos
                    contraste,     #Contraste
                    grupos = NULL, #Variable donde están los elementos del contraste
                    lotes = NULL,  #Si lo hay, variable para el efecto batch
                    trend = FALSE) {
  
  grupos <- as.factor(grupos) #Los grupos tienen que ser un factor
  
  if (is.null(lotes)) {
    design <- model.matrix(~0+grupos)
    colnames(design) <- levels(grupos)
  } else {
    lotes <- as.factor(lotes)
    design <- model.matrix(~0+grupos+lotes)
    colnames(design) <- c(levels(grupos), levels(lotes)[-2])
    #colnames(design) <- c(levels(grupos), levels(lotes))
  }
  set.seed(123)
  cont.matrix <- makeContrasts(contrasts = contraste, levels = design)
  
  fit <- lmFit(GSEmat, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  genes_limma.res <- eBayes(fit2, trend = trend)
  
  return(genes_limma.res)
}


### Realizamos las comparaciones (AD vs NC y AD en mujeres vs AD en hombres):


ED_Case_Ctrl <- dif_exp(matriz, "(AD - NC)", datos$casos)
# Guardamos el objeto MArrayLM
saveRDS(ED_Case_Ctrl, file=paste(args[1],"_ED_Case_Ctrl.RDS",sep=""))
tt_c_ctrl <- topTable(ED_Case_Ctrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default


ED_2 <- dif_exp(matriz, "(AD.female - NC.female) - (AD.male - NC.male)", datos$interaccion)
# Guardamos el objeto MArrayLM
saveRDS(ED_2, file=paste(args[1],"_ED_SexDiff.RDS",sep=""))
tt_sexdiff <- topTable(ED_2, number = Inf, adjust.method="fdr" ) # Ajuste de fdr


ED_FCase_FCtrl <- dif_exp(matriz, "(AD.female - NC.female)", datos$interaccion)
# Guardamos el objeto MArrayLM
saveRDS(ED_FCase_FCtrl, file=paste(args[1],"_ED_FCase_FCtrl.RDS",sep=""))
tt_FCase_FCtrl <- topTable(ED_FCase_FCtrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default

ED_MCase_MCtrl <- dif_exp(matriz, "(AD.male - NC.male)", datos$interaccion)
# Guardamos el objeto MArrayLM
saveRDS(ED_MCase_MCtrl, file=paste(args[1],"_ED_MCase_MCtrl.RDS",sep=""))
tt_MCase_MCtrl <- topTable(ED_MCase_MCtrl, number = Inf, adjust.method = "fdr") # Ajuste de BH por default


# Guardamos la tabla con los resultados del analisis de expresion diferencial
saveRDS(file="tt_case_ctrl.RDS", tt_c_ctrl)
saveRDS(file="tt_sex_case.RDS", tt_sexdiff)
saveRDS(file="tt_FCase_FCtrl.RDS", tt_FCase_FCtrl)
saveRDS(file="tt_MCase_MCtrl.RDS", tt_MCase_MCtrl)


rm(list=ls())
