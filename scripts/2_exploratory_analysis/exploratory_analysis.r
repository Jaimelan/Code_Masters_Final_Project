#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: Exploratory analysis
##
## Purpose of script: Plots information and basic statistics of a given data matrix (passed as an argument)
##
## Author: José F. Catalá, modified by Jaime Llera
##
## Date Created: 05-04-2022
##
## Email: jlleraoy@gmail.com
##
## ---------------------------
##
## Notes: The argument [1] must contain a data matrix (relative route from the folder of execution of this script)
##
##
## ---------------------------


require(tidyverse)
require(data.table)
require(knitr)
require(ggplot2)
require(ggdendro)
require(dplyr)
require(ggpubr)

## ---------------------------

args = commandArgs(trailingOnly = T)

# Evaluamos que en los argumentos se aporte una carpeta y una matriz de datos
if (length(args) != 1) {
  stop("Only 1 must be supplied: 1. experiment_id", call. = FALSE)
} else if (length(args) == 1) {
  print("Procesando datos...")
}

# elegimos la carpeta de trabajo, cargamos los datos
setwd(args[1])

pDataname <- paste0(args[1], "_pData.RDS", sep = "")
pDa <- readRDS(pDataname)

matname <- paste0(args[1], "_mat.RDS", sep = "")
datamatrix <- readRDS(matname)


## Funciones para los plots ##

# Plots de clustering por distancia de correlación y distancia euclídea
clustering_ <- function(GSE_object, GSEdata_object, tmc = 2) {
  # Clustering with correlation distance
  factor_casos <- GSEdata_object$casos
  hc <- hclust_corr_dist(GSE_object)
  print(
    hclust_plot(
      hc,
      factor_casos,
      samples_order = GSEdata_object$geo_accession,
      size = tmc,
      main = "Clustering con distancia de correlación"
    )
  ) +
    theme(
      plot.title = element_text(size=25),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  ggsave(filename = "clust_corr.png",
         width = 16,
         height = 10) # Se añade guardado de los plots
  he <- hclust_euc_dist(GSE_object)
  print(
    hclust_plot(
      he,
      factor_casos,
      samples_order = GSEdata_object$geo_accession,
      size = tmc,
      main = "Clustering con distancia euclídea"
    )
  ) +
    theme(
      plot.title = element_text(size=25),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  ggsave(filename = "clust_eucl.png",
         width = 16,
         height = 10) # Se añade guardado de los plots
}

# Función que genera los gráficos de clustering
hclust_plot <-
  function (hclt,
            group = NULL,
            samples_order = NULL,
            rotate = TRUE,
            segments = TRUE,
            size = 5,
            main = "Clustering",
            xlab = "Muestras",
            ylab = "",
            main_legend = "Condición")
  {
    dend <- as.dendrogram(hclt)
    dend_data <- dendro_data(dend, type = "rectangle")
    plot <- ggplot(dend_data$segments) + labs(title = main, x = xlab,
                                              y = ylab) + scale_colour_discrete(name = main_legend, type =
                                                                                  c("#440154", "#35B779")) +
      geom_segment(aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ), size=1.5) + theme(legend.text = element_text(size=17),
                 legend.title = element_text(size=17),
                 axis.title=element_text(size=17))
    if (is.null(group) == FALSE) {
      if (is.null(samples_order) == TRUE) {
        stop("samples_order is needed")
      }
      if (length(samples_order) != length(group) | length(group) !=
          length(dend_data$labels$label)) {
        stop("lengths are not equal")
      }
      group <- group[base::match(as.vector(dend_data$labels$label),
                                 samples_order)]
      if (rotate) {
        #IN USE
        plot <- plot + scale_y_continuous(expand = c(0.2, 0)) + 
          geom_text(data = dend_data$labels, aes(x,y, label = label, color = group),
                                                       hjust = 1,
                                                       size = 5,
                                                       angle = 90
                                                     )
      }
      else {
        plot <- plot + coord_flip() + scale_y_reverse(expand = c(0.15,
                                                                 0)) + geom_text(
                                                                   data = dend_data$labels,
                                                                   aes(x,
                                                                       y, label = label, color = group),
                                                                   hjust = 0,
                                                                   size = size
                                                                 )
      }
      if (segments) {
        #IN USE
        plot <- plot + geom_segment(
          data = dend_data$segments %>%
            dplyr::filter(yend == 0) %>% dplyr::left_join(dend_data$labels,
                                                          by = "x"),
          aes(x = x, y = y.x, xend = xend, yend = yend, color = group),
          size = 1.5
        )
      }
    }
    else {
      if (rotate) {
        plot <- plot + scale_y_continuous(expand = c(0.2,
                                                     0)) + geom_text(
                                                       data = dend_data$labels,
                                                       aes(x,
                                                           y, label = label),
                                                       hjust = 1,
                                                       size = size,
                                                       angle = 90
                                                     )
      }
      else {
        plot <- plot + coord_flip() + scale_y_reverse(expand = c(0.15,
                                                                 0)) + geom_text(
                                                                   data = dend_data$labels,
                                                                   aes(x,
                                                                       y, label = label),
                                                                   hjust = 0,
                                                                   size = size
                                                                 )
      }
    }
    plot
  }

# Clustering por distancia de correlación
hclust_corr_dist <- function (GSE)
{
  if ((class(GSE) %in% c("ExpressionSet", "matrix")) == FALSE) {
    stop("GSE is not an ExpressionSet or a matrix")
  }
  if (class(GSE) == "ExpressionSet") {
    correlacion <- cor(Biobase::exprs(GSE))
  }
  else {
    correlacion <- cor(GSE)
  }
  distancia <- as.dist((1 - correlacion) / 2)
  hc <- hclust(distancia)
  if (all(hc$labels == colnames(GSE)) == FALSE) {
    stop("hc$labels != colnames(GSE)")
  }
  return(hc)
}

# Clustering por distancia euclídea
hclust_euc_dist <- function (GSE)
{
  if ((class(GSE) %in% c("ExpressionSet", "matrix")) == FALSE) {
    stop("GSE is not an ExpressionSet or a matrix")
  }
  if (class(GSE) == "ExpressionSet") {
    distancia <- dist(t(Biobase::exprs(GSE)), method = "euclidean")
  }
  else {
    distancia <- dist(t(GSE), method = "euclidean")
  }
  he <- hclust(distancia)
  if (all(he$labels == colnames(GSE)) == FALSE) {
    stop("hc$labels != colnames(GSE)")
  }
  return(he)
}

# Función que genera el gráfico del PCA. Necesita un factor (g1) para colorear 
# las muestras. Admite un segundo factor que representará en forma de "shape" de
# los puntos.
plot_pca_genes <-
  function (mi.pca,
            g1,
            g2 = NULL,
            main = NULL,
            g1_name = NULL,
            g2_name = NULL,
            size = 3,
            addNames = FALSE)
  {
    if (addNames) {
      PC1_PC2 <- ggplot() + geom_text(aes(
        x = mi.pca$scores[,
                          1],
        y = mi.pca$scores[, 2],
        colour = g1,
        label = colnames(mi.pca$Xoff)
      ),
      size = size) + xlab(paste(
        "PC1: ",
        round(mi.pca$var.exp[1, 1] * 100),
        "% variación explicada",
        sep = ""
      )) + ylab(paste(
        "PC2: ",
        round(mi.pca$var.exp[2, 1] * 100),
        "% explained variance",
        sep = ""
      )) + scale_colour_discrete(name = g1_name, type = c("#440154", "#1F9E89"))
      PC3_PC2 <- ggplot() + geom_text(aes(
        x = mi.pca$scores[,
                          3],
        y = mi.pca$scores[, 2],
        colour = g1,
        label = colnames(mi.pca$Xoff)
      ),
      size = size) + xlab(paste(
        "PC3: ",
        round(mi.pca$var.exp[3,
                             1] * 100),
        "% variación explicada",
        sep = ""
      )) + ylab(paste(
        "PC2: ",
        round(mi.pca$var.exp[2, 1] * 100),
        "% explained variance",
        sep = ""
      )) + scale_colour_discrete(name = g1_name, type = c("#440154", "#1F9E89"))
    }
    else {
      if (is.null(g2) == TRUE) {
        PC1_PC2 <- ggplot() + geom_point(aes(
          x = mi.pca$scores[,
                            1],
          y = mi.pca$scores[, 2],
          colour = g1
        ),
        size = size) +
          xlab(paste(
            "PC1: ",
            round(mi.pca$var.exp[1, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          ylab(paste(
            "PC2: ",
            round(mi.pca$var.exp[2, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          scale_colour_discrete(name = g1_name)
        PC3_PC2 <- ggplot() + geom_point(aes(
          x = mi.pca$scores[,
                            3],
          y = mi.pca$scores[, 2],
          colour = g1
        ),
        size = size) +
          xlab(paste(
            "PC3: ",
            round(mi.pca$var.exp[3, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          ylab(paste(
            "PC2: ",
            round(mi.pca$var.exp[2, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          scale_colour_discrete(name = g1_name)
      }
      else {
        # IN USE
        PC1_PC2 <- ggplot() + geom_point(aes(
          x = mi.pca$scores[, 1],
          y = mi.pca$scores[, 2],
          colour = g1,
          shape = g2
        ),
        size = size) + xlab(paste(
          "PC1: ",
          round(mi.pca$var.exp[1,
                               1] * 100),
          "% variación explicada",
          sep = ""
        )) +
          ylab(paste(
            "PC2: ",
            round(mi.pca$var.exp[2, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          scale_colour_discrete(name = g1_name, type = c("#440154", "#1F9E89")) + scale_shape_discrete(name = g2_name)
        PC3_PC2 <- ggplot() + geom_point(aes(
          x = mi.pca$scores[,
                            3],
          y = mi.pca$scores[, 2],
          colour = g1,
          shape = g2
        ),
        size = size) + xlab(paste(
          "PC3: ",
          round(mi.pca$var.exp[3,
                               1] * 100),
          "% variación explicada",
          sep = ""
        )) +
          ylab(paste(
            "PC2: ",
            round(mi.pca$var.exp[2, 1] *
                    100),
            "% variación explicada",
            sep = ""
          )) +
          scale_colour_discrete(name = g1_name, type = c("#440154", "#1F9E89")) + scale_shape_discrete(name = g2_name)
      }
    }
    figure <- ggarrange(
      PC1_PC2,
      PC3_PC2,
      ncol = 2,
      nrow = 1,
      common.legend = TRUE,
      legend = "bottom"
    )
    annotate_figure(figure, top = text_grob(main))
  }

# Calcula el PCA de una matriz o ExpressionSet
pcaGenes <- function (ob)
{
  if ((class(ob) %in% c("ExpressionSet", "matrix")) == FALSE) {
    stop("X is not an ExpressionSet or a matrix")
  }
  if (class(ob) == "ExpressionSet") {
    X <- (t(Biobase::exprs(ob)))
  }
  else {
    X <- t(ob)
  }
  X <- as.matrix(X)
  n <- ncol(X)
  p <- nrow(X)
  offset <- apply(X, 2, mean)
  Xoff <- X - (cbind(matrix(1, p, 1)) %*% rbind(offset))
  eigen <- eigen(Xoff %*% t(Xoff) / (p - 1))
  var <-
    cbind(eigen$values / sum(eigen$values), cumsum(eigen$values / sum(eigen$values)))
  loadings2 <- eigen$vectors
  scores2 <- t(Xoff) %*% loadings2
  normas2 <- sqrt(apply(scores2 ^ 2, 2, sum))
  scores1 <- loadings2 %*% diag(normas2)
  loadings1 <- scores2 %*% diag(1 / normas2)
  Xoff <- t(Xoff)
  output <- list(eigen, var, scores1, loadings1, Xoff)
  names(output) <- c("eigen", "var.exp", "scores", "loadings",
                     "Xoff")
  if (all(colnames(output$Xoff) == colnames(ob)) == FALSE) {
    stop("No coinciden los nombres de las columnas")
  }
  return(output)
}



## Case-sex info & plot
table(pDa$casos, pDa$sex)#, col.names = c("Sex","N"))
ggplot(data = pDa, aes(x = casos, fill = sex)) +
  geom_bar() +
  xlab("Condición") +
  ylab("N") +
  scale_fill_discrete(name = "Sexo", type = c("#440154", "#1F9E89"))

ggsave(filename = "sex_info.png",
       width = 4,
       height = 5)

#Boxplot

png(filename = "boxplot.png")
par(mar = c(10, 4.1, 4.1, 2.1))
par(cex.axis = 1.3)
boxplot(
  datamatrix,
  main = paste("Boxplot de", args[1]),
  las = 2,
  col = "#fde725"
)
dev.off()


# Clustering
clustering_(datamatrix, pDa)

# PCA
print("Realizar PCA...")

plot_pca_genes(
  pcaGenes(datamatrix),
  pDa$casos,
  pDa$sex,
  g1_name = "Cases",
  g2_name = "Sex",
  main = "Análisis de componentes principales (PCA)"
)
ggsave(filename = "PCA.png",
       width = 10,
       height = 5)
