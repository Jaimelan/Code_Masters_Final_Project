#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: Volcano plots for DE results
##
## Purpose of script: Makes volcanoplots of MArrayLM DE objects
##
## Author: Jaime Llera
##
##
## Email: jlleraoy@gmail.com
##
## ---------------------------

library(argparse)
library(limma)
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-d", "--dir", type="character", default="Datasets", 
                    help="Set the directory containing the GSE studied",
                    metavar="dir name")
parser$add_argument("-e", "--exprmat", type="character", default=NULL, 
                    help="Differential expression matrix to make the plot of (MArrayLM object)",
                    metavar="expression matrix")



# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()



## Generacion de volcanoplots

fit1 <- readRDS(paste(args$d, "/", args$e ,sep=""))

dir.create("DE_Plots", recursive=T)
setwd("DE_Plots")

# Name of the contrast (contr)
contrname <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", args$e)

fit1$p.value <- p.adjust(fit1$p.value, method="BH")

png(filename=paste("volcano_",contrname,".png",sep=""))

volcanoplot(fit1, cex=.7) 

lp <- -log10(fit1$p.value)

sig_up <- lp > -log10(0.05) & fit1$coefficients > 0
sig_down <- lp > -log10(0.05) & fit1$coefficients < 0

points(fit1$coefficients[sig_up], lp[sig_up], pch = 16, cex = 0.7, col = "#55D500")
points(fit1$coefficients[sig_down], lp[sig_down], pch = 16, cex = 0.7, col = "#D500FF")
legend("bottomleft" , legend = c("Up", "Down"), col = c("#55D500", "#D500FF"), pch = 16)

dev.off()

