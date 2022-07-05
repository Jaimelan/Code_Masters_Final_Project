## ---------------------------
##
## Script name: microRNA annotation
##
## Purpose of script: fundamentally, receives a list of microRNAs 
##                    and gets the Gene they target
## Author: Jaime Llera
##
## Date Created: 27-05-2022
##
## ---------------------------
# get_multimir from the multiMiR library allows for an enrichment using several 
# microRNA-target databases and returns a table with target genes for the 
# microRNAs given as well as information related with the annotation.

# change the directory to the output folder
setwd("~/Results/MA_SexDiff")

library(multiMiR)

# first we read the file containing all top and down microRNAs as a vector
topmirs <- scan(file="topmirs_sexdiff.txt", what = character())
botmirs <- scan(file="botmirs_sexdiff.txt", what = character())

# we make a query with get_multimir in the databases used for the function (
# mirecords, mirtarbase and tarbase. The return is a mmquery_bioc object.
top_enrichment <- get_multimir(mirna=topgenes)
bot_enrichment <- get_multimir(mirna=botgenes)

# we get a list of unique gene symbols associated with our microRNAs
top_unique_genes <- unique(top_enrichment@data$target_symbol)
bot_unique_genes <- unique(bot_enrichment@data$target_symbol)


write(top_unique_genes, file="top_unique_genes.txt")
write(bot_unique_genes, file="bot_unique_genes.txt")
