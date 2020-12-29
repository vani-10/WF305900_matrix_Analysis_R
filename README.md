# WF305900_matrix_Analysis_R
Analysis of WF305900 using R
This is a R script to analyse Quant Matrix:
Load all the required packages:-
library(data.table)
library(ggplot2)
library(pheatmap)
library(limma)
raw_matrix <- read.delim('/Users/kalaivanikarthik/Downloads/WF305900_matrix.txt')
head(raw_matrix)#There are Peptide ,Protein Scores RTs and Intensities 
raw_matrix_intensities = raw_matrix[, grep("^Intensity", names(raw_matrix))]
head(raw_matrix_intensities)#There are only intensities now but should include peptide and protein column (analyte id columns)as well 
raw_matrix_intensities = cbind(raw_matrix[, 1:2], raw_matrix_intensities)
head(raw_matrix_intensities)
#Clean up the run ids
names(raw_matrix_intensities) = gsub("Intensity_|_scored", "", names(raw_matrix_intensities))
names(raw_matrix_intensities) = sapply(names(raw_matrix_intensities), function(x){unlist(strsplit(x, split = "\\."))[1]})
names(raw_matrix_intensities)
# check intensity sums and distributions - are there outliers to be treated with care/removed?
# par(mar = c(7,7,7,7))
barplot(colSums(raw_matrix_intensities[, 3:ncol(raw_matrix_intensities)]), main = "Raw intensity sums per run", las = 2)
# Check intensity distributions.
# That's best done on a clean matrix, with analyte annotation in the row.names
int_log10 = as.matrix(log10(raw_matrix_intensities[, 3:ncol(raw_matrix_intensities)]))
