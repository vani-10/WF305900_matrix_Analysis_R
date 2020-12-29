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
# you wanna remove any NAs or Inf values (replace with 0s)
int_log10_clean = copy(int_log10)
int_log10_clean[is.infinite(int_log10_clean)] <- 0 
int_log10_clean[is.na(int_log10_clean)] <- 0 
# Now this should be clean for plotting in different ways
boxplot(int_log10_clean, main = "Raw intensity (log10) distributions per run", las = 2)
# at least two runs look funky
# you could filter them out e.g. my a criterion on the minimum median log10 intensity:
boxplot(int_log10_clean, main = "Raw intensity (log10) distributions per run", las = 2)
cutoff = 5.6
abline(cutoff,0, col = "red", lty = 2)
goodruns = apply(int_log10_clean, 2, median) > cutoff
int_log10_clean_filt = int_log10_clean[, goodruns]
# note which ones were removed
message("removed ", paste(colnames(int_log10_clean)[!goodruns]), " from the analysis due to low log intensity mean below ", cutoff)
# re-plot after removing funky runs
boxplot(int_log10_clean_filt, main = "Raw intensity (log10) distributions per run", las = 2)
# Clustered heatmap
pheatmap(int_log10_clean_filt, show_rownames = F, main = "Precursor level heatmap, log10(Precursor.Intensity)")
dev.copy(pdf, "Heatmap_precursorintensities_log10.pdf", height = 12, width = 18)
dev.off()
## How does this look with cross-run intensity normalization?
int_log10_clean_filt_norm = limma::normalizeQuantiles(int_log10_clean_filt)
par(mfrow = c(2,1))
boxplot(int_log10_clean_filt, main = "Raw intensity (log10) per run", las = 2)
boxplot(int_log10_clean_filt_norm, main = "Quantile-normalized intensity (log10) per run", las = 2)
pheatmap(int_log10_clean_filt_norm, show_rownames = F, main = "Precursor level heatmap, log10(Precursor.Intensity)")
dev.copy(pdf, "Heatmap_precursorintensities_log10_normalized.pdf", height = 12, width = 18)
dev.off()
## Write out a table with normalized quantities
write.table(int_log10_clean_filt_norm, "WF305900_matrix_clean_normalized_log10.tsv", sep = "\t", row.names = T, quote = F)