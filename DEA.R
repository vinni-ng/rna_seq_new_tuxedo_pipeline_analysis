#Install and load relevant R packages

install.packages("RSkittleBrewer")
install.packages("genefilter")
install.packages("dplyr")
install.packages("devtools")

library(RSkittleBrewer)
library(genefilter
library(dplyr)
library(devtools)

# Load phenotype data
pheno_data <- read.csv("data.csv")

# Read in expression data using ballgown
gtf_files <- ballgown(dataDir = "ballgown", samplePattern = "SRR", pData = pheno_data)

# Filter to remove non-abundance genes
gtf_files_filt <- subset(gtf_files, "rowVars(texpr(gtf_files)) > 1", genomesubset = TRUE)

# Perform differential expression analysis for transcripts and genes
results_transcripts <- stattest(gtf_files_filt, feature = "transcript", covariate = "Condition", getFC = TRUE, meas = "FPKM")
results_genes <- stattest(gtf_files_filt, feature = "gene", covariate = "Condition", getFC = TRUE, meas = "FPKM")

# Combine gene names and IDs with the results for transcripts
results_transcripts <- data.frame(geneNames = ballgown::geneNames(gtf_files_filt),
                                  geneIDs = ballgown::geneIDs(gtf_files_filt), results_transcripts)

# Sort the results based on p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <- arrange(results_genes, pval)

# Write the results to CSV files
write.csv(results_transcripts, "Overian_transcript_results.csv", row.names = FALSE)
write.csv(results_genes, "Overian_gene_results.csv", row.names = FALSE)

# Write differentially expressed transcripts to a separate CSV file
write.csv(subset(results_transcripts, results_transcripts$qval < 0.05), "dif_trans.csv", row.names = FALSE)

