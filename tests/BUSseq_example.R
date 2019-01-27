#######################################
# Apply BUSseq to the simulation data #
#######################################
library(BUSseq)
ObservedCounts <- BUSseqfits_example$CountData_raw
BUSseqfits_res <- BUSseq_MCMC(ObservedCounts, n.celltypes = 4, n.iterations = 500, 
                          seed = 123)

################################################
# Extract Estimates from the BUSseqfits Object #
################################################

#return cell type indicators
w.est <- celltypes(BUSseqfits_res)
table(w.est)

#return the intercept and odds ratio the logistic regression
#for the dropout events
gamma.est <- dropout_coefficient_values(BUSseqfits_res)

#return the baseline expression values
alpha.est <-  baseline_expression_values(BUSseqfits_res)

#return the cell-type effects
beta.est <- celltype_effects(BUSseqfits_res)

#return the cell-specific global effects
delta.est <- cell_effect_values(BUSseqfits_res)

#return the location batch effects
nu.est <- location_batch_effects(BUSseqfits_res)

#return the scale batch effects
phi.est <- overdispersions(BUSseqfits_res)

#return the intrinsic genes
D.est <- intrinsic_genes_BUSseq(BUSseqfits_res)

#return the BIC score
BIC <- BIC_BUSseq(BUSseqfits_res)

#return the raw read count matrix
CountData_raw <- raw_read_counts(BUSseqfits_res)

#return the corrected read count matrix
CountData_imputed <- imputed_read_counts(BUSseqfits_res)

#return the corrected read count matrix
CountData_corrected <- corrected_read_counts(BUSseqfits_res)

