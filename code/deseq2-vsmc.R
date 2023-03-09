library("DESeq2")


#----------------------------------------------------------------------------------#
################################ miRCTRL VS miR323: ################################
#----------------------------------------------------------------------------------#

gencode_v26_gft_table <- read.delim(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/gencode_v26_gtf_table.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
annotation <- gencode_v26_gft_table[, c(4, 7, 6, 1:3, 5)]
colnames(annotation) <- c("ENSEMBL", "name", "type", "chr", "start", "end", "str")

#----------------------------------------------------------------------------------#
#--------------------------- Metadata & running DESeq2: ---------------------------#
#----------------------------------------------------------------------------------#

type <- c(rep("paired-end", 6))
patient <- rep(c("p1", "p2", "p3"), 2)
condition <- c(rep("miRCTRL", 3), rep("miR323", 3))
# Count data + meta data:
countTable <- read.csv(
  "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/HSVSMC_miRNAOE_counts.csv",
  header = TRUE
)
countTableDESeq <- sapply(countTable[, c(17, 28, 39, 11, 22, 33)], as.integer)
row.names(countTableDESeq) <- countTable[, 1]
metadata <- data.frame(
  row.names = colnames(countTableDESeq),
  condition = condition,
  type = type,
  patient = patient
)

# Making DESeqDataSet object - adding design (correting for "patient"):
dds <- DESeqDataSetFromMatrix(
  countData = countTableDESeq,
  colData = metadata,
  design = ~ patient + condition
)

# Setting FBS02 as ref (can't use "contrast" due to running apeglm shrinkage):
dds$condition <- relevel(dds$condition, ref = "miRCTRL")
# Running DESeq2 with a Wald test:
dds <- DESeq(dds, test = "Wald")



#----------------------------------------------------------------------------------#
#-------------------------- Making DESEqResults objects: --------------------------#
#----------------------------------------------------------------------------------#

### WITHOUT LFCTHRESHOLD: ###
# With indep.filt, BH, alpha matched to FDR cuttoff
res <- results(
  dds,
  independentFiltering = TRUE,
  pAdjustMethod = "BH", # default
  alpha = 0.01 # To match the intended FDR cutoff
)

# Performing LFC shrinkage - WITHOUT LFC:
ddsShrink <- lfcShrink(
  dds = dds,
  coef = "condition_miR323_vs_miRCTRL",
  res = res,
  type = "apeglm",
  svalue = FALSE,
  returnList = FALSE
)

# Filtering for significant FDR:
ddsShrinkTibble <- ddsShrink %>%
  data.frame() %>%
  rownames_to_column(var = "EnsID") %>%
  as_tibble()

sigShrink <- ddsShrinkTibble %>%
  filter(padj < 0.01) %>%
  merge(annotation, by = 1, all.x = FALSE)

write.csv(
  sigShrink,
  file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/DESeq2_sig_noLFCthreshold/miR323vsmiRCTRL_noLFCthreshold.csv"
)



### WITH LFCTHRESHOLD: ###
# With indep.filt, BH, alpha matched to FDR cuttoff:
resLFCthresh <- results(
  dds,
  independentFiltering = TRUE,
  pAdjustMethod = "BH", # default
  alpha = 0.01, # To match the intended FDR cutoff
  lfcThreshold = 1,
  altHypothesis = "greaterAbs"
)
# Performing LFC shrinkage - WITH LFC:
ddsShrinkLFCthresh <- lfcShrink(
  dds = dds,
  coef = "condition_miR323_vs_miRCTRL",
  res = resLFCthresh,
  type = "apeglm",
  svalue = FALSE,
  returnList = FALSE
)

# Filtering for significant FDR:
ddsShrinkTibbleLFCthresh <- ddsShrinkLFCthresh %>%
  data.frame() %>%
  rownames_to_column(var = "EnsID") %>%
  as_tibble()

sigShrinkLFCthresh <- ddsShrinkTibbleLFCthresh %>%
  filter(padj < 0.01) %>%
  merge(annotation, by = 1, all.x = FALSE)

write.csv(
  sigShrinkLFCthresh,
  file = "/Users/effieklimi/Documents/PhD/miRNA screening paper/HSVSMC RNA sequencing/DESeq2_sig_withLFCthreshold/miR323vsmiRCTRL_withLFCthreshold.csv"
)


####################################################################################
###################### NOTES, EXPLANATIONS & DATA EXPLORATION ######################
####################################################################################
# Above: if independentFiltering	= TRUE:                                        #
# The threshold chosen (the vertical line) is the lowest quantile (theta) for      #
# which the number of rejections is within 1 residual standard deviation to the    #
# peak of the curve.                                                               #
#                                                                                  #
# DESeq2 helps reduce the number of genes tested by removing those genes unlikely  #
# to be significantly DE prior to testing, such as those with low number of counts #
# and outlier samples (gene-level QC).                                             #
#                                                                                  #
# Genes with a mean expression value under a certain threshold are removed.        #
# Such filtering is permissible only if the filter criterion is independent of     #
# the actual test statistic, otherwise, the filtering would invalidate the test    #
# and consequently the assumptions of the FDR procedure. This is why filtering     #
# is done on the average expression over all samples, irrespective of biological   #
# condition: this filter is blind to the assignment of samples to the treatment    #
# and control group and hence independent.                                         #
#                                                                                  #
# Read comment by Bioconductor core member James W. MacDonald on this topic here:  #
# https://support.bioconductor.org/p/101504/#101507                                #
#                                                                                  #
# alpha was set to match the intended (future) FDR cutoff of 0.01 based on advice  #
# coming from ?results (details on the alpha argument):                            #
# "alpha: the significance cutoff used for optimizing the independent filtering    #
# (by default 0.1). If the adjusted p-value cutoff (FDR) will be a value other     #
# than 0.1, alpha should be set to that value."                                    #
####################################################################################
####################################################################################
# NOTE: Shrinking the log2 fold changes will not change the total number of        #
# genes that are identified as significantly differentially expressed. The         #
# shrinkage of fold change is to help with downstream assessment of results.       #
# For example, if you wanted to subset your significant genes based on fold        #
# change for further evaluation, you may want to use shruken values.               #
# Additionally, for functional analysis tools such as GSEA which require fold      #
# change values as input you would want to provide shrunken values.                #
####################################################################################
####################################################################################
# What does FDR < 0.01 mean? By setting the FDR cutoff to < 0.01, we’re saying     #
# that the proportion of false positives we expect amongst our differentially      #
# expressed genes is 1%. For example, if you call 500 genes as differentially      #
# expressed with an FDR cutoff of 0.01, you expect 5 of them to be false positives.#
####################################################################################
#                                                                                  #
# % of genes filtered out + mean threshold for filtering out lowly expr genes:     #
metadata(res)$filterThreshold # <~5.46 basemean reject = 72.8% of genes            #
#                                                                                  #
# Visualise quantiles (theta) and number of rejected genes:                        #
as_tibble(metadata(res)$filterNumRej) %>% #
  ggplot(aes(x = theta, y = numRej)) + #
  geom_point() + #
  geom_vline( #
    xintercept = metadata(res)$filterTheta, #
    color = "red" #
  ) #
#                                                                                  #
####################################################################################

