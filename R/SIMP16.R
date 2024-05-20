library(ape) #install.packages("ape")
library(dplyr) #install.packages("dplyr")
library(ggplot2) #install.packages("ggplot2")
library(phangorn) #BiocManager::install("phangorn")
library(phyloseq) #BiocManager::install("phyloseq")
library(plotly) # install.packages("plotly")
library(tidyr)
library(vegan)
library(RCM) #BiocManager::install("RCM")
library(ANCOMBC) #BiocManager::install("ANCOMBC")
#install.packages("usethis")
#install.packages("devtools")
library("devtools")
library(ShortRead)
library(dada2)

#devtools::install_github("benjjneb/dada2", ref="v1.16")

# Save the current workspace to .RData
save.image()
file.edit("/media/aya/One Touch/Study 2/Data/.Rprofile")
savehistory("/media/aya/One Touch/Study 2/Data/file.Rhistory")




# 1st step is loading the DADA2 package and other packages to read the fastq files
path <- file.path("/media/aya/One Touch/Study 2/Data")  # CHANGE ME to the directory containing the fastq files after unzipping.
#list.files(path)


# 2nd step is to sort fastq files to make sure they are in right order to make sure that that forward and reverse trees are actually corresponding to each other.Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.names
# Define the path to the F folder
folder_path <- "/media/aya/One Touch/Study 2/Outputs"

# Combine the folder path with the file name
file_path <- file.path(folder_path, "Samples_names.csv")

# Write the data frame to the CSV file
write.csv(sample.names, file = file_path, row.names = FALSE)


#make list of number of sequences
readfiles_length <- countLines(fnFs) / 4
sum_readfiles_length <- sum(readfiles_length)


plot <- plotQualityProfile(fnFs)
data <- plot$data

#aggregate data for each sequencing cycle
df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
cycles <- sort(unique(data$Cycle))
for (cycle in cycles) {
    subdata <- data[data[, "Cycle"] == cycle, ]
    score <- list()
    #convert to list to calculate median
    for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata$Score[j], subdata$Count[j])))}
    temp = data.frame(Cycle=cycle, Count=sum(subdata$Count), Median=median(score), stringsAsFactors=FALSE)
    df <- rbind(df, temp)
}

#write output
write.table( t(df), file = paste0("qual_statsF",".tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)


#make list of number of sequences
readfiles_length <- countLines(fnRs) / 4
sum_readfiles_length <- sum(readfiles_length)


plot <- plotQualityProfile(fnRs)
data <- plot$data

#aggregate data for each sequencing cycle
df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
cycles <- sort(unique(data$Cycle))
for (cycle in cycles) {
    subdata <- data[data[, "Cycle"] == cycle, ]
    score <- list()
    #convert to list to calculate median
    for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata$Score[j], subdata$Count[j])))}
    temp = data.frame(Cycle=cycle, Count=sum(subdata$Count), Median=median(score), stringsAsFactors=FALSE)
    df <- rbind(df, temp)
}

#write output
write.table( t(df), file = paste0("qual_statsR",".tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

reticulate::repl_python()
reticulate::source_python("/home/aya/Documents/trim.py")

library(reticulate)

# Load pandas module
pd <- import("pandas")

# Load the data from the pickle file
data <- pd$read_pickle("/home/aya/trunc.pkl")

# Now you can use the data in R
print(data[1])


#4th step Filter and trim
#to keep sure that we keep only high quality data we will trim the final 10-15 base pairs of reads
# Place filtered files in filtered/ sub directory
#create new dir for the filtered files called "filtered" and rename the filtered file with f for forward and r for reverse
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names






















#Quality filter and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#head(out)



#5th step Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]


#6th step
#We now merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#7th step Construct an amplicon sequence variant table (ASV) table a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
file_path <- file.path("/media/aya/One Touch/Study 2/Outputs", "seqtabAllSamples_20.csv")
write.csv(seqtab, file = file_path, row.names = FALSE)
dim(seqtab)



# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#8th step Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
write.table(seqtab.nochim, "/media/aya/One Touch/Study 2/Outputs/asv_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)



#Track reeds through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "/media/aya/One Touch/Study 2/Outputs/Summary.csv", row.names = FALSE)

#Average Length
#seq_lengths <- nchar(getSequences(asv_table_filtered))
#average_length <- sum(seq_lengths) / length(seq_lengths)
#print(average_length)

# 9th step Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/media/aya/One Touch/Study 2/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa, file="/media/aya/One Touch/Study 2/Outputs/taxonomy.csv")
saveRDS(taxa, "ASVs_taxonomy138.rds")


# ---------------for shared---------------------
# Read the CSV file
asv <- read.csv("/media/aya/One Touch/Study 2/Outputs/asv.csv")
asv <- read.csv("/home/aya/seqtab.csv")

# Rename columns
new_column_names <- paste0("ASV", seq_along(asv))
colnames(asv) <- new_column_names

group_data <- read.csv("/home/aya/Downloads/Samples_names.csv", header = TRUE)
label_column <- rep(0.03, nrow(asv))
numAsvs_column <- rep(ncol(asv), nrow(asv))
asvshared <- cbind(label = label_column, Group = group_data[,1], numAsvs = numAsvs_column, asv)

# Write data frame to CSV with tab-separated values
write.table(asvshared, file = "/home/aya/ASVshared.shared", sep = "\t", row.names = FALSE, quote = FALSE)

#--------------------------------for taxonomy-----------------
# Read the tax file

new_data <- read.csv("/home/aya/Downloads/ASVs_taxonomy.csv")

new_data <- read.csv("/media/aya/One Touch/Study 2/Outputs/taxonomy.csv")

# Rename the first column to 'ASV'
names(new_data)[1] <- "ASV"

# Rename the rows of the first column
new_data$ASV <- paste0("ASV", 1:nrow(new_data))

# Sum each row and create a new column with the sums
data <- colSums(asvshared[4:ncol(asvshared)])

# Prepare a new data frame with only the sum column
size_data <- data.frame(size = data)

# If you want to combine the sum column with the original data and write it to a new file:
# Prepare a new data frame with the sum column added as the second column
contaxonomy <- cbind(ASV=new_data[, 1], size = size_data$size, new_data[, 2:ncol(new_data)])


# Merge columns from the second to the last, using ";" as the delimiter
merged_column <- apply(contaxonomy[, 3:ncol(contaxonomy)], 1, paste, collapse = ";")

# Create a new data frame with the merged column under the header "Taxonomy"
contax <- data.frame(Taxonomy = merged_column)

# Write the combined data to a new CSV file
write.table(cbind(contaxonomy[, 1:2], contax), file = "/media/aya/One Touch/Study 2/Outputs/contaxonomy.taxonomy", sep = "\t", row.names = FALSE, quote = FALSE)

taxoncom <- cbind(contaxonomy[, 1:2], contax)

#-------------------------------------------------------------------------
contaxfile="/home/aya/contaxonomy.taxonomy"
mapfile = "/home/aya/Metadata1.csv"
constaxfile = "/media/aya/One Touch/Study 2/Outputs/contaxonomy.taxonomy"
sharedfile = "/media/aya/One Touch/Study 2/Outputs/ASVshared.shared"
sharedfile = "/home/aya/ASVshared.shared"



#############Mothur to phyloseq####################

mothur_data <- import_mothur(
    mothur_shared_file = sharedfile,
    mothur_constaxonomy_file = constaxfile
)
map <- read.csv(mapfile,sep = "\t")
# transform dataframe into phyloseq class object.
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$sampleID

#checking the object "map"
str(map)
# Merge mothurdata object with sample metadata
mothur_merge <- merge_phyloseq(mothur_data, map)

phylo = mothur_merge
######################################## RCM ##################################################

# choose the order of the target columns
phylo@sam_data[["Host_disease"]] <- factor(phylo@sam_data[["Host_disease"]], levels = c("Null", "Low", "Medium", "High"))


asv = RCM(phylo)

plot(asv, samColour = "ants_binary")
plot(asv, samColour = "host_sex")
asv_confounder = RCM(phylo, confounders = "host_sex", prevCutOff = 0.5)
plot(asv_confounder, samColour = "Host_disease")
#plot(asv_confounder, samColour = "FinalComment")

######################################## Alpha Diversity ##################################################


# Check the available metadata columns in your phyloseq object
sample_data(phylo)

# Assuming 'Category' is the metadata column with four groups
# You may need to replace 'Category' with the actual column name in your data
metadata <- sample_data(phylo)
category <- metadata$Category

# Extract the OTU table from the phyloseq object
otu_table <- as(otu_table(phylo), "matrix")

# Perform Kruskal-Wallis test for each OTU
kruskal_results <- apply(otu_table, 1, function(otu_row) {
    kruskal.test(otu_row ~ category)
})

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(kruskal_results, function(x) x$p.value)

# Adjust p-values for multiple testing if necessary
# For example, using Bonferroni correction
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Identify significant OTUs
significant_otus <- which(adjusted_p_values < 0.05)

# Display significant OTUs and their adjusted p-values
significant_results <- data.frame(OTU = names(significant_otus), Adjusted_p_value = adjusted_p_values[significant_otus])
print(significant_results)




# Calculate Chao1 alpha diversity
alpha_data <- estimate_richness(phylo, measures = "Chao1")

# Convert the data to a format suitable for plotting
alpha_data_df <- as.data.frame(alpha_data)
alpha_data_df$Host_disease <- as.factor(sample_data(phylo)$Host_disease)

# Create boxplot
p <- ggplot(alpha_data_df, aes(x = Host_disease, y = Chao1, fill = Host_disease)) +
    geom_boxplot(outlier.shape = NA) +
    ylab("Chao1 Alpha Diversity") +
    xlab("Host_disease") +
    theme_minimal()

print(p)



# Calculate InvSimpson alpha diversity
alpha_data <- estimate_richness(phylo, measures = "InvSimpson")

# Convert the data to a format suitable for plotting
alpha_data_df <- as.data.frame(alpha_data)
alpha_data_df$Host_disease <- as.factor(sample_data(phylo)$Host_disease)

# Create boxplot
p <- ggplot(alpha_data_df, aes(x = Host_disease, y = InvSimpson, fill = Host_disease)) +
    geom_boxplot(outlier.shape = NA) +
    ylab("InvSimpson Alpha Diversity") +
    xlab("Host_disease") +
    theme_minimal()

print(p)


# Calculate Shannon alpha diversity
alpha_data <- estimate_richness(phylo, measures = "Shannon")

# Convert the data to a format suitable for plotting
alpha_data_df <- as.data.frame(alpha_data)
alpha_data_df$Host_disease <- as.factor(sample_data(phylo)$Host_disease)

# Create boxplot
p <- ggplot(alpha_data_df, aes(x = Host_disease, y = Shannon, fill = Host_disease)) +
    geom_boxplot(outlier.shape = NA) +
    ylab("Shannon Alpha Diversity") +
    xlab("Host_disease") +
    theme_minimal()

print(p)



#################################R 3.5 (Beta Diversity)

#permanova
install.packages("gplots")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(dplyr)
library(ape)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(pairwiseAdonis)
library(microbiome)


#beta diversity

#clr transformation
pseq.rel <- microbiome::transform(phylo, "clr")

otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

#clr distances
clr_dist_matrix <- phyloseq::distance(pseq.rel, method = "euclidean")
#clr permanova
h1 <- with(meta,how(nperm=10000,
                    within = Within(type = "free")))

adonis2(formula = clr_dist_matrix ~Category,data=meta, permutations = h1)
#adonis2(formula = clr_dist_matrix ~Comment,data=meta, permutations = h1)
adonis2(formula = clr_dist_matrix ~Category+FinalComment,data=meta, permutations = h1)

pairwise.adonis2(clr_dist_matrix ~ Category,data=meta)
pairwise.adonis2(clr_dist_matrix ~ FinalComment,data=meta)
pairwise.adonis2(clr_dist_matrix ~ Category+FinalComment,data=meta)

#########################################################################################################################################

#ANCOM BC
install.packages("gmp")
install.packages("nloptr")
install.packages("Rmpfr")
install.packages("lme4")
install.packages("lmerTest")
install.packages("Cairo")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")  # Replace with the version compatible with your R version

BiocManager::install("microbiome")

install.packages('DT')
library(DT)
install.packages('gt')
library(gt)
BiocManager::install("ANCOMBC", dependencies = TRUE)
library(ANCOMBC)

BiocManager::install("mia")
library(mia)
BiocManager::install("DECIPHER")





out = ancombc(phyloseq = phylo, formula = "host_sex + Host_disease",
              p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000,
              group = "Host_disease", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
              max_iter = 100, conserve = TRUE, alpha = 0.1, global = TRUE)
#ALL files in res
res = out$res
res_global = out$res_global

## ANCOMBC primary result
#Coefficients

tab_coef = res$lfc
col_name = colnames(tab_coef)
tab_coef %>% datatable(caption = "Coefficients from the Primary Result") %>%
    formatRound(col_name, digits = 2)
tab_html <-data.frame(tab_coef) %>%
    gt(rownames_to_stub = TRUE,caption = "Coefficients from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("Coeff.html"))

#SEs
tab_se = res$se
col_name = colnames(tab_se)
tab_se %>% datatable(caption = "SEs from the Primary Result") %>%
    formatRound(col_name, digits = 2)
tab_html <-data.frame(tab_se) %>%
    gt(rownames_to_stub = TRUE,caption = "SEs from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("SEs.html"))


#Test statistics
tab_w = res$W
col_name = colnames(tab_w)
tab_w %>% datatable(caption = "Test Statistics from the Primary Result") %>%
    formatRound(col_name, digits = 2)
tab_html <-data.frame(tab_w) %>%
    gt(rownames_to_stub = TRUE,caption = "Test Statistics from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("t.statistics.html"))

#P-values
tab_p = res$p_val
col_name = colnames(tab_p)
tab_p %>% datatable(caption = "P-values from the Primary Result") %>%
    formatRound(col_name, digits = 2)
tab_html <-data.frame(tab_p) %>%
    gt(rownames_to_stub = TRUE,caption = "P-values from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("p-values.html"))


#Adjusted p-values
tab_q = res$q_val
col_name = colnames(tab_q)
tab_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>%
    formatRound(col_name, digits = 2)
tab_html <-data.frame(tab_q) %>%
    gt(rownames_to_stub = TRUE,caption = "Adjusted p-values from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("q-values.html"))

#Differentially abundant taxa
tab_diff = res$diff_abn
col_name = colnames(tab_diff)
tab_diff %>%
    datatable(caption = "Differentially Abundant Taxa from the Primary Result")
tab_html <-data.frame(tab_diff) %>%
    gt(rownames_to_stub = TRUE,caption = "Differentially Abundant Taxa from the Primary Result") %>%
    as_raw_html()
write(tab_html,paste0("Differentially abundant taxa.html"))

#Bias-adjusted abundances
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0

# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(microbiome::abundances(phylo) + 1)
# Adjust the log observed abundances
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_obs_abn_adj[, 1:6], 2) %>%
    datatable(caption = "Bias-adjusted log observed abundances")
tab_html <-data.frame(round(log_obs_abn_adj, 2)) %>%
    gt(rownames_to_stub = TRUE,caption = "Bias-adjusted log observed abundances") %>%
    as_raw_html()
write(tab_html,paste0("Log observed abundances.html"))

