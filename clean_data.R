library("dndscv")
# Steps to find inference for dnds model provided in the paper 
data("dataset_simbreast", package="dndscv")
dndsout = dndscv(mutations)
my_data <- read.delim("/Users/atishnasamantaray/Downloads/Team_2_HNSC/TCGA.HNSC.mutations.txt")
#cleaning data to bring it to desired format
clean_data <- my_data[,c("patient_id", "Chromosome","Start_Position", "Reference_Allele","Tumor_Seq_Allele2","Variant_Type")]
clean_data <- clean_data[clean_data$Variant_Type == "SNP", ]
drop <- c('Variant_Type')
clean_data <- clean_data[,!(names(clean_data) %in% drop)]
colnames(clean_data) <- c("sampleID", "chr", "pos","ref","mut")
clean_data <- clean_data[!clean_data$ref == "-", ]
clean_data <- clean_data[!clean_data$mut == "-", ]
#calculating dnds and comparing it with that of paper to check if we are doing it correctly
dndsout = dndscv(clean_data)
dndsout = dndscv(mutations)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
#Checking the most significant genes
signif_genes = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes)
print(dndsout$globaldnds)
head(dndsout$annotmuts)
#Checking if dN/dS model is actually suitable for our mutations data
print(dndsout$nbreg$theta)
#Finding the most significantly mutated genes by dnds
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)

# Work for our own model
clean_data <- my_data[,c("patient_id", "Chromosome","Start_Position", "Reference_Allele","Tumor_Seq_Allele2","Variant_Type","Variant_Classification","Hugo_Symbol")]
clean_data <- clean_data[clean_data$Variant_Type == "SNP", ]

library('dplyr')
clean_summary = clean_data %>%
  group_by(Hugo_Symbol,Variant_Classification)  %>%
  summarize(count = n())

clean_summary <- clean_summary[!clean_summary$Variant_Classification == "3'UTR", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "5'UTR", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "Intron", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "RNA", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "5'Flank", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "3'Flank", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "Translation_Start_site", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "Translation_Start_Site", ]
clean_summary <- clean_summary[!clean_summary$Variant_Classification == "Nonstop_Mutation", ]
unique(clean_summary[,2])
head(mutations)

library("tidyr")
clean_summary <- clean_summary %>%
  pivot_wider(names_from = Variant_Classification, values_from = count, values_fill = 0)
#find the most 500 mutated data: BASELINE
clean_summary$most_mutated <- clean_summary$Missense_Mutation + clean_summary$Nonsense_Mutation + clean_summary$Silent + clean_summary$Splice_Site
clean_summary <-clean_summary %>% arrange(desc(most_mutated)) %>% head(500) 
clean_summary$Nonsyno = clean_summary$Nonsense_Mutation + clean_summary$Missense_Mutation + clean_summary$Splice_Site
clean_summary$dnds = clean_summary$Nonsyno / clean_summary$Silent
clean_summary = clean_summary[!clean_summary$dnds == "Inf",]
clean_summary$missense_dnds =clean_summary$Missense_Mutation / clean_summary$Silent
clean_summary$nonsense_dnds =clean_summary$Nonsense_Mutation / clean_summary$Silent
clean_summary$Splice_Site_dnds =clean_summary$Splice_Site / clean_summary$Silent

#clean data using chi-square test  
p_values <- apply(clean_summary[,8:11], 1, function(x) {
  chisq.test(x)$p.value})

# Perform FDR correction on the p-values
q_values <- p.adjust(p_values, method="fdr")

# Add the q-values to the dataframe
clean_summary$q_values <- q_values

# Filter the dataframe to keep only the significant results (q-value < 0.05, for example)
sig_results <- clean_summary[clean_summary$q_values < 0.01, ]
sig_results 


# possion on the columns bc of the rare events resample rare events
# resample: improvement
set.seed(123)
new_missense <- sapply(clean_summary$Missense_Mutation, function(x) rpois(10, lambda = x))
new_nonsense <- sapply(clean_summary$Nonsense_Mutation, function(x) rpois(10, lambda = x))
new_splice <- sapply(clean_summary$Splice_Site, function(x) rpois(10, lambda = x))
new_silent <- sapply(clean_summary$Silent, function(x) rpois(10, lambda = x))

resample <- data.frame(Gene =  rep(clean_summary$Hugo_Symbol, each = 10), 
                           Missense = c(new_missense), 
                           Nonsense = c(new_nonsense),
                           Splice = c(new_splice),
                           Silent = c(new_silent))
clean_resampled <- resample %>%
  group_by(Gene) %>%
  summarize(total_mis = sum(Missense),
            total_non = sum(Nonsense),
            total_spl = sum(Splice),
            total_silent = sum(Silent))

#get dnds from resampled
clean_resampled$Nonsyno = clean_resampled$total_mis +  clean_resampled$total_non +  clean_resampled$total_spl
clean_resampled$dnds = clean_resampled$Nonsyno /  clean_resampled$total_silent
clean_resampled$missense_dnds =clean_resampled$total_mis / clean_resampled$total_silent
clean_resampled$nonsense_dnds =clean_resampled$total_non / clean_resampled$total_silent
clean_resampled$Splice_Site_dnds =clean_resampled$total_spl / clean_resampled$total_silent
clean_resampled = clean_resampled[!clean_resampled$dnds == "Inf",]

# get significant genes using resample data
p_values <- apply(clean_resampled[,7:10], 1, function(x) {
  chisq.test(x)$p.value})

# Perform FDR correction on the p-values
q_values <- p.adjust(p_values, method="fdr")

# Add the q-values to the dataframe
clean_resampled$q_values <- q_values

# Filter the dataframe to keep only the significant results (q-value < 0.05, for example)
sig_results_resample <- clean_resampled[clean_resampled$q_values < 0.01, ]
sig_results_resample


# count of each selection
count_neg = sum(clean_summary$dnds < 1)
count_neg
count_neutral = sum(clean_summary$dnds == 1)
count_neutral
count_pos = sum(clean_summary$dnds > 1)
count_pos


# t test z test
library(BSDA)
dnds_res = clean_summary$Nonsyno / clean_summary$Silent
dnds_res <- dnds_res[!is.na(dnds_res) & !is.infinite(dnds_res)]
dnds_ref = (sel_cv$n_mis + sel_cv$n_non + sel_cv$n_spl) / sel_cv$n_syn 
dnds_ref <- dnds_ref[!is.na(dnds_ref) & !is.infinite(dnds_ref)]
t.test(dnds_res, dnds_ref)
z.test(dnds_res, dnds_ref, sigma.x=sd(dnds_res), sigma.y=sd(dnds_ref))
                     
                     
                     
#hypermutator
install.packages("org.Hs.eg.db")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

clean_sum <- clean_summary[!clean_summary$dnds == "Inf", ]
clean_summary$most_mutated <- clean_summary$Missense_Mutation + clean_summary$Nonsense_Mutation +clean_summary$Silent  + clean_summary$Splice_Site
clean_summary <- clean_summary %>% arrange(desc(most_mutated)) %>% head(500)
# Load required packages
library(BiocManager)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(biomaRt)
# define the Ensembl dataset and the mart
ensembl_dataset <- "hsapiens_gene_ensembl"
ensembl_mart <- useMart(biomart = "ensembl", dataset = ensembl_dataset)
# define the list of genes
gene_list <- clean_resampled[[1]]
gene_list <- head(gene_list, 500)
# retrieve the coding sequence length for each gene
coding_seq_length_list <- list()
for (gene in gene_list) {
  exon_positions <- getBM(attributes = c("exon_chrom_start", "exon_chrom_end"),
                          filters = "hgnc_symbol",
                          values = gene,
                          mart = ensembl)
  exon_lengths <- exon_positions$exon_chrom_end - exon_positions$exon_chrom_start + 1
  total_length <- sum(exon_lengths)
  coding_seq_length_list[[gene]] <- total_length
  cat("Total length of protein-coding exons in", gene, "is", total_length, "bp.")}
                     
library(readxl)
#already combine the gene length with clean data
clean_with_lenth <- read_excel("C:/Users/panbo/OneDrive/桌面/project2/clean_with_lenth.xlsx")
clean_with_lenth = clean_with_lenth[! clean_with_lenth$gene_length == 0,]
clean_with_lenth$gene_length = as.numeric(clean_with_lenth$gene_length)
clean_with_lenth$mis_normal = clean_with_lenth$Missense_Mutation / clean_with_lenth$gene_length
clean_with_lenth$non_normal = clean_with_lenth$Nonsense_Mutation / clean_with_lenth$gene_length
clean_with_lenth$splice_normal = clean_with_lenth$Splice_Site / clean_with_lenth$gene_length
clean_with_lenth$silent_normal = clean_with_lenth$Silent/ clean_with_lenth$gene_length
                     
                     
                     
#reference for dndscv
clean_data_ref <- my_data[,c("patient_id", "Chromosome","Start_Position", "Reference_Allele","Tumor_Seq_Allele2","Variant_Type")]
clean_data_ref <- clean_data_ref[clean_data_ref$Variant_Type == "SNP",]
drop <- c('Variant_Type')
clean_data_ref <- clean_data_ref[,!(names(clean_data_ref) %in% drop)]
colnames(clean_data_ref) <- c("sampleID", "chr", "pos","ref","mut")
dndsout = dndscv(clean_data_ref)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)

#significant genes from reference
signif_genes = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes) 
                     
#get dnds from reference                    
sel_cv$dnds = (sel_cv$n_mis + sel_cv$n_non + sel_cv$n_spl) / sel_cv$n_syn
sel_cv_clean = sel_cv[c('gene_name', 'dnds')]
sel_cv_clean = sel_cv_clean[!sel_cv_clean$dnds == "Inf",]
any(is.na(sel_cv_clean$dnds))
sel_cv_clean <- na.omit(sel_cv_clean, cols = "dnds")
sel_cv_clean = sel_cv_clean[sel_cv_clean$gene_name %in% clean_summary$Hugo_Symbol, ]

# selection type count for reference
count_neg = sum(sel_cv_clean$dnds < 1)
count_neg
count_neutral = sum(sel_cv_clean$dnds == 1)
count_neutral
count_pos = sum(sel_cv_clean$dnds > 1)
count_pos
                     
