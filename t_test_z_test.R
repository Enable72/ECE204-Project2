library("dndscv")
data("dataset_simbreast", package="dndscv")

my_data <- read.delim("TCGA.HNSC.mutations.txt")
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

clean_summary$Nonsyno = clean_summary$Nonsense_Mutation + clean_summary$Missense_Mutation + clean_summary$Splice_Site
clean_summary$dnds = clean_summary$Nonsyno / clean_summary$Silent

# t test z test
library(BSDA)
dnds_res = clean_summary$Nonsyno / clean_summary$Silent
dnds_res <- dnds_res[!is.na(dnds_res) & !is.infinite(dnds_res)]
dnds_ref = (sel_cv$n_mis + sel_cv$n_non + sel_cv$n_spl) / sel_cv$n_syn 
dnds_ref <- dnds_ref[!is.na(dnds_ref) & !is.infinite(dnds_ref)]
t.test(dnds_res, dnds_ref)
z.test(dnds_res, dnds_ref, sigma.x=sd(dnds_res), sigma.y=sd(dnds_ref))
