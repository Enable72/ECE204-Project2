library("dndscv")
data("dataset_simbreast", package="dndscv")
dndsout = dndscv(mutations)

head(mutations)

my_data <- read.delim("/Users/atishnasamantaray/Downloads/Team_2_HNSC/TCGA.HNSC.mutations.txt")
clean_data <- my_data[,c("patient_id", "Chromosome","Start_Position", "Reference_Allele","Tumor_Seq_Allele2","Variant_Type")]
clean_data <- clean_data[clean_data$Variant_Type == "SNP", ]
drop <- c('Variant_Type')
clean_data <- clean_data[,!(names(clean_data) %in% drop)]
colnames(clean_data) <- c("sampleID", "chr", "pos","ref","mut")
clean_data <- clean_data[!clean_data$ref == "-", ]
clean_data <- clean_data[!clean_data$mut == "-", ]
dndsout = dndscv(clean_data)
dndsout = dndscv(mutations)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)