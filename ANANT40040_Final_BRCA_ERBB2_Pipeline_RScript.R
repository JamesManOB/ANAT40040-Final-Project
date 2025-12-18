
#ANAT40040 Assignment 2 - BRCA ERBB2 analysis
#0 load libraries
library(BiocManager)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(pheatmap)
library(ggplot2)
library(survival)
library(glmnet)

#1 Untar and load clinical, RNA, CNA (cBioPortal)
downloads_path = #secret
tar_filename = "brca_tcga_pan_can_atlas_2018.tar.gz"
file_path = paste(downloads_path, tar_filename, sep = "/")
untar(file_path)
#1.1 folder name after untar
folder_path = paste0(getwd(), "/brca_tcga_pan_can_atlas_2018")
#quick check
print(folder_path)
list.files(folder_path)

#2 Clinical Data
data_patient_path = paste(folder_path, "data_clinical_patient.txt", sep = "/")
data_patient_raw = read.delim(data_patient_path, stringsAsFactors =  FALSE,
                              check.names = FALSE)
data_patient = data_patient_raw[5:nrow(data_patient_raw), ]

#Patient Identifier
col_patient_id = 1
pat_ids = data_patient[, col_patient_id]

#survival columns
col_overallsurvival_status = which(colnames(data_patient) == "Overall Survival Status")
col_overallsurvival_time = which(colnames(data_patient) == "Overall Survival (Months)")

#################3 RNA Sequence data#######################
path_RNA = paste(folder_path, "data_mrna_seq_v2_rsem.txt", sep = "/")
data_RNA = read.delim(path_RNA, stringsAsFactors = FALSE,
                      check.names = FALSE)
#pattern in practical: col1 = gene ID, col2 = gene symbol, rest are samples
assay = round(as.matrix(data_RNA[,-c(1,2)]))
rownames(assay) = data_RNA[,1]

#clean counts
assay[is.na(assay)] = 0
assay[assay < 0] = 0

#filter genes woth too few counts
smallestgroupsize = 3
keep = rowSums(assay >= 10) >= smallestgroupsize
assay = assay[keep,]

############4 CNA data + ERBB2 amplification#############

path_cna = paste(folder_path, "data_cna.txt", sep = "/")
data_cna = read.delim(path_cna, stringsAsFactors = FALSE,
                      check.names = FALSE)


#build ERBB2 amplification metadata
erbb2_row = which(data_cna$Hugo_Symbol == "ERBB2")
if (length(erbb2_row) !=1){
stop("Could not find ERBB2 in data_cna, check Hugo_Symbol column")}

#extract ERBB2 cna values across samples - drop both gene ID columns
erbb2_cna = as.numeric(data_cna[erbb2_row, -c(1,2)])
names(erbb2_cna) = colnames(data_cna)[-c(1,2)]

#################5 Align samples between RNA(assay) and CNA(erbb2_cna)#######
#sample IDs from rna and cna
rna_samples = colnames(assay)
cna_samples = names(erbb2_cna)
            
#intersection of samples to use for joint analysis
common_samples = intersect(rna_samples, cna_samples)
length(common_samples)
head(common_samples)

#restrict both matrices/vectores to shared samples
assay = assay[, common_samples]
erbb2_cna = erbb2_cna[common_samples]

#quick check
dim(assay)
length(erbb2_cna)

########6 ERBB2 amplification status (cna > 0 = amplified)#####
erbb2_status = ifelse(erbb2_cna>0, "Amplified","NotAmplified")
erbb2_status = factor(erbb2_status,
                      levels = c("NotAmplified", "Amplified"))
table(erbb2_status)

############7 Match smaples to clinical data - overall survival######

#function standardise tcga barcodes
standardise_barcode = function(x) { x = gsub("\\.","-",x)
substr(x,1,12)}

#patient ids and barcodes
pat_ids = data_patient[, col_patient_id]
pat_barcodes = standardise_barcode(pat_ids)

#sample barcodes in same order as clmns of assay
sample_barcodes = standardise_barcode(colnames(assay))

#vectors survival time and status
OS_time = rep(NA_real_, length(sample_barcodes))
OS_status = rep(NA_integer_, length(sample_barcodes))
                
for(i in seq_along(sample_barcodes)) {bc = sample_barcodes[i]
idx = which(pat_barcodes == bc)
#os time months
if(length(idx) == 1){OS_time[i] = as.numeric(data_patient[idx, col_overallsurvival_time])
#os status 1= deceased 0= living
status_str = data_patient[idx, col_overallsurvival_status]
OS_status[i] = ifelse(status_str == "1:DECEASED", 1L, 0L)}}

#checks
summary(OS_time)
table(OS_status, useNA = "ifany")

######8 Build metadata for DESeq2####
metadata = data.frame(erbb2_status = erbb2_status, OS_time = OS_time, 
                      OS_status = OS_status, row.names = colnames(assay))
###check rownames match column names assay
all(rownames(metadata) == colnames(assay))
#####9 DESeq2 Differential expression (Amplified vs Not Amplified)#####

dds = DESeqDataSetFromMatrix(
  countData = assay,
  colData = metadata,
  design  = ~ erbb2_status)

dds = DESeq(dds)

# Differential expression results: Amplified vs NotAmplified
res = results(dds, contrast = c("erbb2_status", "Amplified", "NotAmplified"))

# Order by adjusted p-value
res_ordered = res[order(res$padj),]

# Top 10 differentially expressed genes
top10_DE = res_ordered[1:10,]

#making a table to be exported
top10_table = as.data.frame(top10_DE)
top10_table$Gene = rownames(top10_DE)
top10_table = top10_table[, c("Gene","baseMean",
                              "log2FoldChange", "lfcSE",
                              "stat","pvalue","padj")]

#out_dir = "C:/Users/James/OneDrive - University College Dublin/Desktop/ANAT40040 - Materials/Final_CBioPortal_Assignment/ANAT40040_Final_Project"
#out_path = file.path(out_dir, "Top10_DE_ERBB2_Amplified_vs_NotAmplified.csv")
#write.csv(top10_table, file = out_path, row.names = FALSE)
#commented out but can be changed for other computers

##########10 Pathway Enrichment################
#selecting signifcant genes - keep genes with padj < 0.05 and non-missing padj
res_sig = res[!is.na(res$padj) & res$padj < 0.05, ]

#split - over/under expressed in ERBB2 amplified tumours
DE_over = rownames(res_sig[res_sig$log2FoldChange > 0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange < 0,])

length(DE_over)
length(DE_under)

##GO biological processes enrichment (clusterprofile::enrichGO)

go_over = enrichGO(gene = DE_over,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

go_under = enrichGO(gene = DE_under,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

#GO terms
head(go_over)
head(go_under)

#example plots
dotplot(go_over, showCategory = 10) + ggtitle("GO BP enrichment ERBB2 Amplified (Over Expressed)")
dotplot(go_under, showCategory = 10) + ggtitle("GO BP enrichment ERBB2 Amplified (Under expressed)")

#Mapping of gene symbols to Entrez IDs for Kegg+Reactome
gene_entrez_over = bitr(DE_over,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
gene_entrez_under = bitr(DE_under,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

#did it map though
nrow(gene_entrez_under)
nrow(gene_entrez_over)

#KEGG pathway enrichment
kegg_over = enrichKEGG(gene = gene_entrez_over$ENTREZID,
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
kegg_under = enrichKEGG(gene = gene_entrez_under$ENTREZID,
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

head(kegg_over)
head(kegg_under)

#KEGG plots
dotplot(kegg_over,  showCategory = 10) + ggtitle("KEGG enrichment – ERBB2 Amplified (Over expressed)")
dotplot(kegg_under, showCategory = 10) + ggtitle("KEGG enrichment – ERBB2 Amplified (Under expressed)")

#reactome pathway enrichment
reactome_over = enrichPathway(gene = gene_entrez_over$ENTREZID,
                              organism = "human",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
reactome_under = enrichPathway(gene = gene_entrez_under$ENTREZID,
                               organism = "human",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)
head(reactome_over)
head(reactome_under)

dotplot(reactome_over,  showCategory = 10) + ggtitle("Reactome enrichment – ERBB2 Amplified (Over-expressed)")
dotplot(reactome_under, showCategory = 10) + ggtitle("Reactome enrichment – ERBB2 Amplified (Under-expressed)")


#############11 Variance Stabilisation VST + PCA Plotting######

#vst transform
vsd = vst(dds, blind = FALSE)
vsd_mat = assay(vsd)

#pca coloured by erbb2
pca_plot = plotPCA(vsd, intgroup = "erbb2_status")
pca_plot

#######12 Heatmap of vst for de genes#####
#taking the top 100 DE genes by adjusted p value because my laptop did not like
#the fact I was trying to use all DE genes
res_sig_ordered = res_sig[order(res_sig$padj), ]
top_genes = rownames(res_sig_ordered)[1:100]

#subset the vst matrx to those genes
vsd_DE = vsd_mat[top_genes,]

#annotation frame
ann_col = data.frame(ERBB2 = metadata$erbb2_status)
rownames(ann_col) = colnames(vsd_DE)

#heatmap
pheatmap(vsd_DE, scale = "row",
         annotation_col = ann_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Top 100 DE genes ERBB Amplified vs Not Amplified")

##########13 Lasso Cox Regression ANALYSISUSUS######

#I believe in having all packages at the start of the script
#but here are the the packages needed
#library(glmnet)
#library(survival)

res_df = as.data.frame(res)
res_sig = subset(res_df, !is.na(padj) & padj < 0.05)

genes_cox = rownames(res_sig)
#as a check
genes_cox = intersect(genes_cox, rownames(vsd_mat))

#how many genes go into model?
length(genes_cox)

#build expression matx for these genes
expr_cox = t(vsd_mat[genes_cox, ]) #rows are samples and cols are genes

#samples with complete survival data
keep_surv = which(!is.na(metadata$OS_time) & !is.na(metadata$OS_status))
expr_cox_cox = expr_cox[keep_surv,]
meta_cox = metadata[keep_surv,]

#survival object
#positive only cox model??????????
meta_cox$OS_time_adj = meta_cox$OS_time
meta_cox$OS_time_adj[meta_cox$OS_time_adj <= 0] = 0.1

y_surv = Surv(time = meta_cox$OS_time_adj, event = meta_cox$OS_status)

#fit lasso cox model with cv + reproducibility
set.seed(123)

cvfit = cv.glmnet(x = expr_cox_cox,
                  y = y_surv,
                  family = "cox",
                  alpha = 1,
                  nfolds = 10)

#cv curve??
plot(cvfit)

lambda_min = cvfit$lambda.min
lambda_1se = cvfit$lambda.1se
lambda_min
lambda_1se

#extraction of genes selected at lamda.min
coef_min = coef (cvfit, s = "lambda.min")
active_idx = which(coef_min != 0)
genes_selected = rownames(coef_min)[active_idx]

cox_gene_table = data.frame(Gene = genes_selected,
                            Coeffiecient = as.numeric(coef_min[active_idx]))

#ordering by absolute effect size#######################
cox_gene_table = cox_gene_table[order(abs(cox_gene_table$Coeffiecient),decreasing = TRUE), ]

cox_gene_table   #survival signature table

#risk score for each patient
risk_score = as.numeric(predict(cvfit, newx = expr_cox_cox,
                                s    = "lambda.min",
                                type = "link"))

meta_cox$risk_score = risk_score

# Split into high/low risk by median
cutoff = median(risk_score, na.rm = TRUE)

meta_cox$risk_group = ifelse(risk_score > cutoff, "High", "Low")
meta_cox$risk_group = factor(meta_cox$risk_group,levels = c("Low", "High"))

table(meta_cox$risk_group)

# Cox model using risk score (for loglikelihood + C index)
cox_model = coxph(y_surv ~ risk_score, data = meta_cox)
summary(cox_model)

# Log likelihood
logLik(cox_model)

# Concordance index
summary(cox_model)$concordance  # first element = C index

#Kaplan/Meier curves for High vs Low risk
fit_km = survfit(y_surv ~ risk_group, data = meta_cox)

plot(fit_km, col = c("blue","red"),xlab = "Overall survival (months)",ylab = "Survival probability",
     main = "Kaplan-Meier curves by LASSO Cox risk group",lwd  = 2)

legend("bottomleft",col = c("blue","red"),legend = levels(meta_cox$risk_group),lwd = 2, bty = "n")


