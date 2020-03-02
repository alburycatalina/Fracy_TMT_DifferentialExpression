# Robinson 2010
library(edgeR)

setwd("~/Desktop/OneDrive - Dalhousie University/MSc/CHI/FC_Global_Proteome/Methods/FC_GP_DE/DE")

# DE Analysis for vitamin B12 ----------------

# Upload protein file that includes accessions and descriptions of B12 treatment comparisons
DE_data_vit <- read.csv('DE_vit_130119.csv')
prot_counts_vit <- DE_data_vit[,2:length(DE_data_vit)]
prot_id_vit <- make.unique(as.character(DE_data_vit$protein_id), sep = '.')

# Make vector with treatment groups
DE_data_groups_vit <- factor(c('B12_4',	'noB12_4',	'noB12_4',	'B12_4'))

# Create DGEList object with protein counts, groups, and protein ID's
DGE_list_vit <- DGEList(counts = prot_counts_vit, 
                   group = DE_data_groups_vit,
                   genes = DE_data_vit$protein_id)

# Create design matrix 
design.mat <- model.matrix(~ 0 + DE_data_groups_vit)

# Get tagwise dispersion of tags
DGE_disp <- estimateDisp(DGE_list_vit, design.mat, robust = TRUE)

# Fit to GLM
fit <- glmQLFit(DGE_disp, design.mat, robust = TRUE)

# Make pairwise comparisons
B12_contrast <- makeContrasts( DE_data_groups_vitB12_4 - DE_data_groups_vitnoB12_4, levels = design.mat)

# Fit test
qlf_1 <- glmQLFTest(fit, contrast = B12_contrast)

# Generate fire top hits mixtape with top 50 differentially expressed proteins
vit_hits <- topTags(qlf_1, n = 500)

# Export greatest hits mixtape to csv 
write.csv(vit_hits, file = "vit_hits_500.csv")





# DE Analysis for temperature 10----------------------

# Upload protein file that includes accessions and descriptions of B12 treatment comparisons
DE_data_temp <- read.csv('DE_temp_130119.csv')

# Index dataset to get only counts (remove protein descriptions)
prot_counts_temp <- DE_data_temp[,2:length(DE_data_temp)]

prot_id_temp <- make.unique(as.character(DE_data_temp$protein_id), sep = '.')

# Make vector with treatment groups
DE_data_groups_temp <- factor(c('B12_4',	'B12_10',	'B12_12',	'B12_12',	'B12_10',	'B12_4'))

# Create DGEList object with protein counts, groups, and protein ID's
DGE_list_temp <- DGEList(counts = prot_counts_temp, 
                        group = DE_data_groups_temp,
                        genes = DE_data_temp$protein_id)

# Create design matrix 
design.mat <- model.matrix(~ 0 + DE_data_groups_temp)

# Get tagwise dispersion of tags
DGE_disp <- estimateDisp(DGE_list_temp, design.mat, robust = TRUE)


# Fit to GLM
fit <- glmQLFit(DGE_disp, design.mat, robust = TRUE)

# Make pairwise comparisons
temp_contrast <- makeContrasts(DE_data_groups_tempB12_4 - DE_data_groups_tempB12_10, levels = design.mat)

# Fit test
qlf_1 <- glmQLFTest(fit, contrast = temp_contrast)

# Generate fire top hits mixtape with top 50 differentially expressed proteins
temp_hits_10 <- topTags(qlf_1, n = 500)

# Export greatest hits mixtape to csv 
write.csv(temp_hits_10, file = "temp_hits_500_10.csv")


# DE Analysis for temperature 12----------------------

# Upload protein file that includes accessions and descriptions of B12 treatment comparisons
DE_data_temp <- read.csv('DE_temp_130119.csv')

# Index dataset to get only counts (remove protein descriptions)
prot_counts_temp <- DE_data_temp[,2:length(DE_data_temp)]

prot_id_temp <- make.unique(as.character(DE_data_temp$protein_id), sep = '.')

# Make vector with treatment groups
DE_data_groups_temp <- factor(c('B12_4',	'B12_10',	'B12_12',	'B12_12',	'B12_10',	'B12_4'))

# Create DGEList object with protein counts, groups, and protein ID's
DGE_list_temp <- DGEList(counts = prot_counts_temp, 
                         group = DE_data_groups_temp,
                         genes = DE_data_temp$protein_id)

# Create design matrix 
design.mat <- model.matrix(~ 0 + DE_data_groups_temp)

# Get tagwise dispersion of tags
DGE_disp <- estimateDisp(DGE_list_temp, design.mat, robust = TRUE)


# Fit to GLM
fit <- glmQLFit(DGE_disp, design.mat, robust = TRUE)

# Make pairwise comparisons
temp_contrast <- makeContrasts(DE_data_groups_tempB12_4 - DE_data_groups_tempB12_12, levels = design.mat)

# Fit test
qlf_1 <- glmQLFTest(fit, contrast = temp_contrast)

# Generate fire top hits mixtape with top 50 differentially expressed proteins
temp_hits_12 <- topTags(qlf_1, n = 500)

# Export greatest hits mixtape to csv 
write.csv(temp_hits_12, file = "temp_hits_500_12.csv")



