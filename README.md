# Differential Expression Analysis for TMT Protein Data with EdgeR

This repo contains code for conducting a differential expression analysis of proteins from _F. cylindrus_ cultures grown under high temperature treatments and vitamin B<sub>12</sub> starvation. Thanks to the brilliant Scott McCain for the skeleton that this tutorial is based on, which can be found [here](https://github.com/jspmccain). 

The goal of this statistical analysis is to determine which proteins are differentially expressed in the treatment cultures when compared to control. 
There are 3 files in the repository - a script with the analysis (), the data (), and this short walkthrough guide. At this point, I assumed your data has been normalized and all missing values have been removed or imputed. Check out [this repo](https://github.com/pwilmart/IRS_normalization) for a tutorial on the IRS normalization procedure I used in the data. You can reference the [full edgeR user guide](http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) for an under-the-hood look at the workings of these functions. egdeR offers a few different types of analysis, but we will use their quasi-likelihood F-tests in this example. 

There are three main steps for completeting the analysis. 


## 1. "Massaging" the Data

We need to get the data into an appropriate format before using EdgeR's function to identify lists of differentially expressed proteins for each comparison we make. This is all working towards putting your protein abundances, treatments, and protein names into something called a   `DGE_list` object that we can conduct the DE analysis on. 

First, call the required packages to your environment. If you're installing  edgeR for the first time, can be done through Bioconductor. [Instructions here!](https://bioconductor.org/packages/release/bioc/html/edgeR.html) 

``` 
library(edgeR) 
library(statmod)
library(tidyverse)
``` 

``` 
# read protein file as .csv
prot_data <- read.csv('PD1_Norm_11072021.csv')

# Select only the columns needed (note I exclude some treatments here)
DE_data <- prot_data %>% select(
, B12_4_1_A, B12_4_2_B, B12_4_3_B, noB12_4_1_A, noB12_4_2_B, B12_12_1_B, B12_12_3_A, noB12_12_2_B, noB12_12_3_A)

```   
Then, we want to produce a data table with your treatment names as column names and only protein counts 

```
# Remove protein descriptions to include only counts
prot_counts <- DE_data[,2:length(DE_data)]
``` 

Next, we will list the treatment group type of each of your treatment replicates. Ex: my first three samples are control temperature, B12 replete treatments.

```
# Make vector with treatment groups
DE_data_groups <- factor(c('B12_4',	'B12_4','B12_4', 'noB12_4',	'noB12_4','B12_12', 'B12_12', 'noB12_12', 'noB12_12'))
``` 


  
## 2. The Main Event
All of this goes into a DGE_list object as below. The `counts` argument accepts the table with just protein abundances, the `group` argument takes the group type assignments we just made, and the `genes` argument is for a list of your protein names or id's (I used the gene product ID's found in the "accession" column; all that matters is that these are unqiue). It is a good idea to check for duplicates in what you're feeding to the `gene` argument here just to be safe. 

``` 
# Create DGEList object with protein counts, groups, and accession
DGE_list <- DGEList(counts = prot_counts, 
                        group = DE_data_groups,
                        genes = DE_data$accession)
``` 


Then, we will make the design matrix. This outlines the structure of the linear model we will be testing on our proteins to see if they are differentially expressed.

``` 
# Create design matrix 
design.mat <- model.matrix(~ 0 + DE_data_groups)

``` 

Along with the design matrix, we need to estimate tag dispersal of the data. This is done with the `estimateDisp` function. 

```
# Get tagwise dispersion of tags
DGE_disp <- estimateDisp(DGE_list, design.mat, robust = TRUE)
```
 Lastly, we fit these to a general linear model. 
 
```
 # Fit to GLM
fit <- glmQLFit(DGE_disp, design.mat, robust = TRUE)

```

Here we lay out the comparisons needed to tell us which proteins are differentially expressed per treatment. Which treatments you choose to compare depends on your question. Here I make three comparisons: one which tells me about the effect: (1) of B<sub>12</sub> starvation, (2) elevated temperature, and (3) the combination of the two stressors. Then, fit tests are conducted on each comparison. The `topTags` step takes the results of the fit tests and produces a list of the top differntially expressed protein per pairwise comparison (`n = 500`). These are produced in the topTags class and so the values need to be extracted. 


```
# Make pairwise comparisons

# B12_4 vs noB12_4
contrast_4_noB12 <- makeContrasts( DE_data_groupsB12_4 - DE_data_groupsnoB12_4, levels = design.mat)

# B12_4 vs B12_12
contrast_12_B12 <- makeContrasts( DE_data_groupsB12_4 - DE_data_groupsB12_12, levels = design.mat)

# B12_4 vs noB12_12
contrast_12_noB12 <- makeContrasts( DE_data_groupsB12_4 - DE_data_groupsnoB12_12, levels = design.mat)

 # Fit tests
 
# B12_4 vs B12_12
qlf_12_B12 <- glmQLFTest(fit, contrast = contrast_12_B12)

# B12_4 vs noB12_4
qlf_4_noB12 <- glmQLFTest(fit, contrast = contrast_4_noB12)

# B12_4 vs noB12_12
qlf_12_noB12 <- glmQLFTest(fit, contrast = contrast_12_noB12)


# Export list of DE proteins per comparison

# B12_4 vs B12_12
hits_12_B12 <- topTags(qlf_12_B12, n = 500)

# B12_4 vs noB12_4
hits_4_noB12 <- topTags(qlf_4_noB12, n = 500)

# B12_4 vs noB12_12
hits_12_noB12 <- topTags(qlf_12_noB12, n = 500)

# Extract dataframes
hits_12_B12 <- hits_12_B12[[1]][,]
hits_4_noB12 <- hits_4_noB12[[1]][,]
hits_12_noB12 <- hits_12_noB12[[1]][,]
```

A bit of cleanup will make the data a little easier to look at. 

```
# Change genes column to say accession
hits_12_B12 <- hits_12_B12 %>% 
  rename(accession = genes)
hits_4_noB12 <- hits_4_noB12 %>% 
  rename(accession = genes)
hits_12_noB12 <- hits_12_noB12 %>% 
  rename(accession = genes)

# Add back descriptions
hits_12_B12 <- merge(hits_12_B12, prot_data, by="accession")
hits_4_noB12 <- merge(hits_4_noB12, prot_data, by="accession")
hits_12_noB12 <- merge(hits_12_noB12, prot_data, by="accession")

# Trim dataset to only DE info
hits_12_B12 <- hits_12_B12 %>% select(accession, Description, logFC, logCPM, F, PValue, FDR)
hits_4_noB12 <- hits_4_noB12 %>% select(accession, Description, logFC, logCPM, F, PValue, FDR)
hits_12_noB12 <- hits_12_noB12 %>% select(accession, Description, logFC, logCPM, F, PValue, FDR)

# Order by PValue 
hits_12_B12 <- hits_12_B12 %>% arrange(PValue)
hits_4_noB12 <- hits_4_noB12 %>% arrange(PValue)
hits_12_noB12 <- hits_12_noB12 %>% arrange(PValue)
```




## 3. Making Sense of It All 

If done correctly, this is what the final data should look like. In the comparisons we've made, proteins with a positive `logFC` corresponds to an enriched protein compared to the control. 
| accession  | Description                                                                                                 | logFC      | logCPM     | F          | PValue     | FDR        |
| ---------- | ----------------------------------------------------------------------------------------------------------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| OEU11144.1 | 5-methyltetrahydropteroyltriglutamate--homocysteine methyltransferase \[Fragilariopsis cylindrus CCMP1102\] | 5.27942766 | 11.133677  | 84.0634102 | 3.85E-06   | 0.00522224 |
| OEU10229.1 | protoporphyrin IX Mg-chelatase subunit D \[Fragilariopsis cylindrus CCMP1102\]                              | 1.98202135 | 8.55656912 | 17.9758636 | 0.00167818 | 0.83348938 |
| OEU11214.1 | hypothetical protein FRACYDRAFT\_246327 \[Fragilariopsis cylindrus CCMP1102\]                               | 2.4411276  | 9.36918024 | 17.6200734 | 0.00184264 | 0.83348938 |
| OEU18387.1 | P-ATPase family transporter: zinc/lead/cadmium/mercury ion \[Fragilariopsis cylindrus CCMP1102\]            | 1.21189564 | 9.7739918  | 15.05875   | 0.00299821 | 0.99994004 |
| OEU08987.1 | T-complex protein 1 subunit gamma \[Fragilariopsis cylindrus CCMP1102\]                                     | 1.01189583 | 10.2914415 | 11.8166049 | 0.00626436 | 0.99994004 |
