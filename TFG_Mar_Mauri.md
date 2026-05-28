Molecular composition of entorhinal synapses in Alzheimer’s disease
================
Mar
2026-05-28

``` r
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(openxlsx)
```

## 1. Library & data

``` r
library(MSstats)
library(dplyr)
library(stringr)
library(tidyverse)

feature_data <- read.csv("181120-MSStats-input-refined-fromCRG.txt", header = TRUE, sep = ",")
```

## 2. Remove contaminants

``` r
if("Proteins" %in% colnames(feature_data)) {
  feature_data <- feature_data %>%
    filter(!str_detect(Proteins, "CON_"))
} else {
  feature_data <- feature_data %>%
    filter(!str_detect(ProteinName, "CON_"))
}
```

## 3. Filter ≥2 unique peptides per protein

``` r
peptide_counts <- feature_data %>%
distinct(ProteinName, PeptideSequence) %>%
group_by(ProteinName) %>%
summarise(n_unique_peptides = n())

proteins_keep_pep <- peptide_counts %>%
filter(n_unique_peptides >= 2) %>%
pull(ProteinName)

feature_data_filt1 <- feature_data %>%
filter(ProteinName %in% proteins_keep_pep)
```

## 4. Filter ≥50% presence in ≥1 group

We are treating CN and PART as a single group; therefore, we must change
the settings so they are read as one condition. To do this, we will
first create a combined group where CN and PART are treated identically,
allowing us to calculate the presence of peptides within both. It is
essential to modify the Condition column directly rather than creating a
new one (such as condition2 for AD vs CN-PART), because MSstats will
automatically default to the original Condition column when processing
the data.

``` r
table(feature_data_filt1$Condition, useNA = "always")
```

    ## 
    ##      AD Control    PART    <NA> 
    ##  152136   95085   95085       0

``` r
feature_data_filt1 <- feature_data_filt1 %>%
  mutate(
    Condition = ifelse(Condition %in% c("Control","PART"), "CN_PART", "AD"))
 
     table(feature_data_filt1$Condition)
```

    ## 
    ##      AD CN_PART 
    ##  152136  190170

``` r
     # Num of cases in each group
addmargins(table(feature_data_filt1$Condition, useNA = "always"))
```

    ## 
    ##      AD CN_PART    <NA>     Sum 
    ##  152136  190170       0  342306

``` r
feature_data_filt1 %>% count(Condition)
```

    ##   Condition      n
    ## 1        AD 152136
    ## 2   CN_PART 190170

``` r
table(feature_data_filt1$Condition, useNA = "always")
```

    ## 
    ##      AD CN_PART    <NA> 
    ##  152136  190170       0

Once the condition is modified , the filter can be applied:

``` r
# Count total samples per group
group_sizes <- feature_data_filt1 %>%
  distinct(BioReplicate, Condition) %>%
  count(Condition, name = "n_group")

# Identify protein presence per sample
protein_sample_presence <- feature_data_filt1 %>%
  filter(!is.na(Intensity)) %>%
  distinct(ProteinName, BioReplicate, Condition)

# Count per group
presence_counts <- protein_sample_presence %>%
  count(ProteinName, Condition, name = "n_present") %>%
  left_join(group_sizes, by = "Condition") %>%
  mutate(prop_present = n_present / n_group)

# Keep proteins present in ≥50% of at least one group
proteins_keep_presence <- presence_counts %>%
  filter(prop_present >= 0.5) %>%
  distinct(ProteinName) %>%
  pull(ProteinName)

feature_data_filtered <- feature_data_filt1 %>%
  filter(ProteinName %in% proteins_keep_presence)
```

## 5. Report filtering

``` r
cat("Proteins after filtering:",
length(unique(feature_data$ProteinName)),"\n")
```

    ## Proteins after filtering: 3170

``` r
cat("After ≥2peptides filter :",
length(unique(feature_data_filt1$ProteinName)),"\n")
```

    ## After ≥2peptides filter : 2291

``` r
cat("After 50% peptides filter :",
length(unique(feature_data_filtered$ProteinName)),"\n")
```

    ## After 50% peptides filter : 2291

## 6. Log2 transform+equalizeMedians normalization

``` r
processed <- dataProcess(
  raw = feature_data_filtered,
  normalization = "equalizeMedians",
  summaryMethod = "TMP",
  censoredInt = "NA",
  MBimpute = FALSE,
  verbose = FALSE
)
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

## 7. Replace FeatureLevelData Inside MSsatats object

``` r
feature_norm <- processed$FeatureLevelData
```

## 8. Custom Synaptic Normalization (FeatureLevel)

To achieve this, we will use the peptide intensities of the 60 proteins
annotated to the enriched pathways of presynaptic active zone or
post-synaptic density (see attached .txt files). For each sample:
Synaptic factorj = median intensity of 60 synaptic proteins in sample j
. Then subtract that value from all peptide intensities in that sample:

● Adjusted intensityij = Intensityij − Synaptic factorj ● Because
MSstats works in log2 scale, subtraction corresponds to division in
linear space.

``` r
 #1. Charge normalizers
synaptic_proteins <- read.delim("SynapticNormalizers.txt", header = FALSE, sep = "\t")$V2

colnames(feature_norm)
```

    ##  [1] "PROTEIN"      "PEPTIDE"      "TRANSITION"   "FEATURE"      "LABEL"       
    ##  [6] "GROUP"        "RUN"          "SUBJECT"      "FRACTION"     "originalRUN" 
    ## [11] "censored"     "INTENSITY"    "ABUNDANCE"    "newABUNDANCE"

``` r
feature_norm <- feature_norm %>%
  mutate(Intensity = as.numeric(INTENSITY))

# Calculate synaptic median per sample
synaptic_factors <- feature_norm %>%
  filter(PROTEIN %in% synaptic_proteins) %>%
  group_by(originalRUN) %>%
  summarise(
    syn_median = median(Intensity, na.rm = TRUE)
  )

#Optional quality check:
print(synaptic_factors)
```

    ## # A tibble: 18 × 2
    ##    originalRUN syn_median
    ##    <fct>            <dbl>
    ##  1 AD.11         38000000
    ##  2 AD.14         42000000
    ##  3 AD.15         48500000
    ##  4 AD.17         30500000
    ##  5 AD.2          20000000
    ##  6 AD.3          34000000
    ##  7 AD.6          30000000
    ##  8 AD.7          32000000
    ##  9 Control.10    43000000
    ## 10 Control.13    31500000
    ## 11 Control.18    32000000
    ## 12 Control.4     34000000
    ## 13 Control.8     41000000
    ## 14 PART.1        44500000
    ## 15 PART.12       64000000
    ## 16 PART.16       46000000
    ## 17 PART.5        41000000
    ## 18 PART.9        30000000

``` r
#You want reasonable n_syn across samples.
table(feature_norm$ProteinName %in% synaptic_proteins)
```

    ## < table of extent 0 >

``` r
#If many are missing in some samples, your normalization factor may become unstable.
#You may want to require at least e.g.:filter(n() >= 20)
#per sample before computing median. We do this pre-summarization because Synaptic abundance shift affects all peptides, adjustment must occur at feature level and otherwise protein summarization will propagate bias

library(dplyr)

synaptic_factors <- feature_norm %>%
  filter(PROTEIN %in% synaptic_proteins) %>%
  group_by(originalRUN) %>%       # group by sample
  filter(n() >= 20) %>%            # keep only replicates with ≥20 proteins
  summarise(
    syn_median = median(Intensity, na.rm = TRUE)
  )

print(synaptic_factors)
```

    ## # A tibble: 18 × 2
    ##    originalRUN syn_median
    ##    <fct>            <dbl>
    ##  1 AD.11         38000000
    ##  2 AD.14         42000000
    ##  3 AD.15         48500000
    ##  4 AD.17         30500000
    ##  5 AD.2          20000000
    ##  6 AD.3          34000000
    ##  7 AD.6          30000000
    ##  8 AD.7          32000000
    ##  9 Control.10    43000000
    ## 10 Control.13    31500000
    ## 11 Control.18    32000000
    ## 12 Control.4     34000000
    ## 13 Control.8     41000000
    ## 14 PART.1        44500000
    ## 15 PART.12       64000000
    ## 16 PART.16       46000000
    ## 17 PART.5        41000000
    ## 18 PART.9        30000000

``` r
# Optional: if u want to see the dropped replicates
feature_norm %>%
  filter(PROTEIN %in% synaptic_proteins) %>%
  group_by(originalRUN) %>%
  summarise(n_proteins = n()) %>%
  filter(n_proteins < 20)
```

    ## # A tibble: 0 × 2
    ## # ℹ 2 variables: originalRUN <fct>, n_proteins <int>

## 9. Subtract from all peptides:

``` r
library(dplyr)
feature_norm <- feature_norm %>%
  left_join(synaptic_factors, by = "originalRUN") %>%
  mutate(
    # SECURITY : If a sample does not have a synaptic factor, we keep its original intensity instead of converting it to NA
    Intensity = ifelse(is.na(syn_median), Intensity, Intensity - syn_median)
  ) %>%
  dplyr::select(-syn_median)
```

## 10. Replace inside MSstats object

``` r
processed$FeatureLevelData <- feature_norm

colnames(processed$FeatureLevelData)
```

    ##  [1] "PROTEIN"      "PEPTIDE"      "TRANSITION"   "FEATURE"      "LABEL"       
    ##  [6] "GROUP"        "RUN"          "SUBJECT"      "FRACTION"     "originalRUN" 
    ## [11] "censored"     "INTENSITY"    "ABUNDANCE"    "newABUNDANCE" "Intensity"

## 11. Run only summarization (Do NOT re-normalize):

``` r
# Rename to match the MSstats input format
data_for_summarization <- processed$FeatureLevelData %>%
  rename(
    ProteinName = PROTEIN ,
    PeptideSequence = PEPTIDE,
    #MSstats requires these two columns even if they are identical after the initial processing
    PeptideModifiedSequence = PEPTIDE, 
    PrecursorCharge = FEATURE, # Use FEATURE that combines peptide and charge
    IsotopeLabelType = LABEL,
    Condition = GROUP,
    BioReplicate = SUBJECT,
    Run = RUN,
    # Personalized'Intensity' (synaptic adjusted one)
    Intensity = Intensity 
  ) %>%
  mutate(
    FragmentIon = "NA",
    ProductCharge = "NA"
  )
```

``` r
processed_final <- dataProcess(
  raw = data_for_summarization,
  normalization = FALSE,
  summaryMethod = "TMP",
  censoredInt = "NA",
  MBimpute = TRUE,
  verbose = FALSE ) 
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%

    ##   |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%

    ##   |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

## 12. Extract summarized protein data:

``` r
protein_data <- processed_final$ProteinLevelData

head(protein_data$GROUP)
```

    ## [1] AD AD AD AD AD AD
    ## Levels: AD CN_PART

``` r
table(protein_data$GROUP)
```

    ## 
    ##      AD CN_PART 
    ##   10110   11671

``` r
table(protein_data$originalRUN)
```

    ## 
    ##    1   10   11   12   13   14   15   16   17   18    2    3    4    5    6    7 
    ## 1315 1155 1205 1240  989 1151 1308 1277  974 1155 1220 1158 1214 1307 1367 1253 
    ##    8    9 
    ## 1276 1217

``` r
table(distinct(protein_data, originalRUN, GROUP)$GROUP)
```

    ## 
    ##      AD CN_PART 
    ##       8      10

## 13. PCA

``` r
# Step 1: Read the untouched raw file again as an independent metadata dataframe
raw_metadata <- read.csv("181120-MSStats-input-refined-fromCRG.txt", header = TRUE, sep = ",") %>%
  # Keep only the sample column and the true original condition column
  # We use 'BioReplicate' as it maps to 'SUBJECT' in the summarized data
  dplyr::select(BioReplicate, Condition) %>%
  mutate(
    SUBJECT = str_trim(as.character(BioReplicate)),
    TRUE_GROUP = str_trim(as.character(Condition)) %>% str_to_upper()
  ) %>%
  # Standardize "CONTROL" string to "CN"
  mutate(TRUE_GROUP = ifelse(TRUE_GROUP == "CONTROL", "CN", TRUE_GROUP)) %>%
  dplyr::select(SUBJECT, TRUE_GROUP) %>%
  distinct()

# Step 2: Extract MSstats summarized data and map true groups using SUBJECT
pca_data_prep <- processed_final$ProteinLevelData %>%
  dplyr::select(Protein, originalRUN, SUBJECT, LogIntensities) %>%
  pivot_wider(
    names_from = Protein,
    values_from = LogIntensities
  ) %>%
  mutate(SUBJECT = str_trim(as.character(SUBJECT))) %>%
  # Inner join with our freshly read raw metadata
  inner_join(raw_metadata, by = "SUBJECT") %>%
  rename(GROUP = TRUE_GROUP) %>%
  dplyr::select(-SUBJECT)

# Step 3: Run group-specific robust imputation and drop zero-variance columns
pca_data <- pca_data_prep %>%
  group_by(GROUP) %>%
  mutate(
    across(
      -originalRUN,
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    )
  ) %>%
  ungroup() %>%
  # Global fallback imputation for group-specific NaNs
  mutate(
    across(
      -c(originalRUN, GROUP),
      ~ ifelse(is.nan(.) | is.na(.), mean(., na.rm = TRUE), .)
    )
  ) %>%
  # Safe Filtering: Keep only numeric columns with variance > 0 and no NAs
  dplyr::select(
    originalRUN, GROUP, 
    where(~ is.numeric(.) && !any(is.na(.)) && var(., na.rm = TRUE) > 0)
  )

# DIAGNOSTIC CHECK 

cat("--- DATASET DIMENSIONS FOR PCA ---\n")
```

    ## --- DATASET DIMENSIONS FOR PCA ---

``` r
cat("Number of Samples (Rows):", nrow(pca_data), "\n")
```

    ## Number of Samples (Rows): 18

``` r
cat("Number of Proteins (Columns):", ncol(pca_data) - 2, "\n") 
```

    ## Number of Proteins (Columns): 1545

``` r
cat("----------------------------------\n")
```

    ## ----------------------------------

``` r
# Step 4: Run the Principal Component Analysis (PCA)
pca_res <- prcomp(
  pca_data %>% dplyr::select(-originalRUN, -GROUP),
  scale. = TRUE
)

# Step 5: Duplicate PART rows to display them both standalone and inside NoNP
pca_coords <- as.data.frame(pca_res$x) %>%
  mutate(GROUP = pca_data$GROUP)

pca_plot <- bind_rows(
  # Subset A: AD and PART remain standalone, CN becomes NoNP
  pca_coords %>%
    mutate(
      VISUAL_GROUP = case_when(
        GROUP == "AD"   ~ "AD",
        GROUP == "PART" ~ "PART",
        GROUP == "CN"   ~ "NoNP",
        TRUE            ~ "CHECK_NAME"
      )
    ),
  
  # Subset B: Duplicate ONLY PART rows and force their visual group to NoNP
  pca_coords %>%
    filter(GROUP == "PART") %>%
    mutate(VISUAL_GROUP = "NoNP")
)

# Step 6: Generate the dynamic text values for the subtitle
n_ad   <- sum(pca_coords$GROUP == "AD")
n_part <- sum(pca_coords$GROUP == "PART")
n_cn   <- sum(pca_coords$GROUP == "CN")
n_nonp <- n_cn + n_part

# Step 7: Plotting (With updated legend title)
pca_plot_output <- ggplot(
  pca_plot,
  aes(x = PC1, y = PC2, color = VISUAL_GROUP)
) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "AD"   = "#F8766D",  # Coral Red
      "PART" = "#00BFC4",  # Turquoise
      "NoNP" = "#B254A5"   # Purple (CN + Duplicated PART)
    )
  ) +
  labs(
    title = "Proteomic PCA - Visualization Groups",
    subtitle = paste0(
      "AD: ", n_ad,
      " | PART: ", n_part,
      " | NoNP (CN+PART): ", n_nonp
    ),
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"),
    color = "Group"  # This changes the legend title from VISUAL_GROUP to Group
  )

# ==========================================
# SAVE THE IMAGE TO YOUR WORKSPACE
# ==========================================
ggsave(
  filename = "Proteomic_PCA_Plot.png", 
  plot = pca_plot_output,               
  device = "png",                       
  width = 8,                            
  height = 5,                           
  dpi = 300                             
)
```

## 14. Differential testing no covariables

``` r
levels(processed_final$ProteinLevelData$GROUP)
```

    ## [1] "AD"      "CN_PART"

``` r
# sale: "AD", "CN_PART

# AD vs CN_PART
contrast_matrix <- matrix(
  c(1, -1),  # AD vs CN_PART
  nrow = 1,
  byrow = TRUE
)

colnames(contrast_matrix) <- c("AD", "CN_PART")

rownames(contrast_matrix) <- "AD-CN_PART"

# Visualise
contrast_matrix
```

    ##            AD CN_PART
    ## AD-CN_PART  1      -1

``` r
comparison <- groupComparison(
  contrast.matrix = contrast_matrix,
  data = processed_final
)
```

    ## INFO  [2026-05-28 21:34:27]  == Start to test and get inference in whole plot ...
    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%

    ##   |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%

    ##   |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
    ## INFO  [2026-05-28 21:34:41]  == Comparisons for all proteins are done.

``` r
# Results
results <- comparison$ComparisonResult
head(results)
```

    ##   Protein      Label      log2FC        SE      Tvalue DF    pvalue adj.pvalue
    ## 1       0 AD-CN_PART  0.47040421 0.5764724  0.81600472 15 0.4272663  0.9999652
    ## 2  A0FGR8 AD-CN_PART  0.40987823 1.1416683  0.35901692  6 0.7318782  0.9999652
    ## 3  A0MZ66 AD-CN_PART  0.94752467 0.7731376  1.22555766 12 0.2438749  0.9999652
    ## 4  A6NCE7 AD-CN_PART -0.35729178 0.3296067 -1.08399439 15 0.2954839  0.9999652
    ## 5  A6NHL2 AD-CN_PART  0.04609414 0.6298085  0.07318755 16 0.9425640  0.9999652
    ## 6  D6RGH6 AD-CN_PART -0.03442352 1.2919909 -0.02664378  7 0.9794875  0.9999652
    ##   issue MissingPercentage ImputationPercentage
    ## 1  <NA>        0.52777778           0.47222222
    ## 2  <NA>        0.75000000           0.19444444
    ## 3  <NA>        0.61111111           0.38888889
    ## 4  <NA>        0.11111111           0.05555556
    ## 5  <NA>        0.08333333           0.08333333
    ## 6  <NA>        0.50000000           0.00000000

``` r
library(ggrepel)
# 1. Prepare data
plot_data <- comparison$ComparisonResult %>%
  mutate(
    log_p = -log10(pvalue),
    status = case_when(
      pvalue < 0.05 & log2FC > 1  ~ "Up-regulated",
      pvalue < 0.05 & log2FC < -1 ~ "Down-regulated",
      TRUE                        ~ "Not significant"
    )
  )

# 2. Top 10 Hits
top_10_hits <- plot_data %>% 
  filter(!is.na(pvalue)) %>% 
  arrange(pvalue) %>% 
  head(10)

library(openxlsx)
write.xlsx(top_10_hits, "Top_10_Proteins_Volcano_unadjusted.xlsx", rowNames = FALSE)


# 3. Volcano Plot 
volcano_final <- ggplot(plot_data, aes(x = log2FC, y = log_p, color = status)) +
  geom_point(alpha = 0.5, size = 1.2) + 
  scale_color_manual(values = c(
    "Down-regulated" = "#0072B2", 
    "Up-regulated" = "#D55E00", 
    "Not significant" = "gray85"
  )) +
  scale_y_continuous(limits = c(0, 17)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray40", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray40", size = 0.3) +
  
  geom_text_repel(
    data = top_10_hits,
    aes(label = Protein),
    size = 3.2, 
    color = "black", 
    fontface = "italic", 
    box.padding = 0.8,      
    segment.size = 0.2,
    max.overlaps = Inf,
    force = 2
  ) +
  labs(
    title = "Proteomic Dysregulation",
    subtitle = "Model unadjusted",
    x = expression(log[2]~"Fold Change"), 
    y = expression(-log[10]~italic(P))
  ) +
  coord_cartesian(clip = "off") + 
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

# Guardar
ggsave("Volcano_Final_1.png", plot = volcano_final, width = 15, height = 15, units = "cm", dpi = 300)

volcano_final
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-19-1.png" alt="" style="display: block; margin: auto;" />

``` r
## VOLCANO PLOT + GENE NAMES (from annotation table)


library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)


## 1. PREPARE DATA

plot_data <- comparison$ComparisonResult %>%
  mutate(
    log_p = -log10(pvalue),
    status = case_when(
      pvalue < 0.05 & log2FC > 1  ~ "Up-regulated",
      pvalue < 0.05 & log2FC < -1 ~ "Down-regulated",
      TRUE                        ~ "Not significant"
    )
  )

## 2. LOAD ANNOTATION TABLE (UniProt + Genes)


annotation_table <- read.xlsx("Uniprot_GeneName_Crack.xlsx")

annotation_table <- annotation_table %>%
  rename(
    Protein = Entry,
    Gene = Genes
  )

## 3. MERGE GENE NAMES

plot_data <- plot_data %>%
  left_join(annotation_table, by = "Protein")

# fallback (important)
plot_data$Gene[is.na(plot_data$Gene)] <- plot_data$Protein


## 4. TOP HITS

top_10_hits <- plot_data %>% 
  filter(!is.na(pvalue)) %>% 
  arrange(pvalue) %>% 
  head(10)

## Export top 10
write.xlsx(top_10_hits,
           "Top_10_Proteins_Volcano_unadjusted.xlsx",
           rowNames = FALSE)

## 5. VOLCANO PLOT


volcano_final <- ggplot(plot_data, aes(x = log2FC, y = log_p, color = status)) +
  
  geom_point(alpha = 0.5, size = 1.2) + 
  
  scale_color_manual(values = c(
    "Down-regulated" = "#0072B2", 
    "Up-regulated" = "#D55E00", 
    "Not significant" = "gray85"
  )) +
  
  scale_y_continuous(limits = c(0, 17)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) + 
  
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray40", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray40", size = 0.3) +
  
  geom_text_repel(
    data = top_10_hits,
    aes(label = Gene),
    size = 3.2, 
    color = "black", 
    fontface = "italic", 
    box.padding = 0.8,
    segment.size = 0.2,
    max.overlaps = Inf,
    force = 2
  ) +
  
  labs(
    title = "Proteomic Dysregulation",
    subtitle = "Model unadjusted",
    x = expression(log[2]~"Fold Change"), 
    y = expression(-log[10]~italic(P))
  ) +
  
  coord_cartesian(clip = "off") + 
  
  theme_bw() + 
  
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

## 6. SAVE FIGURE


ggsave(
  "Volcano_Final_GeneNames1.png",
  plot = volcano_final,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 300
)

volcano_final
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-20-1.png" alt="" style="display: block; margin: auto;" />

``` r
## 15. Heatmap 
library(pheatmap)
library(tidyr)
library(dplyr)
library(tibble)


# matrix
sig_proteins <- processed_final$ProteinLevelData %>%
  distinct(Protein) %>%
  pull(Protein)

heatmap_matrix <- processed_final$ProteinLevelData %>%
  filter(Protein %in% sig_proteins) %>%
  dplyr::select(Protein, SUBJECT, LogIntensities) %>%
  pivot_wider(names_from = SUBJECT, values_from = LogIntensities) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

# columns secuence
heatmap_matrix <- heatmap_matrix[, order(as.numeric(colnames(heatmap_matrix)))]

# cleaning
heatmap_matrix <- na.omit(heatmap_matrix)
heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, var) > 0, ]

label_map <- processed_final$ProteinLevelData %>%
  distinct(SUBJECT, GROUP) %>%
  mutate(SampleName = paste0("Intensity.", GROUP, "_", SUBJECT)) %>%
  filter(SUBJECT %in% colnames(heatmap_matrix)) %>%
  arrange(match(SUBJECT, colnames(heatmap_matrix)))

# rename
colnames(heatmap_matrix) <- label_map$SampleName

annotation_col <- data.frame(
  Condition = label_map$GROUP
)
rownames(annotation_col) <- label_map$SampleName

# 5. Colors
ann_colors <- list(
  Condition = c(
    AD = "#B57FA4",
    CN_PART = "#EBC4E1"
  )
)

# 6. Heatmap
pheatmap(
  heatmap_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#1C6CAA", "white", "firebrick3"))(100),
  border_color = "grey80",
  cellwidth = 15,
  cellheight = 10,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  main = ""
)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-21-1.png" alt="" style="display: block; margin: auto;" />

``` r
# 1. File generation
png("Heatmap_Zscore_Final.png", width = 2200, height = 2200, res = 300)

# 2. Legend
pheatmap(
  heatmap_matrix, 
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#1C6CAA", "white", "firebrick3"))(100),
  border_color = "grey80",
  cellwidth = 15, 
  cellheight = 10,
  show_colnames = TRUE, 
  fontsize_row = 7, 
  fontsize_col = 8,
  name = "Z-score", 
  main = "" 
)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-22-1.png" alt="" style="display: block; margin: auto;" />

``` r
# 3. Close and save the file
dev.off()
```

    ## png 
    ##   3

## 15. Differential testing adjusting AGE , PMI , SEX

MSstats requires covariates in ProteinLevelData.

### 15.1 Merge metadata:

``` r
metadata <- read.table(
  "experimental_annotationfull.txt",
  header = TRUE,
  sep = "\t"
)
```

### 15.2 Data Transformation

Limma expects an “Expression Matrix” where rows are proteins and columns
are samples.

``` r
library(dplyr)
library(stringr)
library(tidyr)
library(limma)

# --- Step A: Clean Metadata ---
metadata_clean <- metadata %>%
  mutate(
    # Extract the number from "Intensity.AD3", "Intensity.PT9", etc.
    join_id = str_extract(run, "\\d+"),
    # Fix commas in PMI if they exist (common in some regions) and force numeric
    PMI_num = as.numeric(str_replace(as.character(PMI), ",", ".")),
    AGE_num = as.numeric(as.character(AGE))
  )

# --- Step B: Clean Protein Data ---
protein_data_clean <- processed_final$ProteinLevelData %>%
  mutate(
    # Ensure originalRUN is a character to match join_id
    join_id = as.character(originalRUN)
  )

# --- Step C: The Join ---
# This merges your protein levels with Sex, AGE, and PMI
library(dplyr)
combined_data <- protein_data_clean %>%
  left_join(metadata_clean %>% dplyr::select(join_id, Sex, AGE_num, PMI_num), 
            by = "join_id")
```

``` r
library(limma)
library(tidyr)

# 1.Convert to wide format (protein matrix).
data_wide <- combined_data %>%
  # We select the key columns
  dplyr::select(Protein, originalRUN, LogIntensities) %>%
  # We convert samples into columns
  pivot_wider(names_from = originalRUN, values_from = LogIntensities) %>%
  as.data.frame()

# Set row names to proteins and clean up
rownames(data_wide) <- data_wide$Protein
data_wide$Protein <- NULL

# 2. Create the aligned metadata object
# We filter one row per sample with its covariates
design_meta <- combined_data %>%
  distinct(originalRUN, GROUP, Sex, AGE_num, PMI_num) %>%
  # IMPORTANT: The row order must be the SAME as the column order in data_wide
  arrange(match(originalRUN, colnames(data_wide)))
```

``` r
# Create the design matrix. 0 + GROUP allows direct comparison of AD vs CN_PART
design <- model.matrix(~ 0 + GROUP + Sex + AGE_num + PMI_num, data = design_meta)

# lean column names to avoid errors in contrasts
colnames(design) <- make.names(colnames(design))
```

``` r
# 1. Fit the model to the data
fit <- lmFit(data_wide, design)

# 2. Define the contrast of interest (AD vs Control)
# Make sure to use the exact names that appear in colnames(design)
cont_matrix <- makeContrasts(
  AD_vs_Control = GROUPAD - GROUPCN_PART, 
  levels = design
)

# 3. Apply the contrast and compute moderated statistics (eBayes)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# 4. GGenerate the final results table
# The BH (Benjamini–Hochberg) method is used to control the FD
top_table <- topTable(fit2, coef = "AD_vs_Control", number = Inf, adjust.method = "BH")

# Show the most significant proteins adjusted for covariates
head(top_table)
```

    ##            logFC  AveExpr         t      P.Value adj.P.Val           B
    ## Q86X76  2.459432 22.61303  6.876893 0.0006053618 0.3472857  0.31828283
    ## P55011  2.321928 22.09253  6.492415 0.0008113258 0.3472857  0.05962508
    ## Q07954  2.000000 22.70554  5.592262 0.0017056669 0.4712304 -0.62078561
    ## Q15173 -6.000000 24.08334 -6.650104 0.0007183205 0.3472857 -1.15364248
    ## Q16762 -5.727920 24.01540 -6.348544 0.0009085304 0.3472857 -1.24783430
    ## O60506 -4.807355 23.78565 -5.328235 0.0021573660 0.4712304 -1.63762510

``` r
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Security filter
results_cleaned <- top_table %>%
  filter(abs(logFC) < 20) %>% 
  mutate(Protein = rownames(.))

# prepare the data
results_adj_plot <- results_cleaned %>%
  mutate(Status = case_when(
    P.Value < 0.05 & logFC > 0.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# Top 10
top_labels <- results_adj_plot %>%
  filter(Status != "Not Significant") %>%
  arrange(P.Value) %>%
  head(10)

# volcano
volcano2 <- ggplot(results_adj_plot, aes(x = logFC, y = -log10(P.Value), color = Status)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(
    "Down-regulated" = "#0072B2", 
    "Up-regulated" = "#D55E00", 
    "Not Significant" = "gray85"
  )) +
  geom_text_repel(
    data = top_labels, aes(label = Protein),
    size = 3, 
    color = "black", 
    fontface = "italic",
    box.padding = 0.5, 
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  
 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray40", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray40", size = 0.3) +
  
 
  scale_x_continuous(limits = c(-6, 6)) + 
  
  theme_bw() +
  labs(
    title = "Proteomic Dysregulation",
    subtitle = "Model adjusted for Age, Sex, and PMI",
    x = expression(log[2]~"Fold Change"), 
    y = expression(-log[10]~italic(P))
  ) +
  theme(
    aspect.ratio = 1,                 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "top",          
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10, color = "black")
  )

# save
ggsave("volcano_final_estilo_limpio.png", volcano2, width = 15, height = 15, units = "cm", dpi = 300)

volcano2
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-28-1.png" alt="" style="display: block; margin: auto;" />

``` r
## VOLCANO PLOT 2 + GENE NAMES (from annotation table)

library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)


## 1. SAFETY FILTERING

results_cleaned <- top_table %>%
  filter(abs(logFC) < 20) %>% 
  mutate(Protein = rownames(.))


## 2. DATA PREPARATION

results_adj_plot <- results_cleaned %>%
  mutate(Status = case_when(
    P.Value < 0.05 & logFC > 0.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))


## 3. LOAD + CLEAN ANNOTATION TABLE


annotation_table <- read.xlsx("Uniprot_GeneName_Crack.xlsx")

annotation_table <- annotation_table %>%
  rename(
    Protein = Entry,
    Gene = Genes
  )


## 4. MERGE GENE NAMES INTO VOLCANO DATA


results_adj_plot <- results_adj_plot %>%
  left_join(annotation_table, by = "Protein")

# fallback (important for missing annotations)
results_adj_plot$Gene[is.na(results_adj_plot$Gene)] <- results_adj_plot$Protein


## 5. TOP LABELS

top_labels <- results_adj_plot %>%
  filter(Status != "Not Significant") %>%
  arrange(P.Value) %>%
  head(10)

## 6. VOLCANO PLOT

volcano2 <- ggplot(results_adj_plot, aes(x = logFC, y = -log10(P.Value), color = Status)) +
  
  geom_point(alpha = 0.5, size = 1.2) +
  
  scale_color_manual(values = c(
    "Down-regulated" = "#0072B2", 
    "Up-regulated" = "#D55E00", 
    "Not Significant" = "gray85"
  )) +
  
  geom_text_repel(
    data = top_labels,
    aes(label = Gene),
    size = 3, 
    color = "black", 
    fontface = "italic",
    box.padding = 0.5, 
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray40", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray40", size = 0.3) +
  
  scale_x_continuous(limits = c(-6, 6)) +
  
  theme_bw() +
  
  labs(
    title = "Proteomic Dysregulation",
    subtitle = "Model adjusted for Age, Sex, and PMI",
    x = expression(log[2]~"Fold Change"), 
    y = expression(-log[10]~italic(P))
  ) +
  
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10, color = "black")
  )

## 7. SAVE FIGURE


ggsave(
  "volcano_final_clean_geneNames.png",
  volcano2,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 300
)

volcano2
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-29-1.png" alt="" style="display: block; margin: auto;" />

``` r
# COUNTS OF REGULATED PROTEINS

# Option 1: Clean summary table using dplyr
status_summary <- results_adj_plot %>%
  count(Status, name = "Count") %>%
  mutate(Percentage = round((Count / sum(Count)) * 100, 1))

print("--- Volcano Plot Regulation Summary ---")
```

    ## [1] "--- Volcano Plot Regulation Summary ---"

``` r
print(status_summary)
```

    ##            Status Count Percentage
    ## 1  Down-regulated    16        1.0
    ## 2 Not Significant  1459       95.7
    ## 3    Up-regulated    49        3.2

``` r
# Option 2: Quick base R table view
table(results_adj_plot$Status)
```

    ## 
    ##  Down-regulated Not Significant    Up-regulated 
    ##              16            1459              49

``` r
library(circlize) 
library(ggplot2)
library(ComplexHeatmap)
```

``` r
processed_final$ProteinLevelData %>%
  distinct(SUBJECT, originalRUN) %>%
  arrange(as.numeric(SUBJECT))
```

    ##    SUBJECT originalRUN
    ## 1        1           9
    ## 2        2           1
    ## 3        3           2
    ## 4        4          10
    ## 5        5          11
    ## 6        6           3
    ## 7        7           4
    ## 8        8          12
    ## 9        9          13
    ## 10      10          14
    ## 11      11           5
    ## 12      12          15
    ## 13      13          16
    ## 14      14           6
    ## 15      15           7
    ## 16      16          17
    ## 17      17           8
    ## 18      18          18

``` r
subject_map <- processed_final$ProteinLevelData %>%
  distinct(SUBJECT, originalRUN) %>%
  mutate(
    SUBJECT = as.character(SUBJECT),
    originalRUN = as.character(originalRUN)
  )

design_meta2 <- design_meta %>%
  mutate(originalRUN = as.character(originalRUN)) %>%
  left_join(subject_map, by = "originalRUN")
```

``` r
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(tibble)


# 1. Building matrix from processed_final
sig_proteins <- top_table %>%
  tibble::rownames_to_column("Protein") %>%
  filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
  pull(Protein)

# 2. matrix
heatmap_mat <- processed_final$ProteinLevelData %>%
  filter(Protein %in% sig_proteins) %>%
  dplyr::select(Protein, SUBJECT, LogIntensities) %>%
  tidyr::pivot_wider(names_from = SUBJECT, values_from = LogIntensities) %>%
  tibble::column_to_rownames("Protein") %>%
  as.matrix()

# colums by SUBJECT
heatmap_mat <- heatmap_mat[, order(as.numeric(colnames(heatmap_mat)))]

# cleaning
heatmap_mat <- na.omit(heatmap_mat)
heatmap_mat <- heatmap_mat[apply(heatmap_mat, 1, var) > 0, ]

# 2. label_map 

label_map <- processed_final$ProteinLevelData %>%
  distinct(SUBJECT, GROUP) %>%
  mutate(SUBJECT = as.character(SUBJECT)) %>%
  filter(SUBJECT %in% colnames(heatmap_mat)) %>%
  arrange(match(SUBJECT, colnames(heatmap_mat))) %>%
  mutate(SampleName = paste0("Intensity.", GROUP, "_", SUBJECT))


colnames(heatmap_mat) <- label_map$SampleName


# 3. Z-score


heatmap_mat_scaled <- t(apply(heatmap_mat, 1, scale))
colnames(heatmap_mat_scaled) <- colnames(heatmap_mat)
heatmap_mat_scaled <- heatmap_mat_scaled[is.finite(rowSums(heatmap_mat_scaled)), ]


# 4. Metadata -->design_meta2

meta_df <- design_meta2 %>%
  mutate(SUBJECT = as.character(SUBJECT)) %>%
  filter(SUBJECT %in% label_map$SUBJECT) %>%
  distinct(SUBJECT, .keep_all = TRUE) %>%
  dplyr::select(SUBJECT, Sex, AGE_num, PMI_num)

# annotation 
anno_df <- label_map %>%
  dplyr::select(SUBJECT, SampleName, GROUP) %>%
  left_join(meta_df, by = "SUBJECT")

# 
anno_df <- anno_df[match(colnames(heatmap_mat_scaled), anno_df$SampleName), ]

# rownames
rownames(anno_df) <- anno_df$SampleName

# final columns
anno_df <- anno_df %>%
  dplyr::select(GROUP, Sex, AGE_num, PMI_num)

# 5. Colors

anno_cols <- list(
  GROUP = c(
    AD = "#B57FA4",
    CN_PART = "#EBC4E1"
  ),
  Sex = c(
    Female = "#FFFF37",
    Male = "#A2FF00"
  ),
  AGE_num = colorRamp2(
    c(min(anno_df$AGE_num, na.rm = TRUE),
      max(anno_df$AGE_num, na.rm = TRUE)),
    c("white", "darkgreen")
  ),
  PMI_num = colorRamp2(
    c(min(anno_df$PMI_num, na.rm = TRUE),
      max(anno_df$PMI_num, na.rm = TRUE)),
    c("white", "darkorange")
  )
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#1C6CAA", "white", "firebrick")
)


# 6. Anotación superior

top_anno <- HeatmapAnnotation(
  df = anno_df,
  col = anno_cols,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 7)
)


# 7. Heatmap

celda_size <- unit(5, "mm")

ht <- Heatmap(
  heatmap_mat_scaled,
  name = "Z-score",
  col = col_fun,
  top_annotation = top_anno,
  width = ncol(heatmap_mat_scaled) * celda_size,
  height = nrow(heatmap_mat_scaled) * celda_size,
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_km = 2,
  row_km = 2,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7),
  column_title = "Análisis Compacto de Proteínas",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)


# 8. Save

png("Heatmap_Final.png", width = 2400, height = 2400, res = 300)
draw(ht, merge_legend = TRUE)
dev.off()
```

    ## png 
    ##   2

``` r
#
ht
```

### Top hits

``` r
library(openxlsx)
library(dplyr)


## 1. BASE TABLE + STATUS
results_for_table <- top_table %>%
  mutate(Protein = rownames(.)) %>%
  mutate(Status = case_when(
    P.Value < 0.05 & logFC > 0.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

## 2. ANNOTATION (UniProt -> Gene)


annotation_table <- annotation_table %>%
  rename(
    Protein = Protein,   
    Gene = Gene
  )

## 3. MERGE GENE NAMES

results_annotated <- results_for_table %>%
  left_join(annotation_table, by = "Protein")

# fallback para genes faltantes
results_annotated$Gene[is.na(results_annotated$Gene)] <- results_annotated$Protein


## 4. TOTAL SIGNIFICANT PROTEINS 

all_significant <- results_annotated %>%
  filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
  mutate(Status = if_else(logFC > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(P.Value)

n_total_significant <- nrow(all_significant)


## 5. TOP 10 PROTEINS (TABLE 3)


top_10 <- all_significant %>%
  arrange(P.Value) %>%
  head(10)

## 6. FINAL EXPORT TABLE (WITH GENE NAMES)

table3_final <- top_10 %>%
  dplyr::select(
    Protein,
    Gene,
    logFC,
    P.Value,
    adj.P.Val,
    AveExpr,
    Status
  )

file_table3 <- paste0("Table3_Top10_with_GeneNames_", Sys.Date(), ".xlsx")
file_all <- paste0("All_Significant_Proteins_", Sys.Date(), ".xlsx")

write.xlsx(table3_final, file = file_table3, rowNames = FALSE)
write.xlsx(all_significant, file = file_all, rowNames = FALSE)


## 7. OUTPUT MESSAGE (FOR PAPER TEXT)

message("Table 3 saved: ", file_table3)
message("Full significant list saved: ", file_all)
message("Total significant proteins (for manuscript text): ", n_total_significant)
```

``` r
library(openxlsx)
library(dplyr)


## 1. BASE TABLE + STATUS

results_for_table <- top_table %>%
  mutate(Protein = rownames(.)) %>%
  mutate(Status = case_when(
    P.Value < 0.05 & logFC > 0.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))


## 2. ANNOTATION TABLE
## (must contain: Protein, Gene, Protein_name)

annotation_table <- annotation_table %>%
  rename(
    Protein = Protein,
    Gene = Gene,
    Protein_name = Protein.names   
  )


## 3. MERGE

results_annotated <- results_for_table %>%
  left_join(annotation_table, by = "Protein")

# fallbacks
results_annotated$Gene[is.na(results_annotated$Gene)] <- results_annotated$Protein
results_annotated$Protein_name[is.na(results_annotated$Protein_name)] <- "Unknown"


## 4. TOTAL SIGNIFICANT PROTEINS 

all_significant <- results_annotated %>%
  filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
  mutate(Status = if_else(logFC > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(P.Value)

n_total_significant <- nrow(all_significant)


## 5. TOP 10 (TABLE 3)


top_10 <- all_significant %>%
  arrange(P.Value) %>%
  head(10)


## 6. FINAL EXPORT TABLE


table3_final <- top_10 %>%
  dplyr::select(
    Protein,
    Gene,
    Protein_name,
    logFC,
    P.Value,
    adj.P.Val,
    AveExpr,
    Status
  )

write.xlsx(
  table3_final,
  file = paste0("Table3_Top10_FULL_", Sys.Date(), ".xlsx"),
  rowNames = FALSE
)

write.xlsx(
  all_significant,
  file = paste0("All_Significant_Proteins_+", Sys.Date(), ".xlsx"),
  rowNames = FALSE
)

## =========================================================
## 7. OUTPUT
## =========================================================

message("Top 10 table saved with Protein name + Gene")
message("Total significant proteins: ", n_total_significant)
```

``` r
library(openxlsx)
library(dplyr)

# 1. Create the Protein column and define the status as in the plot.
results_for_excel <- top_table %>%
  mutate(Protein = rownames(.)) %>%
  mutate(Status = case_when(
    P.Value < 0.05 & logFC > 0.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# 2. Filter to obtain exactly the ‘Top Hits’ from the plot.

top_10_hits <- results_for_excel %>%
  filter(Status != "Not Significant") %>%
  arrange(P.Value) %>%
  head(10)

# 3. Excel
file_name <- paste0("Top10_Labels_Volcano_", Sys.Date(), ".xlsx")

# top_10_hits
write.xlsx(top_10_hits, file = file_name, rowNames = FALSE)

message("¡Solucionado! Se ha guardado el Top 10 del gráfico en: ", file_name)
```

``` r
library(openxlsx)
library(dplyr)

# 1. We identify all colored proteins (Up and Down)
all_significant_colored <- top_table %>%
  mutate(Protein = rownames(.)) %>%
  filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
  mutate(Status = if_else(logFC > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(P.Value)

# 2.  Excel
file_name_all <- paste0("All_Colored_Proteins_", Sys.Date(), ".xlsx")

# 
write.xlsx(all_significant_colored, file = file_name_all, rowNames = FALSE)

# 
message("¡Éxito! Se han guardado ", nrow(all_significant_colored), 
        " proteínas coloreadas en el archivo: ", file_name_all)
```

## 16. Correlation with pathology

``` r
library(readxl)
library(limma)
library(dplyr)
library(tidyr)
library(reshape2)
library(openxlsx)
library(stringr)

# -Data
metadata_patho <- read_excel("experimental_annotationfull.xlsx")

muestras_ordenadas <- c(
  "Intensity.AD3", "Intensity.AD6", "Intensity.AD11", "Intensity.AD17", 
  "Intensity.AD2", "Intensity.AD7", "Intensity.PT9", "Intensity.AD14", 
  "Intensity.AD15", "Intensity.PT1", "Intensity.CN4", "Intensity.PT5", 
  "Intensity.CN8", "Intensity.CN10", "Intensity.PT12", "Intensity.CN13", 
  "Intensity.PT16", "Intensity.CN18"
)

design_patho <- metadata_patho %>%
  mutate(
    join_id = as.character(match(run, muestras_ordenadas)),
    Thal = as.numeric(Thal),
    Braak.Tau = as.numeric(Braak.Tau),
    AGE_num = as.numeric(as.character(AGE)),
    PMI_num = as.numeric(as.character(PMI)),
    Sex = as.factor(Sex)
  ) %>%
  filter(!is.na(Thal), !is.na(Braak.Tau), !is.na(join_id))

# 2: Preparation of the Protein Matrix

#IMPORTANT: We use originalRUN because these are the IDs (1–18) that MSstats stores internally
data_matrix_patho <- acast(processed_final$ProteinLevelData, 
                           Protein ~ originalRUN, 
                           value.var = "LogIntensities")
data_matrix_patho <- as.matrix(data_matrix_patho)

# 3: : Final Synchronization -Merge by sample numbers (1, 2, 3, …)
common_samples <- intersect(design_patho$join_id, colnames(data_matrix_patho))

if(length(common_samples) == 0) {
  stop("ERROR: Sigue habiendo 0 coincidencias. Revisa si 'originalRUN' son números o nombres en tu objeto processed_final.")
}

design_patho_final <- design_patho %>% 
  filter(join_id %in% common_samples) %>%
  arrange(match(join_id, common_samples))

data_matrix_final <- data_matrix_patho[, design_patho_final$join_id]

# 4: Models

# 1. MODEL BRAAK (Tau)
design_mat_braak <- model.matrix(~ Braak.Tau + AGE_num + Sex + PMI_num, data = design_patho_final)
fit_braak <- lmFit(data_matrix_final, design_mat_braak)
fit_braak <- eBayes(fit_braak)
res_braak <- topTable(fit_braak, coef = 2, number = Inf)

# 2. MODEL THAL (Amiloide)
design_mat_thal <- model.matrix(~ Thal + AGE_num + Sex + PMI_num, data = design_patho_final)
fit_thal <- lmFit(data_matrix_final, design_mat_thal)
fit_thal <- eBayes(fit_thal)
res_thal <- topTable(fit_thal, coef = 2, number = Inf)

# 3. MODEL INTERACCTION (Braak x Thal)
design_mat_int <- model.matrix(~ Thal * Braak.Tau + AGE_num + Sex + PMI_num, data = design_patho_final)
fit_int <- lmFit(data_matrix_final, design_mat_int)
fit_int <- eBayes(fit_int)
res_int <- topTable(fit_int, coef = ncol(design_mat_int), number = Inf)

# 5: SAVE
wb <- createWorkbook()
addWorksheet(wb, "Braak_Tau")
writeData(wb, "Braak_Tau", res_braak, rowNames = TRUE)
addWorksheet(wb, "Thal_Amyloid")
writeData(wb, "Thal_Amyloid", res_thal, rowNames = TRUE)
addWorksheet(wb, "Interaccion_Thal_Braak")
writeData(wb, "Interaccion_Thal_Braak", res_int, rowNames = TRUE)

saveWorkbook(wb, "Resultados_Correlacion_Patologia_Final.xlsx", overwrite = TRUE)

message("Análisis completado exitosamente con ", length(common_samples), " muestras.")
```

## S1: Enrichment validation

``` r
# 1. mARKERS
lista_marcadores <- data.frame(
  Protein = c("P08247", "P78352", "P17600", "P60880", "P14136", "P09543", "P16401", "P00338"),
  Symbol = c("SYP", "DLG4", "SYN1", "SNAP25", "GFAP", "CNP", "HIST1H1B", "LDHA"),
  Type = c("Synaptic", "Synaptic", "Synaptic", "Synaptic", "Non-Synaptic", "Non-Synaptic", "Non-Synaptic", "Non-Synaptic"),
  Category = c("Presynaptic", "Postsynaptic", "Presynaptic", "Membrane", "Glia", "Myelin", "Nuclear", "Cytosolic")
)

# 2. PREPARE DATAFRAME
df_nature <- as.data.frame(data_matrix_final) %>%
  filter(rownames(.) %in% lista_marcadores$Protein) %>%
  mutate(Protein = rownames(.)) %>%
  left_join(lista_marcadores, by = "Protein") %>%
  pivot_longer(cols = -c(Protein, Symbol, Type, Category), names_to = "Sample", values_to = "Intensity")

# 3. PLOT
fig_s1 <- ggplot(df_nature, aes(x = reorder(Symbol, -Intensity), y = Intensity, fill = Type)) +
  geom_boxplot(outlier.shape = NA, color = "#222222", size = 0.4, width = 0.6) +
  geom_jitter(color = "#444444", width = 0.15, size = 0.8, alpha = 0.4) +
  facet_grid(. ~ Type, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("Synaptic" = "#87CEFA", "Non-Synaptic" = "#CD919E")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    legend.position = "none"
  ) +
  labs(
    title = "Validation of Synaptosomal Enrichment",
    x = "Protein Marker",
    y = "Log2 Intensity"
  )

print(fig_s1)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-39-1.png" alt="" style="display: block; margin: auto;" />

``` r
ggsave("Figure_S1_Validation.pdf", fig_s1, width = 7, height = 5)
```

## S2. Global protein intensity distributions for quality control.

``` r
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(stringr)


mis_muestras <- colnames(data_matrix_final)
```

``` r
prep_data_qc <- function(obj, name) {
  temp <- as.data.frame(obj)
  
  # 1. CASE: MSstats format (long dataframe)
  if("originalRUN" %in% colnames(temp)) {
    # Detect intensity and protein columns
    int_col <- intersect(colnames(temp), c("LogIntensities", "Intensity"))[1]
    prot_col <- intersect(colnames(temp), c("Protein", "PROTEIN", "protein"))[1]
    
    df_plot <- temp %>%
      # Remove as.numeric to avoid NA warnings
      mutate(Sample = as.character(originalRUN)) %>% 
      select(Protein = !!sym(prot_col), Sample, Intensity = !!sym(int_col))
  } 
  
  # 2. Matrix format (wide)
  else {
    # Filter only numeric columns (ignore the ‘PROTEIN’ column if it is present)
    mat_pure <- temp[, sapply(temp, is.numeric), drop = FALSE]
    
    # Transform matrix to long format
    df_melted <- reshape2::melt(as.matrix(mat_pure))
    
    df_plot <- df_melted
    # Enforce standard column names
    colnames(df_plot) <- c("Protein", "Sample", "Intensity")
    df_plot$Sample <- as.character(df_plot$Sample)
  }
  
  # 3. CLEANING AND LOG2 TRANSFORMATION
  df_plot <- df_plot %>% filter(!is.na(Intensity) & is.finite(Intensity))
  
  # If there are no rows after cleaning, the object was empty.
  if(nrow(df_plot) == 0) stop(paste("El objeto", name, "no tiene datos válidos."))

  # Apply log2 only if the data is on a linear scale (values greater than 1000).
  if(max(df_plot$Intensity, na.rm=TRUE) > 1000) {
    df_plot$Intensity <- log2(df_plot$Intensity)
  }
  
  df_plot$Step <- name
  return(df_plot)
}
```

``` r
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
df1 <- prep_data_qc(feature_norm, "a | Raw Intensities")
df2 <- prep_data_qc(processed_final$ProteinLevelData, "b | MSstats Normalization")
df3 <- prep_data_qc(data_matrix_final, "c | Post 60-Protein Refinement")
```

``` r
# 

theme_qc <- theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text.x = element_blank(), 
  axis.ticks.x = element_blank(),
  plot.title = element_text(size = 10)
)

# Colors
p1 <- ggplot(df1, aes(x=Sample, y=Intensity)) + 
      geom_boxplot(fill="#C6E2FF", alpha=0.8, outlier.size = 0.2) + 
      theme_qc + labs(title=unique(df1$Step), y="log2 Intensity")

p2 <- ggplot(df2, aes(x=Sample, y=Intensity)) + 
      geom_boxplot(fill="#7EC0EE", alpha=0.8, outlier.size = 0.2) + 
      theme_qc + labs(title=unique(df2$Step), y="log2 Intensity")

p3 <- ggplot(df3, aes(x=Sample, y=Intensity)) + 
      geom_boxplot(fill="#4682B4", alpha=0.8, outlier.size = 0.2) + 
      theme_bw() + theme(panel.grid = element_blank(), 
                         axis.text.x = element_text(angle=90, size=8, vjust=0.5, hjust=1)) +
      labs(title=unique(df3$Step), x="Samples (n=18)", y="log2 Intensity")

#

fig_s2_qc <- p1 / p2 / p3 + 
  plot_annotation(title = "Figure S2: Proteomics Data Normalization Workflow",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

# 
print(fig_s2_qc)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-43-1.png" alt="" style="display: block; margin: auto;" />

``` r
# 
ggsave("Figure_S2_QC_Distributions.png", fig_s2_qc, width = 8.5, height = 11, dpi = 300)
```

## 17. Enrichment analysis

### 17.0 Analysis of the filtered proteins

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot"
))
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
```

Uni ProtIDs?

``` r
head(feature_data_filtered$ProteinName)
```

    ## [1] "Q8N3E9" "P23471" "Q5TF21" "Q6H8Q1" "Q14C86" "Q9P1U1"

``` r
proteins <- unique(feature_data_filtered$ProteinName)
```

Convert Convertir UniProt → ENTREZ

``` r
protein_mapping <- bitr(
  proteins,
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)
```

List for enrichment analysis

``` r
entrez_ids <- unique(protein_mapping$ENTREZID)
```

#### BiologicalProcesses

``` r
ego_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
```

#### Cellular Components

``` r
ego_cc <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
```

``` r
ego_cc_df <- as.data.frame(ego_cc)

top10_CC <- ego_cc_df %>%
  arrange(p.adjust) %>%
  head(10)
```

``` r
# SYNAPTIC PROTEOME OVER-REPRESENTATION ANALYSIS
# Using ALL detected proteins after QC filtering
# Goal:
# Demonstrate enrichment of synaptic-related pathways

# -----------------------------
# 1. Libraries
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

library(dplyr)
library(ggplot2)
library(stringr)

# -----------------------------
# 2. Extract protein list
# -----------------------------
proteins <- unique(feature_data_filtered$ProteinName)

# Clean UniProt IDs
proteins_clean <- gsub("sp\\||tr\\||\\|.*", "", proteins)

# -----------------------------
# 3. Convert UniProt -> ENTREZ
# -----------------------------
protein_mapping <- bitr(
  proteins_clean,
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# ENTREZ IDs for enrichment
entrez_ids <- unique(protein_mapping$ENTREZID)

# -----------------------------
# 4. GO Cellular Component enrichment
# -----------------------------
ego_cc <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# -----------------------------
# 5. Select TOP10 enriched terms
# -----------------------------
top10_CC <- as.data.frame(ego_cc) %>%
  arrange(p.adjust) %>%
  dplyr::slice(1:10) %>%

  # Calculate Fold Enrichment
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio)),

    BgRatio_num = as.numeric(sub("/.*", "", BgRatio)) /
      as.numeric(sub(".*/", "", BgRatio)),

    FoldEnrichment = GeneRatio_num / BgRatio_num
  )

# -----------------------------
# 6. Nature-inspired palette
# -----------------------------
nature_palette <- c(
  "#08306B",
  "#08519C",
  "#2171B5",
  "#4292C6",
  "#6BAED6"
)

# -----------------------------
# 7. Generate enrichment plot
# -----------------------------
p_synapse <- ggplot(
  top10_CC,
  aes(
    x = reorder(Description, FoldEnrichment),
    y = FoldEnrichment,
    fill = p.adjust
  )
) +

  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.2
  ) +

  # Add protein counts
  geom_text(
    aes(label = Count),
    hjust = -0.2,
    size = 3
  ) +

  coord_flip() +

  scale_fill_gradientn(
    colours = nature_palette,
    trans = "reverse",
    name = "FDR",
    labels = function(x)
      format(x, scientific = TRUE, digits = 1)
  ) +

  labs(
    title = "Synaptic Proteome Composition",
    subtitle = "GO Cellular Component over-representation analysis",
    x = NULL,
    y = "Fold Enrichment"
  ) +

  theme_classic(base_size = 11) +

  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

# -----------------------------
# 8. Save figure
# -----------------------------
ggsave(
  "Synaptic_Proteome_Composition_GOCC.png",
  p_synapse,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300
)

# -----------------------------
# 9. Display plot
# -----------------------------
p_synapse
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-52-1.png" alt="" style="display: block; margin: auto;" />

``` r
print("Synaptic enrichment figure successfully generated.")
```

    ## [1] "Synaptic enrichment figure successfully generated."

Reactome enrichment

``` r
reactome_res <- enrichPathway(
  gene         = entrez_ids,
  organism     = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable     = TRUE
)
```

Results (top pathways)

``` r
head(as.data.frame(ego_cc))
```

    ##                    ID                                  Description GeneRatio
    ## GO:0005925 GO:0005925                               focal adhesion  183/2265
    ## GO:0030055 GO:0030055                      cell-substrate junction  183/2265
    ## GO:0098978 GO:0098978                        glutamatergic synapse  172/2265
    ## GO:0098800 GO:0098800 inner mitochondrial membrane protein complex   98/2265
    ## GO:0098798 GO:0098798     mitochondrial protein-containing complex  134/2265
    ## GO:0008021 GO:0008021                             synaptic vesicle  110/2265
    ##              BgRatio       pvalue     p.adjust       qvalue
    ## GO:0005925 421/19886 4.398006e-64 2.858704e-61 1.226812e-61
    ## GO:0030055 431/19886 3.887198e-62 1.263339e-59 5.421618e-60
    ## GO:0098978 407/19886 4.955312e-58 1.073651e-55 4.607571e-56
    ## GO:0098800 158/19886 1.254012e-52 2.037770e-50 8.745086e-51
    ## GO:0098798 300/19886 1.150710e-48 1.495923e-46 6.419752e-47
    ## GO:0008021 213/19886 7.313133e-48 7.922561e-46 3.399965e-46
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            geneID
    ## GO:0005925 FLOT2/SCARB2/AHNAK/CYFIP1/PTK2B/TNS3/ANXA6/PPIB/ANXA5/PLEC/MAPK3/EZR/MSN/ACTN4/CLASP2/TLN1/CAPN1/PRKAR2A/MAP4K4/DST/LRP1/TLN2/MYH9/YWHAB/PDIA3/ARHGEF2/CASK/ARPC2/HACD3/CAPN5/EPB41L2/GAK/RPS3/CLASP1/ITGAV/WASF1/RAB21/HSPA1A/HSPA1B/TWF1/MAP2K2/CORO1C/HSPA5/CDH2/RPL12/CTTN/BSG/HSP90B1/GDI2/AKAP12/RDX/NCKAP1/FLNA/SNTB1/PTPRA/ACTN1/HYOU1/CTNNB1/EHD3/CLTC/RPS14/TPM4/P4HB/ACTC1/ANXA1/RPLP0/PRUNE1/RPL7/PI4KA/DNM2/ADD1/CAPN2/GNA13/ACTN2/VIM/CNN3/SH3KBP1/FLII/RPL6/TNC/RPS5/VCL/CAT/RAC1/MCAM/RPL4/NCSTN/RPS19/CDH13/CD44/PFN1/CORO1B/PDCD6IP/MAP2K1/ACTR2/PDPK1/ARHGAP26/RPL10A/YWHAQ/GJA1/GIT1/RPS7/HSPA9/SRP68/CHP1/FERMT2/HSPA8/PIP5K1C/YWHAG/CORO2B/LASP1/RAB10/CAP1/CD99L2/ARF6/RALA/HNRNPK/LAP3/PAK1/ITGB1/CPNE3/ARHGEF7/RPL3/SDCBP/PTK2/RPS4X/RPS8/ACTB/EVL/RPS2/FLOT1/GSN/RPL18/MAPRE1/DCTN4/MAPK1/L1CAM/ATAT1/ARPC3/YES1/NFASC/CFL1/RPLP2/RPS16/SORBS1/RHOB/PPP1R12A/RPS11/SRC/NPM1/YWHAZ/ACTR3/ALCAM/ARPC5L/ARL2/ITGB8/YWHAE/RPS9/RPL27/CALR/LAMTOR3/RHOG/PPIA/SLC9A1/RHOA/RPS13/HSPB1/CD99/B2M/KRAS/MAPRE2/ARPC5/PABPC1/MARCKS/THY1/PPP1CB/CD9/CSRP1/CDC42/FHL1/PCBP2/CD81/ATP6V0C
    ## GO:0030055 FLOT2/SCARB2/AHNAK/CYFIP1/PTK2B/TNS3/ANXA6/PPIB/ANXA5/PLEC/MAPK3/EZR/MSN/ACTN4/CLASP2/TLN1/CAPN1/PRKAR2A/MAP4K4/DST/LRP1/TLN2/MYH9/YWHAB/PDIA3/ARHGEF2/CASK/ARPC2/HACD3/CAPN5/EPB41L2/GAK/RPS3/CLASP1/ITGAV/WASF1/RAB21/HSPA1A/HSPA1B/TWF1/MAP2K2/CORO1C/HSPA5/CDH2/RPL12/CTTN/BSG/HSP90B1/GDI2/AKAP12/RDX/NCKAP1/FLNA/SNTB1/PTPRA/ACTN1/HYOU1/CTNNB1/EHD3/CLTC/RPS14/TPM4/P4HB/ACTC1/ANXA1/RPLP0/PRUNE1/RPL7/PI4KA/DNM2/ADD1/CAPN2/GNA13/ACTN2/VIM/CNN3/SH3KBP1/FLII/RPL6/TNC/RPS5/VCL/CAT/RAC1/MCAM/RPL4/NCSTN/RPS19/CDH13/CD44/PFN1/CORO1B/PDCD6IP/MAP2K1/ACTR2/PDPK1/ARHGAP26/RPL10A/YWHAQ/GJA1/GIT1/RPS7/HSPA9/SRP68/CHP1/FERMT2/HSPA8/PIP5K1C/YWHAG/CORO2B/LASP1/RAB10/CAP1/CD99L2/ARF6/RALA/HNRNPK/LAP3/PAK1/ITGB1/CPNE3/ARHGEF7/RPL3/SDCBP/PTK2/RPS4X/RPS8/ACTB/EVL/RPS2/FLOT1/GSN/RPL18/MAPRE1/DCTN4/MAPK1/L1CAM/ATAT1/ARPC3/YES1/NFASC/CFL1/RPLP2/RPS16/SORBS1/RHOB/PPP1R12A/RPS11/SRC/NPM1/YWHAZ/ACTR3/ALCAM/ARPC5L/ARL2/ITGB8/YWHAE/RPS9/RPL27/CALR/LAMTOR3/RHOG/PPIA/SLC9A1/RHOA/RPS13/HSPB1/CD99/B2M/KRAS/MAPRE2/ARPC5/PABPC1/MARCKS/THY1/PPP1CB/CD9/CSRP1/CDC42/FHL1/PCBP2/CD81/ATP6V0C
    ## GO:0098978                     FLOT2/APOE/CUL3/RTN4/EEA1/PTK2B/PPP2R1A/MAPK3/CADPS/DNM1/VCP/PLXNA4/GABBR1/ARFGEF2/CNKSR2/PPP2R2A/SPTBN1/DLG1/ARF4/HOMER1/SPTB/PCLO/EPHA4/HPCA/ARPC2/GRIA1/NF1/NRCAM/BSN/CPLX1/SYN3/GUCY1B1/SORT1/SCN2A/AP2B1/AP3D1/SPTBN2/BAIAP2/LGI1/ATP2B4/SNX27/FLNA/NAPA/ACTN1/TNR/OGT/CTNNB1/BIN1/PPFIA3/PHB2/CLSTN1/SLC1A2/ATP2B3/PPFIA2/ACTC1/SPARCL1/SNAP25/SH3GL3/GRM5/APBA1/GRIPAP1/SNX6/PLCB1/VGF/CTNND1/WASL/ATP2B1/SYNGAP1/TRIO/YWHAH/SLC30A3/APPL1/NRXN2/ACTN2/CNN3/ITPKA/SEPTIN11/PPP1CA/CADPS2/NPTX1/RAC1/DLG4/GRM2/DBNL/PTPRD/ERC2/EPS15/PFN1/AP2M1/ATP2B2/GRM3/ICAM5/ABR/PPP3R1/PAFAH1B1/SH3GL2/FBXO2/HIP1R/VPS35/VPS18/MAP1LC3A/ADCY1/DOC2A/CNTNAP1/STXBP1/PPM1H/PPP3CB/KCNA2/RAP1A/PRKAR2B/DNM3/HSPA8/GPM6A/PTPRS/PIN1/ARF6/ASAP1/SRCIN1/BCAN/ITGB1/SLC16A7/BCR/NPTXR/CAMKV/SV2A/ACTB/CACNB4/ADAM23/HCN1/PAK3/PLXNC1/GSK3B/FLOT1/FARP1/PPP3CA/SH3GL1/PAK2/DBN1/STK38L/SLC6A17/DLG3/SYN2/RAB8A/DLGAP1/EIF4E/CTBP1/SRC/CPLX2/YWHAZ/PRRT2/SYT1/ARHGAP44/ARPC5L/ATAD1/DNAJB1/PRKAR1A/DLGAP3/RAP1B/PRKACA/HRAS/STX4/RHOA/HNRNPD/CORO1A/NAPB/PURA/ADGRL3/MAL2/FXYD6/CDC42/CNR1/PFN2
    ## GO:0098800                                                                                                                                                                                                                                                                                                                                                                                                                                 BCS1L/NDUFS4/MTX1/NDUFB10/NDUFA12/DNAJC11/AFG3L2/MTX2/TIMM44/NDUFA5/NDUFA13/NDUFB7/ATP5F1B/PHB2/COX1/COX5A/NDUFS1/NDUFB4/MPC2/COX6B1/NDUFA9/ATP5PB/COX7A2L/NDUFA3/UQCRC1/IMMT/UQCRB/NDUFS3/CYTB/CHCHD6/SAMM50/PHB1/ATP5PO/HSPA9/NDUFS2/ND5/COX3/AGK/NDUFB6/SDHB/NDUFV1/COX6C/NDUFB3/SDHA/UQCRC2/COX4I1/SPG7/ND1/COX5B/ATP5F1A/NDUFS8/MICOS13/GRPEL1/MCU/NDUFB9/APOOL/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/NDUFA8/ATP5MF/ATP5F1C/NDUFC2/DNAJC19/TIMM10/ATP5PD/NDUFA10/CHCHD3/NDUFAB1/NDUFA11/NDUFS7/UQCRQ/SLC25A6/COX7A1/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/TIMM9/APOO/NDUFB1/ATP5MK/COX7A2/TIMM50/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/ATP5MJ/NDUFB11
    ## GO:0098798                                                                                                                                                                                                   BCS1L/DBT/IDH3B/NDUFS4/MTX1/NDUFB10/PDHX/NDUFA12/HADHB/DNAJC11/MCCC1/AFG3L2/HADHA/MTX2/TOMM40L/PPIF/TIMM44/NFS1/TIMM13/NDUFA5/NDUFA13/NDUFB7/ATP5F1B/PDHA1/PHB2/COX1/TOMM22/COX5A/SUCLG2/HSD17B10/NDUFS1/IDH3A/DLD/NDUFB4/MPC2/COX6B1/NDUFA9/ATP5PB/IDH3G/COX7A2L/NDUFA3/UQCRC1/IMMT/UQCRB/NDUFS3/CYTB/SUCLG1/CHCHD6/SAMM50/PHB1/ATP5PO/HSPA9/FXN/NDUFS2/MRPS36/VDAC1/ND5/COX3/AGK/TOMM70/NDUFB6/SDHB/NDUFV1/TOMM5/COX6C/NDUFB3/SDHA/DLAT/UQCRC2/COX4I1/LYRM4/SPG7/PDK2/ND1/COX5B/ATP5F1A/NDUFS8/MICOS13/GRPEL1/MCU/PMPCB/BCKDHA/NDUFB9/APOOL/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/NDUFA8/ATP5MF/ATP5F1C/NDUFC2/DNAJC19/TIMM10/ATP5PD/NDUFA10/CHCHD3/SLC25A5/PDHB/NDUFAB1/ISCU/TOMM6/NDUFA11/NDUFS7/TOMM40/UQCRQ/SLC25A6/COX7A1/TOMM20/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/TIMM9/APOO/NDUFB1/SLC25A4/ATP5MK/COX7A2/TIMM50/MRPL12/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/ATP5MJ/NDUFB11
    ## GO:0008021                                                                                                                                                                                                                                                                                                                                                                         DNM1L/MCTP1/COPS5/GABBR1/UNC13A/TMEM163/APP/MFF/DMXL2/ATP6V1H/ATP6V1C1/SYNPR/GRIA1/SNAP91/ATP6V1B2/BSN/SYN3/ATP6V1A/COPS4/SEPTIN6/SEPTIN4/RPH3A/PTPRN2/ATP6V0D1/STX1A/PICALM/BIN1/STX6/ICA1/PPFIA2/SNAP25/APBA1/SYT7/RAB5B/ATP2B1/SLC30A3/RAB11B/ATP6V0A1/ROGDI/MTMR2/VPS45/SYT5/SYT12/ATP6V1D/DLG4/AMPH/RAB12/ATP8A1/NCSTN/IQSEC1/ATP6AP1/BRSK1/DOC2A/SVOP/PPT1/BTBD8/GAD2/PTPRS/PICK1/GRIN1/SV2B/SLC17A6/WFS1/CLTB/SCAMP1/SEPTIN5/TPRG1L/TMED9/SEPTIN2/SEPTIN8/VTI1A/SV2A/SYNGR3/SLC32A1/STXBP5/VAMP2/SLC6A17/RAB3B/SYN2/RAB8A/CPLX3/ATP6V1F/CLTA/RAB27B/STX7/PRRT2/SYT1/RAB7A/NDEL1/ATP6V1E1/ATP6V1G2/TRAPPC4/STX12/SYNGR1/VAMP1/LAMP1/SYN1/SYP/RAB5A/RAB8B/SLC17A7/LAMP5/SYPL1/DNAJC5/MAL2/RAB3C/SCAMP5/RAB3A/SNCA/ATP6V0C
    ##            Count
    ## GO:0005925   183
    ## GO:0030055   183
    ## GO:0098978   172
    ## GO:0098800    98
    ## GO:0098798   134
    ## GO:0008021   110

``` r
head(as.data.frame(ego_bp))
```

    ##                    ID                                         Description
    ## GO:0099003 GO:0099003               vesicle-mediated transport in synapse
    ## GO:0099504 GO:0099504                              synaptic vesicle cycle
    ## GO:0045333 GO:0045333                                cellular respiration
    ## GO:0009060 GO:0009060                                 aerobic respiration
    ## GO:0015980 GO:0015980 energy derivation by oxidation of organic compounds
    ## GO:0006163 GO:0006163                 purine nucleotide metabolic process
    ##            GeneRatio   BgRatio       pvalue     p.adjust       qvalue
    ## GO:0099003  127/2222 220/18870 4.592482e-61 2.652159e-57 1.848112e-57
    ## GO:0099504  116/2222 198/18870 8.665329e-57 2.502114e-53 1.743555e-53
    ## GO:0045333  127/2222 243/18870 2.611563e-54 5.027258e-51 3.503159e-51
    ## GO:0009060  113/2222 197/18870 5.249047e-54 7.578312e-51 5.280818e-51
    ## GO:0015980  144/2222 337/18870 1.270618e-47 1.467564e-44 1.022647e-44
    ## GO:0006163  170/2222 463/18870 3.059672e-45 2.944934e-42 2.052127e-42
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      geneID
    ## GO:0099003                                                                                                                                                                                                                                                       RAPGEF4/CYFIP1/CADPS/DNM1/STX1B/UNC13A/SNX9/ATP2A2/PCLO/ATP6V1H/HPCA/CASK/ATP6V1C1/PLAA/SYNJ1/GAK/SNAP91/ATP6V1B2/BSN/CPLX1/ITSN1/SYN3/CANX/AP2A1/ATP6V1A/AP2B1/CDH2/AP3D1/RPH3A/NAPA/ATP6V0D1/STX1A/PICALM/AP3M2/CTNNB1/BIN1/PRKCB/PPFIA3/CLSTN1/PPFIA2/SNAP25/SH3GL3/GRIPAP1/SYT7/AP2A2/DNM2/ATP6V0A1/RAB3GAP1/VAC14/RIMS1/SYT5/SYT12/CADPS2/ATP6V1D/AMPH/AP3B2/ERC2/EPS15/AP2M1/ATP6AP1/PPP3R1/BRSK1/SH3GL2/VPS35/SNCB/VPS18/DOC2A/GIT1/STXBP1/PPP3CB/BRAF/RAP1A/DNM3/AP1G1/PIP5K1C/ARF6/RALA/PRKCG/SNCG/SV2B/SLC17A6/SEPTIN5/TBC1D24/PREPL/TPRG1L/SCRN1/SV2A/ACTB/CACNB4/SNPH/PTEN/BRSK2/SH3GL1/AP3S1/SLC32A1/AAK1/DNAJC6/STXBP5/VAMP2/RAB3B/SYN2/RAB8A/CPLX3/ATP6V1F/CTBP1/RAB27B/PACSIN1/CPLX2/PRRT2/SYT1/ATAD1/ATP6V1E1/ARHGDIA/ATP6V1G2/CDK5/SNAP29/RAP1B/NAPB/SYN1/SYP/AP2S1/RAB5A/SLC17A7/DNAJC5/RAB3A/PFN2/SNCA
    ## GO:0099504                                                                                                                                                                                                                                                                                                                       RAPGEF4/CYFIP1/CADPS/DNM1/STX1B/UNC13A/SNX9/ATP2A2/PCLO/ATP6V1H/CASK/ATP6V1C1/PLAA/SYNJ1/GAK/SNAP91/ATP6V1B2/BSN/CPLX1/ITSN1/SYN3/CANX/ATP6V1A/AP2B1/CDH2/AP3D1/RPH3A/NAPA/ATP6V0D1/STX1A/PICALM/AP3M2/CTNNB1/BIN1/PRKCB/PPFIA3/PPFIA2/SNAP25/SH3GL3/GRIPAP1/SYT7/DNM2/ATP6V0A1/RAB3GAP1/RIMS1/SYT5/SYT12/CADPS2/ATP6V1D/AMPH/AP3B2/ERC2/AP2M1/ATP6AP1/PPP3R1/BRSK1/SH3GL2/SNCB/VPS18/DOC2A/GIT1/STXBP1/PPP3CB/BRAF/RAP1A/DNM3/AP1G1/PIP5K1C/ARF6/PRKCG/SNCG/SV2B/SLC17A6/SEPTIN5/TBC1D24/PREPL/TPRG1L/SCRN1/SV2A/ACTB/CACNB4/SNPH/PTEN/BRSK2/SH3GL1/AP3S1/SLC32A1/DNAJC6/STXBP5/VAMP2/RAB3B/SYN2/CPLX3/ATP6V1F/CTBP1/RAB27B/PACSIN1/CPLX2/PRRT2/SYT1/ATP6V1E1/ARHGDIA/ATP6V1G2/CDK5/SNAP29/RAP1B/NAPB/SYN1/SYP/AP2S1/RAB5A/SLC17A7/DNAJC5/RAB3A/PFN2/SNCA
    ## GO:0045333                                                                                                                                                                                                                                                          MTFR1L/IDH3B/NDUFS4/NDUFB10/NDUFA12/PLEC/VCP/OGDHL/SLC25A13/SUCLA2/PPIF/MDH1/ACLY/NDUFA5/NDUFA13/NDUFB7/UQCC2/ATP5F1B/LYRM7/PDHA1/DLST/GHITM/COX1/ME3/COX5A/SUCLG2/SLC25A12/SOD2/NDUFS1/MDH2/IDH3A/CAT/DLD/NDUFB4/IDH2/COX6B1/NDUFA9/ATP5PB/IDH3G/NNT/COX7A2L/NDUFA3/UQCRC1/IDH1/ACO1/TRAP1/UQCRB/NDUFS3/CYTB/SUCLG1/COQ9/ATP5PO/OGDH/FXN/GPD2/STOML2/NDUFS2/GBA1/MRPS36/ND5/COX3/AK4/NDUFB6/SDHB/NDUFV1/COX6C/NDUFB3/SDHA/CS/ACO2/DLAT/ETFDH/UQCRC2/COX4I1/ND1/PARK7/COX5B/ATP5F1A/NDUFS8/ETFB/FH/SLC25A22/SLC25A25/NDUFB9/ARL2/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/NIPSNAP2/SIRT3/ETFA/NDUFA8/ATP5MF/CISD1/ATP5F1C/NDUFC2/RHOA/ATP5PD/NDUFA10/PDHB/NDUFAB1/ISCU/NDUFA11/NDUFS7/CYCS/UQCRQ/COX7A1/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/COX7A2/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/SNCA/NDUFB11
    ## GO:0009060                                                                                                                                                                                                                                                                                                                                                    MTFR1L/IDH3B/NDUFS4/NDUFB10/NDUFA12/VCP/OGDHL/SUCLA2/PPIF/MDH1/ACLY/NDUFA5/NDUFA13/NDUFB7/UQCC2/ATP5F1B/PDHA1/DLST/GHITM/COX1/ME3/COX5A/SUCLG2/NDUFS1/MDH2/IDH3A/CAT/DLD/NDUFB4/IDH2/COX6B1/NDUFA9/ATP5PB/IDH3G/NNT/COX7A2L/NDUFA3/UQCRC1/IDH1/ACO1/UQCRB/NDUFS3/CYTB/SUCLG1/COQ9/ATP5PO/OGDH/FXN/STOML2/NDUFS2/MRPS36/ND5/COX3/AK4/NDUFB6/SDHB/NDUFV1/COX6C/NDUFB3/SDHA/CS/ACO2/DLAT/UQCRC2/COX4I1/ND1/PARK7/COX5B/ATP5F1A/NDUFS8/FH/NDUFB9/ARL2/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/NIPSNAP2/SIRT3/NDUFA8/ATP5MF/ATP5F1C/NDUFC2/RHOA/ATP5PD/NDUFA10/PDHB/NDUFAB1/ISCU/NDUFA11/NDUFS7/CYCS/UQCRQ/COX7A1/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/COX7A2/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/SNCA/NDUFB11
    ## GO:0015980                                                                                                                                                   MTFR1L/IDH3B/NDUFS4/NDUFB10/NDUFA12/PLEC/VCP/OGDHL/SLC25A13/PFKM/SUCLA2/PYGB/PPIF/GABARAPL1/MDH1/ACLY/ACADVL/PTGES3/NDUFA5/NDUFA13/NDUFB7/GSK3A/UQCC2/PPP1R2/ATP5F1B/UGP2/LYRM7/PDHA1/DLST/GHITM/COX1/ME3/COX5A/SUCLG2/SLC25A12/SOD2/NDUFS1/MDH2/PYGM/GYS1/IDH3A/CAT/DLD/PPP1CA/NDUFB4/IDH2/COX6B1/NDUFA9/ATP5PB/IDH3G/NNT/COX7A2L/NDUFA3/UQCRC1/IDH1/ACO1/TRAP1/UQCRB/NDUFS3/CYTB/SUCLG1/COQ9/AGL/ATP5PO/OGDH/FXN/GPD2/STOML2/NDUFS2/GBA1/PHKB/MRPS36/ND5/COX3/AK4/NDUFB6/SDHB/NDUFV1/COX6C/NDUFB3/SDHA/CS/ACO2/GSK3B/DLAT/ETFDH/UQCRC2/COX4I1/ND1/PARK7/COX5B/ATP5F1A/SORBS1/NDUFS8/ETFB/FH/SLC25A22/SLC25A25/NDUFB9/ARL2/ATP5ME/ATP5MG/NDUFV3/GBE1/ND4/CYC1/NDUFV2/NIPSNAP2/SIRT3/ETFA/NDUFA8/ATP5MF/CISD1/ATP5F1C/NDUFC2/RHOA/ATP5PD/NDUFA10/PDHB/NDUFAB1/ISCU/NDUFA11/NDUFS7/CYCS/UQCRQ/COX7A1/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/COX7A2/PPP1CB/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/SNCA/NDUFB11
    ## GO:0006163 GMPS/AK3/TJP2/NDUFS4/NME2/NDUFB10/NAXD/PDHX/PRPSAP2/NDUFA12/GART/DNM1L/VCP/SLC25A13/NADK2/DLG1/AK5/PC/NT5E/CASK/GDA/SUCLA2/PANK4/ME1/MTHFD1/ATP6V1B2/HSPA1A/HSPA1B/KCNAB2/GUCY1B1/MDH1/ACLY/ATP6V1A/NUDT2/TKT/NDUFA5/GUK1/RAN/NDUFA13/NDUFB7/ADSS2/OPA1/ATP5F1B/ACSL6/PDHA1/DLST/ACSS1/HSD17B4/DLG2/OLA1/PGAM1/IMPDH2/AK1/SUCLG2/ATIC/PGLS/ACSBG1/NDUFS1/RPTOR/MDH2/PRPSAP1/ATP1A2/PMVK/DLD/PDK3/NDUFB4/MPC2/IDH2/ACSF3/NDUFA9/DCXR/ATP5PB/NNT/NDUFA3/AMPD2/G6PD/ALDH1L1/APRT/IDH1/NDUFS3/CD38/SUCLG1/ADCY9/PDE2A/ACOT1/PAICS/NUDT3/ADCY1/ATP5PO/OGDH/ACSL3/PPT1/STOML2/NDUFS2/ME2/HSPA8/PDE1A/ACAT1/PGD/MPP1/ND5/HSD17B12/TECR/ACACA/AK4/NDUFB6/SDHB/NDUFV1/NDUFB3/ACOT7/NAMPT/CACNB4/SDHA/PPCS/DLAT/ACOT9/ENO1/PDK2/ND1/COASY/NME1/TALDO1/ATP5F1A/HINT1/NDUFS8/RAB23/ADCY2/GPD1L/SLC25A25/NDUFB9/ACSL1/ATP5ME/ATP5MG/NDUFV3/ND4/NDUFV2/NDUFA8/ATP5MF/ALDOA/ATP5F1C/NDUFC2/LDHB/ATP5PD/NDUFA10/PDHB/NDUFAB1/HPRT1/NDUFA11/SLC25A1/NDUFS7/PRPS1/NDUFB8/NDUFS6/ATP1B1/ATP5PF/NDUFS5/NDUFA7/NDUFB5/NDUFB1/ATP5MK/MMUT/ATP6/ATP5F1D/NDUFA2/ATP5IF1/ATP5MJ/SNCA/NDUFB11/FIS1/ATP6V0C
    ##            Count
    ## GO:0099003   127
    ## GO:0099504   116
    ## GO:0045333   127
    ## GO:0009060   113
    ## GO:0015980   144
    ## GO:0006163   170

``` r
head(as.data.frame(reactome_res))
```

    ##                          ID
    ## R-HSA-1428517 R-HSA-1428517
    ## R-HSA-112315   R-HSA-112315
    ## R-HSA-163200   R-HSA-163200
    ## R-HSA-611105   R-HSA-611105
    ## R-HSA-112316   R-HSA-112316
    ## R-HSA-6798695 R-HSA-6798695
    ##                                                                                                                       Description
    ## R-HSA-1428517                                                      The citric acid (TCA) cycle and respiratory electron transport
    ## R-HSA-112315                                                                                Transmission across Chemical Synapses
    ## R-HSA-163200  Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.
    ## R-HSA-611105                                                                                       Respiratory electron transport
    ## R-HSA-112316                                                                                                      Neuronal System
    ## R-HSA-6798695                                                                                            Neutrophil degranulation
    ##               GeneRatio   BgRatio       pvalue     p.adjust       qvalue
    ## R-HSA-1428517  113/1733 178/11009 8.188571e-48 1.119378e-44 6.481901e-45
    ## R-HSA-112315   131/1733 270/11009 1.898306e-37 1.297492e-34 7.513295e-35
    ## R-HSA-163200    80/1733 127/11009 9.743671e-34 4.439866e-31 2.570961e-31
    ## R-HSA-611105    67/1733 103/11009 1.117194e-29 3.818011e-27 2.210868e-27
    ## R-HSA-112316   155/1733 410/11009 1.366200e-28 3.735192e-26 2.162911e-26
    ## R-HSA-6798695  172/1733 480/11009 1.918157e-28 4.370200e-26 2.530621e-26
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            geneID
    ## R-HSA-1428517                                                                                                                                                                                                                                                                                                                                                  IDH3B/NDUFS4/NDUFB10/PDHX/NDUFA12/NDUFAF3/GLO1/SUCLA2/LRPPRC/ME1/LDHA/BSG/NDUFA5/NDUFA13/NDUFB7/ATP5F1B/PDHA1/DLST/COX1/HAGH/ME3/COX5A/SUCLG2/NDUFS1/MDH2/IDH3A/DLD/PDK3/NDUFB4/NDUFAF4/MPC2/IDH2/COX6B1/NDUFA9/ATP5PB/IDH3G/NNT/COX7A2L/NDUFA3/UQCRC1/TRAP1/UQCRB/NDUFS3/CYTB/SUCLG1/ATP5PO/OGDH/NDUFS2/ME2/VDAC1/ND5/COX3/NDUFB6/SDHB/NDUFAF2/NDUFV1/COX6C/NDUFB3/SLC16A1/SDHA/CS/ACO2/DLAT/ETFDH/UQCRC2/COX4I1/PDK2/ND1/SLC25A27/ACAD9/SCO1/COX5B/ATP5F1A/NDUFS8/ETFB/FH/FAHD1/NDUFB9/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/ETFA/NDUFA8/ATP5MF/ATP5F1C/NDUFC2/LDHB/ATP5PD/NDUFA10/PDHB/NDUFAB1/NDUFA11/NDUFS7/CYCS/UQCRQ/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/NDUFB11
    ## R-HSA-112315                                                                                                                                                                                                                                 MYO6/GAD1/GNG5/RPS6KA2/MAPK3/CAMK4/GABBR1/PRKAR2A/ALDH5A1/DLG1/GLS/CASK/CACNA2D1/GRIA1/CACNA2D2/PRKCA/CAMK2D/CPLX1/SYN3/AP2A1/AP2B1/CAMK2G/MAPT/STX1A/NBEA/PRKCB/PPFIA3/EPB41L1/SLC1A2/GABBR2/PPFIA2/DLG2/SNAP25/APBA1/MAOA/PLCB1/AKAP5/AP2A2/GABRA3/ACTN2/GRIA3/GABRB2/GNB5/RIMS1/RAC1/DLG4/GNG13/CACNA2D3/ABAT/AP2M1/PRKAA2/PRKAA1/TUBA8/PDPK1/ADCY9/ADCY1/COMT/GIT1/STXBP1/GRIA2/GABRB1/PRKAR2B/GAD2/CAMK1/HSPA8/PICK1/GNG10/RASGRF2/NEFL/GABRA1/GRIN1/PRKAB2/PRKCG/ARHGEF7/SLC6A11/GNAI2/SLC1A3/CACNB4/PRKACB/NSF/NCALD/GNG12/GABRG2/MAPK1/SLC32A1/CALM1/CALM2/CALM3/CAMKK2/GRIN2B/TUBA1A/CAMKK1/VAMP2/TUBB4A/CACNG8/DLG3/GNB1/SYN2/PPFIA4/LIN7A/ADCY2/SRC/SLC6A1/TSPAN7/SYT1/GABRB3/GLUL/TUBB3/PRKAR1A/CAMK2B/GNG3/TUBA4A/PRKACA/CACNA1B/HRAS/ALDH2/NPTN/SYN1/KRAS/LIN7B/AP2S1/GNAI1/GNG7/GNG2/TUBB6/SLC17A7/DNAJC5/CAMK2A/ARL6IP5/TUBAL3/RAB3A
    ## R-HSA-163200                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             NDUFS4/NDUFB10/NDUFA12/NDUFAF3/LRPPRC/NDUFA5/NDUFA13/NDUFB7/ATP5F1B/COX1/COX5A/NDUFS1/NDUFB4/NDUFAF4/COX6B1/NDUFA9/ATP5PB/COX7A2L/NDUFA3/UQCRC1/TRAP1/UQCRB/NDUFS3/CYTB/ATP5PO/NDUFS2/ND5/COX3/NDUFB6/SDHB/NDUFAF2/NDUFV1/COX6C/NDUFB3/SDHA/ETFDH/UQCRC2/COX4I1/ND1/SLC25A27/ACAD9/SCO1/COX5B/ATP5F1A/NDUFS8/ETFB/NDUFB9/ATP5ME/ATP5MG/NDUFV3/ND4/CYC1/NDUFV2/ETFA/NDUFA8/ATP5MF/ATP5F1C/NDUFC2/ATP5PD/NDUFA10/NDUFAB1/NDUFA11/NDUFS7/CYCS/UQCRQ/NDUFB8/NDUFS6/ATP5PF/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/NDUFA4/ATP6/ATP5F1D/NDUFA2/UQCR10/NDUFB11
    ## R-HSA-611105                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            NDUFS4/NDUFB10/NDUFA12/NDUFAF3/LRPPRC/NDUFA5/NDUFA13/NDUFB7/COX1/COX5A/NDUFS1/NDUFB4/NDUFAF4/COX6B1/NDUFA9/COX7A2L/NDUFA3/UQCRC1/TRAP1/UQCRB/NDUFS3/CYTB/NDUFS2/ND5/COX3/NDUFB6/SDHB/NDUFAF2/NDUFV1/COX6C/NDUFB3/SDHA/ETFDH/UQCRC2/COX4I1/ND1/ACAD9/SCO1/COX5B/NDUFS8/ETFB/NDUFB9/NDUFV3/ND4/CYC1/NDUFV2/ETFA/NDUFA8/NDUFC2/NDUFA10/NDUFAB1/NDUFA11/NDUFS7/CYCS/UQCRQ/NDUFB8/NDUFS6/UQCRFS1/NDUFS5/COX2/NDUFA7/NDUFB5/NDUFB1/NDUFA4/NDUFA2/UQCR10/NDUFB11
    ## R-HSA-112316                                                                             FLOT2/MYO6/GAD1/GNG5/RPS6KA2/MAPK3/CAMK4/GABBR1/PRKAR2A/ALDH5A1/DLG1/GLS/HOMER1/CASK/CACNA2D1/NRXN3/GRIA1/EPB41L2/CACNA2D2/PRKCA/CAMK2D/CPLX1/SYN3/KCNAB2/AP2A1/EPB41L3/AP2B1/CAMK2G/MAPT/NRXN1/STX1A/NBEA/PRKCB/PPFIA3/EPB41L1/SLC1A2/GABBR2/PPFIA2/DLG2/SNAP25/GRM5/APBA1/SYT7/MAOA/PLCB1/AKAP5/AP2A2/GABRA3/NRXN2/ACTN2/GRIA3/GABRB2/GNB5/RIMS1/SYT12/RAC1/DLG4/GNG13/DBNL/CACNA2D3/PTPRD/ABAT/AP2M1/PRKAA2/PRKAA1/TUBA8/PDPK1/KCNMA1/ADCY9/ADCY1/COMT/GIT1/STXBP1/GRIA2/GABRB1/RTN3/KCNA2/PRKAR2B/GAD2/CAMK1/HSPA8/PTPRS/PICK1/GNG10/RASGRF2/NEFL/GABRA1/GRIN1/PRKAB2/PRKCG/ARHGEF7/SLC6A11/KCNN3/GNAI2/SLC1A3/CACNB4/HCN1/PRKACB/NSF/FLOT1/NCALD/GNG12/GABRG2/MAPK1/SLC32A1/CALM1/CALM2/CALM3/CAMKK2/GRIN2B/TUBA1A/CAMKK1/VAMP2/TUBB4A/CACNG8/DLG3/GNB1/SYN2/PPFIA4/LIN7A/DLGAP1/ADCY2/SRC/SLC6A1/TSPAN7/SYT1/GABRB3/GLUL/TUBB3/PRKAR1A/CAMK2B/GNG3/DLGAP3/TUBA4A/PRKACA/CACNA1B/HRAS/ALDH2/KCND3/KCNA4/NPTN/SYN1/KRAS/LIN7B/AP2S1/GNAI1/GNG7/GNG2/TUBB6/SLC17A7/DNAJC5/CAMK2A/ARL6IP5/TUBAL3/RAB3A
    ## R-HSA-6798695 PFKL/DYNC1H1/NME2/NME1-NME2/SPTAN1/CYFIP1/KRT1/VCP/ANXA2/CAPN1/RAB6A/TOM1/ACAA1/PSMD7/GCA/PYGB/CST3/PSMD13/TMEM30A/SLC44A2/ITGAV/UBR4/HSPA1A/HSPA1B/KCNAB2/PSMD12/APEH/ACLY/CD47/PSMA5/GDI2/PADI2/CAND1/PTPRN2/LTA4H/EEF1A1/PSMD3/SNAP25/ALDOC/DYNC1LI1/RAB5C/RAB5B/CCT8/VPS35L/AP2A2/PGAM1/IMPDH2/STOM/ACTR10/HUWE1/PSMD1/HMOX2/C3/SLC2A5/ATP6V0A1/TOLLIP/VCL/CAT/RAC1/ATP6V1D/ATG7/DBNL/ATP8A1/NCSTN/PA2G4/CD44/ERP44/CAB39/PRKCD/PKM/FTH1/SERPINA1/APRT/IDH1/ACTR2/LGALS3/PSMD2/AGL/DNAJC13/VAT1/LAMTOR1/TMT1A/EEF2/RAP1A/CCT2/PSMD11/VAPA/PRDX6/GPI/KPNB1/HMGB1/HSPA8/PGM1/RAB10/S100A9/CAP1/IST1/HSP90AA1/PSMC2/ARMC8/CPNE3/CYB5R3/SCAMP1/S100A8/AP1M1/SDCBP/PSMD14/COPB1/PSMA2/RAB4B/RAP2B/PSAP/GSN/FLG2/SLC2A3/SERPINB6/LAMP2/S100A11/SIRPA/RAB14/MAPK1/COTL1/FABP5/NFASC/RAB18/PAFAH1B2/MLEC/RAB24/PSMC3/DIAPH1/ALAD/SERPINA3/CSNK2B/RAB7A/ACTR1B/PSMB1/PDXK/CD63/TRAPPC1/ASAH1/SYNGR1/LAMTOR3/GSTP1/ARL8A/LAMP1/SNAP29/RAP1B/RHOG/ALDOA/PPIA/PSMD6/HSP90AB1/PTGES2/NDUFC2/RHOA/PGRMC1/CTSB/NIT2/B2M/FTL/RAB9B/TUBB/ARPC5/FAF2/CPNE1/MIF/CTSD/DNAJC5/RAB3A/SVIP/HBB/ATP6V0C
    ##               Count
    ## R-HSA-1428517   113
    ## R-HSA-112315    131
    ## R-HSA-163200     80
    ## R-HSA-611105     67
    ## R-HSA-112316    155
    ## R-HSA-6798695   172

#### Dotplots

``` r
dotplot(ego_cc, showCategory = 15) +
  ggtitle("GO Cellular Component Enrichment")
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-55-1.png" alt="" style="display: block; margin: auto;" />

``` r
dotplot(reactome_res, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment")
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-55-2.png" alt="" style="display: block; margin: auto;" />

``` r
ggsave(
  "GO_CC_synaptic_enrichment.png",
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "Reactome_synaptic_enrichment.png",
  width = 10,
  height = 8,
  dpi = 300
)
```

### 17.1 PathfindR Installation and Input Preparation

``` r
# 1. Install necessary packages if not already present
# install.packages("remotes")
# remotes::install_github("noriakis/ggkegg")
# remotes::install_github("egeulgen/pathfindR")

library(pathfindR)
library(dplyr)

# 2. Prepare input from the limma results (top_table)
# We use the p-values and logFC adjusted for Age, Sex, and PMI
pathfindr_input <- top_table %>%
  mutate(Protein = rownames(.)) %>%
  dplyr::select(Gene.symbol = Protein, 
                logFC = logFC,   
                p_val = P.Value,
                adj_p_val = adj.P.Val)

# Optional: Filter for significant proteins to focus the subnetwork search
# pathfindr_input <- pathfindr_input %>% filter(p_val < 0.05)
```

### 17.2 Gene Mapping and Active Subnetwork Processing

``` r
# Load Mapping file: UniProt Accession → Gene Symbols (HGNC)
gene_mapping <- read.delim("Uniprot_GeneName_Crack.txt", header = TRUE, stringsAsFactors = FALSE)

# Process input for pathfindR compatibility
pathfindr_input_net <- pathfindr_input %>%
  # Clean p-values and handle NAs
  mutate(P_VALUE = as.numeric(trimws(p_val))) %>%
  mutate(P_VALUE = ifelse(is.na(P_VALUE), 1, P_VALUE)) %>%
  
  # Join with mapping file to convert UniProt IDs to Symbols
  left_join(gene_mapping, by = c("Gene.symbol" = "Entry")) %>%
  # Use mapped Gene Name, fallback to UniProt ID if not found
  mutate(Gene.symbol = ifelse(!is.na(Genes), sub(" .*", "", Genes), Gene.symbol)) %>%
  
  dplyr::select(Gene.symbol, logFC, P_VALUE)

# Final verification of P-values
print(paste("Missing P-values:", sum(is.na(pathfindr_input_net$P_VALUE))))
```

    ## [1] "Missing P-values: 0"

``` r
# Processing with pathfindR
# We set a higher p-val threshold to allow more genes to connect STRING sub-networks
processed_pathfindR <- input_processing(
  input = pathfindr_input_net,
  p_val_threshold = 0.1,      # Genes with p < 0.1 will be used for subnetwork search
  convert2alias = TRUE       # Helps find genes with outdated nomenclature
)

# Run Enrichment (GO-All and GO-CC)
output_df_all <- run_pathfindR(processed_pathfindR[,c(2,4)], gene_sets = "GO-All")
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-57-1.png" alt="" style="display: block; margin: auto;" />

``` r
output_df_CC <- run_pathfindR(processed_pathfindR[,c(2,4)], gene_sets = "GO-CC")
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-57-2.png" alt="" style="display: block; margin: auto;" />

### 17.3 g:Profiler Enrichment and Semantic Similarity

``` r
library(gprofiler2)

# Extract gene names for query
query_genes <- pathfindr_input$Gene.symbol

# Run g:Profiler
enr_results <- gost(
    query = query_genes, 
    organism = "hsapiens",
    ordered_query = TRUE,       
    significant = TRUE,
    evcodes = TRUE, 
    correction_method = "fdr",
    sources = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC")
)

# Visualization
if(!is.null(enr_results)){
    gostplot(enr_results, interactive = FALSE)
}
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-58-1.png" alt="" style="display: block; margin: auto;" />

``` r
# Extract Biological Process results
res_table <- enr_results$result
go_bp_results <- res_table[res_table$source == "GO:BP" & res_table$p_value < 0.05, ]
```

### 17.4 RRVGO: Semantic Similarity and Redundancy Reduction

``` r
# install.packages("BiocManager")
# BiocManager::install("rrvgo")
library(rrvgo)

# Prepare scores (Negative Log10 P-value)
go_ids <- go_bp_results$term_id
go_scores <- -log10(go_bp_results$p_value)
names(go_scores) <- go_ids

# 1. Calculate Semantic Similarity Matrix
sim_matrix <- calculateSimMatrix(
  go_ids,
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel"
)

# 2. Reduce redundant terms (Threshold 0.7 for clustering)
reduced_terms <- reduceSimMatrix(
  sim_matrix,
  go_scores,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)

# 3. Visualizations
treemapPlot(reduced_terms)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-59-1.png" alt="" style="display: block; margin: auto;" />

``` r
scatterPlot(sim_matrix, reduced_terms)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-59-2.png" alt="" style="display: block; margin: auto;" />

``` r
# 4. Extract representative Non-Redundant Pathways
final_representative_terms <- reduced_terms[reduced_terms$parentTerm == reduced_terms$term, ]
final_representative_terms <- final_representative_terms[order(final_representative_terms$score, decreasing = TRUE), ]

head(final_representative_terms[, c("term", "score", "size")], 10)
```

    ##                                                     term    score size
    ## GO:0051179                                  localization 91.11221 5320
    ## GO:0051641                         cellular localization 69.78343 3433
    ## GO:0016192                    vesicle-mediated transport 65.61130 1536
    ## GO:0006996                        organelle organization 62.59412 3613
    ## GO:0008104                          protein localization 62.59412 2499
    ## GO:0032879                    regulation of localization 54.04750 2147
    ## GO:0009060                           aerobic respiration 53.73621  197
    ## GO:0051128 regulation of cellular component organization 51.44899 2451
    ## GO:0065008              regulation of biological quality 43.60404 2985
    ## GO:0061024                         membrane organization 41.28572  815

### 17.5. clusterProfiler: Network Visualization

``` r
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Convert UniProt to Entrez ID for clusterProfiler
gene_conversion <- bitr(
  query_genes,
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO:BP Enrichment
ego_bp <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Semantic similarity for plotting
ego_bp_sim <- pairwise_termsim(ego_bp)
emapplot(ego_bp_sim)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-60-1.png" alt="" style="display: block; margin: auto;" />

``` r
# Simplify GO terms (automated reduction)
ego_bp_simplified <- clusterProfiler::simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
top10_BP_table <- as.data.frame(ego_bp_simplified)[1:10, ]

# GO:CC Enrichment
ego_cc <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

ego_cc_sim <- pairwise_termsim(ego_cc)
emapplot(ego_cc_sim)
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-60-2.png" alt="" style="display: block; margin: auto;" />

``` r
# Top 10 Cellular Component terms
ego_cc_simplified <- clusterProfiler::simplify(
  ego_cc,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
top10_CC_table <- as.data.frame(ego_cc_simplified)[1:10, ]
```

### 17.6. Obtain Reference Gene Sets

``` r
# Fetch gene sets for downstream custom analysis
reference_go_list <- fetch_gene_set(
  gene_sets = "GO-All", 
  min_gset_size = 10, 
  max_gset_size = 300
)

go_gene_sets <- reference_go_list[[1]]
go_descriptions <- reference_go_list[[2]]
```

### 17.8 Combined Enrichment and Abundance Visualization

``` r
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

# 1. Calculate Fold Enrichment 
# Converting GeneRatio and BgRatio from strings (e.g., "5/100") to numeric values
top10_CC <- top10_CC_table %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
                    as.numeric(sub(".*/", "", GeneRatio)),

    BgRatio_num = as.numeric(sub("/.*", "", BgRatio)) /
                  as.numeric(sub(".*/", "", BgRatio)),

    FoldEnrichment = GeneRatio_num / BgRatio_num
  )

# --- 2. Color Palette Selection (Nature-inspired) ---
nature_palette <- c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6")

# --- 3. Gene ID Mapping (ENTREZ to UNIPROT) ---
# Expand grouped gene IDs into individual rows for mapping
top10_expanded <- top10_CC %>%
  separate_rows(geneID, sep = "/")

# Convert ENTREZID to UNIPROT for synchronization with protein data
mapping_table <- bitr(top10_expanded$geneID, 
                      fromType = "ENTREZID", 
                      toType = "UNIPROT", 
                      OrgDb = org.Hs.eg.db)

# Join mapping with ontology data
top10_with_uniprot <- top10_expanded %>%
  inner_join(mapping_table, by = c("geneID" = "ENTREZID"))

# --- 4. Prepare Abundance Data for Boxplots ---
# Aggregate protein intensities to get a single mean value per protein
abundance_summary <- processed_final$FeatureLevelData %>%
  group_by(PROTEIN) %>%
  summarise(Mean_Abundance = mean(ABUNDANCE, na.rm = TRUE)) %>%
  ungroup()

# Merge ontology terms with summarized abundance values
data_to_plot <- top10_with_uniprot %>%
  inner_join(abundance_summary, by = c("UNIPROT" = "PROTEIN")) %>%
  filter(Mean_Abundance > 0) # Remove zero-values to clean the plot

# --- 5. Generate Top Panel: Fold Enrichment Bar Chart (p1) ---
p1 <- ggplot(top10_CC, aes(x = reorder(Description, -as.numeric(p.adjust)), 
                            y = FoldEnrichment, fill = p.adjust)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  # Labels showing the count of proteins per term
  geom_text(aes(label = Count), vjust = -1, size = 3) +
  scale_fill_gradientn(
    colours = nature_palette,
    trans = "reverse",
    name = "FDR", 
    labels = function(x) format(x, scientific = TRUE, digits = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), 
                     breaks = c(0, 2.5, 5, 7.5)) +
  labs(y = "Fold Enrichment", x = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),    # Remove X text (it will be shown in p2)
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm"),
    plot.margin = margin(b = -2, t = 10, r = 10, l = 10) # Connect with p2
  )

# --- 6. Generate Bottom Panel: Protein Abundance Boxplot (p2) ---
p2 <- ggplot(data_to_plot, aes(x = reorder(Description, -as.numeric(p.adjust)), 
                               y = Mean_Abundance, fill = p.adjust)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black", linewidth = 0.3) +
  geom_jitter(color = "#4682B4", size = 0.6, alpha = 0.3, width = 0.2) +
  theme_minimal(base_size = 11) +
  scale_fill_gradientn(
    colours = nature_palette,
    trans = "reverse",
    name = "FDR"
  ) +
  scale_y_continuous(limits = c(0, max(data_to_plot$Mean_Abundance) + 5)) +
  # Wrap long pathway names for better readability
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "GO Cellular Component Terms", y = "Protein Abundance") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
    axis.title.x = element_text(size = 10, vjust = -1),
    axis.title.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    legend.position = "none",        # FDR legend is already in p1
    plot.margin = margin(t = -2, b = 10, r = 10, l = 10) 
  )

# --- 7. Final Assembly and Export ---
# Using patchwork to stack plots vertically
plot_final <- (p1 / p2) + plot_layout(heights = c(1, 1.2))

# Save the final publication-ready figure
ggsave("GO_Enrichment_Abundance_Combined.png", plot_final, 
       width = 200, height = 180, units = "mm", dpi = 300)
plot_final
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-62-1.png" alt="" style="display: block; margin: auto;" />

``` r
print("Figure successfully generated and saved as GO_Enrichment_Abundance_Combined.png")
```

    ## [1] "Figure successfully generated and saved as GO_Enrichment_Abundance_Combined.png"

## 18. Data Preparation and Group Logic

This block consolidates the group classification and synchronization
with your previous ontology mapping.

``` r
library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)

# 1. Aseguramos que usamos los datos por PROTEÍNA (más limpio para medias)
# Usamos tu 'top10_with_uniprot' para filtrar los datos resumidos
data_para_grafico <- processed_final$ProteinLevelData %>%
  # Pasamos Protein a carácter por si acaso
  mutate(Protein = as.character(Protein)) %>%
  # Filtramos solo muestras Control y PART
  filter(GROUP %in% c("Control", "PART", "CN_PART")) %>%
  # Unimos con tu tabla de conversión
  inner_join(top10_with_uniprot, by = c("Protein" = "UNIPROT"))

# 2. Generamos el gráfico de barras por muestra
ggplot(data_para_grafico, aes(x = originalRUN, y = LogIntensities, fill = GROUP)) +
  # Calculamos la media de todas las proteínas de esa vía en esa muestra
  stat_summary(fun = mean, geom = "bar", color = "black", size = 0.2, alpha = 0.7) +
  # Añadimos puntos pequeños para ver la dispersión de las proteínas
  geom_jitter(size = 0.4, alpha = 0.3, width = 0.2) +
  # Un panel por cada una de las 10 categorías CC
  facet_wrap(~ str_wrap(Description, width = 30), scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Control" = "#00BFC4", "PART" = "#F8766D")) +
  theme_bw(base_size = 10) +
  labs(
    title = "Abundancia Media por Muestra: Control vs PART",
    subtitle = "Análisis basado en Componente Celular (CC)",
    x = "Muestras (Samples)",
    y = "Abundancia Media (Log2 Intensity)",
    fill = "Grupo"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom"
  )
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-63-1.png" alt="" style="display: block; margin: auto;" />

## 19. Detailed abundance profiles for the top 10 cellular component ontologies

``` r
library(ggplot2)
library(stringr)
library(dplyr)

# 1. We define the list in the exact order (from 1 to 18)

#This list associates the row/run number with its actual name
muestras_lista <- c(
  "Intensity.AD3",   # 1
  "Intensity.AD6",   # 2
  "Intensity.AD11",  # 3
  "Intensity.AD17",  # 4
  "Intensity.AD2",   # 5
  "Intensity.AD7",   # 6
  "Intensity.PT9",   # 7
  "Intensity.AD14",  # 8
  "Intensity.AD15",  # 9
  "Intensity.PT1",   # 10
  "Intensity.CN4",   # 11
  "Intensity.PT5",   # 12
  "Intensity.CN8",   # 13
  "Intensity.CN10",  # 14
  "Intensity.PT12",  # 15
  "Intensity.CN13",  # 16
  "Intensity.PT16",  # 17
  "Intensity.CN18"   # 18
)

# 2. We prepare the data using the numeric index from originalRUN.
data_para_grafico <- processed_final$ProteinLevelData %>%
  mutate(
    # convert originalRUN to a numeric value to use it as a list index
    ID_numerico = as.numeric(as.character(originalRUN)),
    # obtain the full name (e.g., Intensity.CN18)
    FullName = muestras_lista[ID_numerico],
    # clean the name for the plot (e.g., Control 18)
    SampleLabel = FullName %>% 
      str_replace("Intensity.CN", "Control ") %>% 
      str_replace("Intensity.PT", "PART ") %>%
      str_replace("Intensity.AD", "AD "),
    # assign the final group.
    Group = case_when(
      str_detect(FullName, "CN") ~ "Control",
      str_detect(FullName, "PT") ~ "PART",
      TRUE ~ "AD"
    )
  ) %>%
  # filter only Control and PART
  filter(Group %in% c("Control", "PART")) %>%
  inner_join(top10_with_uniprot, by = c("Protein" = "UNIPROT"))

# 3. plot
ggplot(data_para_grafico, aes(x = SampleLabel, y = LogIntensities, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.3, alpha = 0.8) +
  geom_jitter(size = 0.4, alpha = 0.4, width = 0.2, color = "#222222") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.4) +
  
  facet_wrap(~ str_wrap(Description, width = 25), scales = "free_y", ncol = 2) +
  
  # colors
  scale_fill_manual(values = c("Control" = "#EBC4E1", "PART" = "#B57FA4")) +
  
  labs(x = NULL, y = "Mean Abundance", fill = "Group") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )
```

<img src="TFG_Mar_Mauri_files/figure-gfm/unnamed-chunk-64-1.png" alt="" style="display: block; margin: auto;" />

``` r
# Assuming you saved your ggplot in an object called ‘p_final’

#If not, just run the line after your ggplot code

ggsave(
  filename = "Figura_Ontologias_Final.png", 
  plot = last_plot(),       
  width = 210,              
  height = 297,             
  units = "mm", 
  dpi = 300,                
  bg = "white"              
)
```
