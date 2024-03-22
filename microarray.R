# 21-03-2024
# RMA normalization of microarray data
# setwd("C:/Users/salon/OneDrive/Desktop/DNA_Microarray")

# Load the libraries
library(tidyverse)
library(GEOquery)
library(BiocManager)
BiocManager::install("affy")
library(affy)

# Get supplementary files from NCBI GEO 
getGEOSuppFiles("GSE148537")

# Can be downloaded manually as well
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148537
# GSE148537_RAW.tar

# Uncompress the tar file
untar("GSE148537_RAW.tar")

# Reading in .CEL files
raw_data <- ReadAffy()
raw_data

# Perform RMA normlization
normalised_data <- rma(raw_data)
normalised_data

# Fetch the RMA normalised data from the object raw data (get expression estimates)
normalised_data_expression <- exprs(normalised_data)

# Comparison
boxplot(raw_data, main = "Raw Data", ylab = "Expression Value")
boxplot(normalised_data_expression, main = "Normalized Data", ylab = "Expression Value")

# Conversion of expression data matrix to data frame
df_normalised_data_expression <- data.frame(normalised_data_expression)
df_normalised_data_expression

# Map probe IDs to gene symbols
gse <- getGEO('GSE148537', GSEMatrix = T)

# Fetch feature data to get Id - gene symbol mapping
feature_data <- gse$GSE148537_series_matrix.txt.gz@featureData@data

# Subset of feature_data to obtain the columns containing probe ids and gene symbols
feature_data <- feature_data[,c(1,11)]

# Merging the data
df_normalised_data_expression <- df_normalised_data_expression %>% 
  rownames_to_column(var = 'ID') %>% 
  inner_join(., feature_data, by = 'ID')
