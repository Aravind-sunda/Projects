# install packages from bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
BiocManager::install("hgu133plus2.db")

# load installed packages
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(tibble)

packageVersion('ggplot2')

# set working directory in scc
setwd("/projectnb/bf528/users/group_2/project-1-group_2/Samples/")


# read CEL files using ReadAffy() in working directory
read_CEL_files <- affy::ReadAffy()

# normalize CEL files
eset_rma <- affy::rma(read_CEL_files) # returns an expression set
head(eset_rma)


# first get output from readaffy()
output_of_ReadAffy <- fitPLM(read_CEL_files, normalize= TRUE, background = TRUE) 

# compute RLE and NUSE from samples
# NUSE and RLE resource -- chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.bioconductor.org/packages/devel/bioc/vignettes/affyPLM/inst/doc/QualityAssess.pdf
nuse_data <- NUSE(output_of_ReadAffy, type="stats") # compute NUSE
#head(nuse_data)
nuse_t <- as_tibble(t(nuse_data)) # or nuse_t <- as.data.frame(t(nuse_data))
#head(nuse_t)
#head(nuse_t$median)

rle_data <- RLE(output_of_ReadAffy, type="stats") # compute RLE
#head(rle_data)
rle_t <- as_tibble(t(rle_data))
rle_t
#head(rle_t$median)

#plot histogram and save jpeg in current directory
jpeg(file="/projectnb/bf528/users/group_2/project-1-group_2/Analysis/Programmer_RLE_histogram.jpeg")
hist(rle_t$median, xlab="median RLE values", main = "Distribution of median RLE values", col="lightblue")
dev.off()

jpeg(file="/projectnb/bf528/users/group_2/project-1-group_2/Analysis/Programmer_NUSE_histogram.jpeg")
hist(nuse_t$median, xlab="median NUSE values", main = "Distribution of median NUSE values", col="lightpink", xlim=c(0.97,1.08))
dev.off()

# Adjusting for batch effects 
# http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html#adjusting-for-batch-effects-with-combat
# https://rdrr.io/bioc/sva/man/ComBat.html AWESOME RESOURCE FOR COMBAT!!
proj_metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv") 
batch_effect <- proj_metadata$normalizationcombatbatch # batch effects
mod_foi = model.matrix(~as.factor(proj_metadata$normalizationcombatmod), data=proj_metadata) # Feature of Interest

expr_matrix <- exprs(eset_rma)
combat_data <- ComBat(dat=expr_matrix, batch=batch_effect, mod=mod_foi) 
head(combat_data)
# write files out to a csv
write.csv(combat_data, "/projectnb/bf528/users/group_2/project-1-group_2/Analysis/corrected_for_batch_effects.csv") 


# Perform PCA
transposed_combat_data <- t(combat_data)
scaled_combat_data_trans <- scale(transposed_combat_data)
scaled_combat_data <- t(scaled_combat_data_trans)
pca_data <- prcomp(scaled_combat_data, scale = FALSE, center = FALSE)
#head(pca_data)

# What is $rotation attribute
rotation <- pca_data$rotation
head(rotation)

# What is $importance attribute 
# https://stats.stackexchange.com/questions/254592/calculating-pca-variance-explained
# https://cmdlinetips.com/2019/04/introduction-to-pca-with-r-using-prcomp/
pca_sum <- summary(pca_data)
head(pca_sum)
pca_importance <- pca_sum$importance
head(pca_importance)
#head(pca_importance)[1]


# Plot PCA plot
library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
library(plotly)
pca_data$x

cancer_subtype <- proj_metadata$SixSubtypesClassification

jpeg(file="/projectnb/bf528/users/group_2/project-1-group_2/Analysis/Programmer_PCA-Plot.jpeg")
pca_data$rotation %>%
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(colour = cancer_subtype)) +
  #theme_bw() + 
  labs(x=paste0("PC1: ",round(pca_importance[2,1]*100,1),"%"),
       y=paste0("PC2: ",round(pca_importance[2,2]*100,1),"%")) +
  ggtitle("PCA Plot of PC1 vs PC2")
dev.off()