# Caroline Sullivan
# BF528
# Project 1 Part 5
# Role: Analyst 

library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

if (!require("amap", quietly = TRUE)){
  install.packages("amap")
}

# ----------------------- read in the data -------------------------------------
# reading in the normalized microarray from part 3, saved in a .csv file
# data is a sample x gene tibble

# now we can read data, providing our own header and skipping the original
data <- readr::read_csv("fully_filtered_data.csv", 
                        show_col_types = FALSE)

# in order to cluster patients, the data must be transposed/made tidy
data_t = pivot_longer(data, cols = -probesets, names_to = 'subject_id') %>%
  pivot_wider(names_from = probesets, values_from=value)

# ----------------------- hierarchical clustering ------------------------------
# 1. Perform hierarchical clustering on your fully filtered data matrix 
# from Part 4.4. Be sure to check that you are clustering the patients 
# and not the genes.

# from paper: - Clustering algorithm: hierarchical clustering
#             - Clustering metrics: (1-Pearson correlation) distance and Ward linkage

# first prep the data for the distance matrix calculation
data_t_df <- as.data.frame(data_t)
rownames(data_t_df) <- data_t_df$subject_id

# do distance calculation using 1 - Pearson correlation as specified
# with Dist function from amap
dist_mat = Dist(data_t_df[-1], method='correlation') 

# now complete hierarchical clustering
h_cluster = hclust(dist_mat, method='ward.D2')

# and cut
cluster_fit <- cutree(h_cluster, k = 2)
print(table(cluster_fit))

# 59 patients in cluster 1, 75 patients in cluster 2

# get metadata to figure out which patients correspond to which cluster
metadata <- readr::read_csv('proj_metadata.csv', show_col_types = FALSE) %>%
  dplyr::rename(cit_coloncancermolecularsubtype = 'cit-coloncancermolecularsubtype') %>%
  dplyr::mutate('color_bar' = if_else(cit_coloncancermolecularsubtype == 'C3', 'red', 'blue'))

# create a new matrix to set the color bar properly
cc <- metadata$color_bar

# now time to make a heatmap
data_matrix <- as.matrix(data_t_df[-1])

png("5_3heatmap_highres3.png", width = 8, height = 8, units = 'cm', res = 600)
heatmap(x = data_matrix, 
        RowSideColors = cc, 
        Rowv = as.dendrogram(h_cluster),
        Colv = NA,
        labCol=FALSE,
        labRow=FALSE)
dev.off()

# ----------------------- Welch t test -----------------------------------------

# import the 4.2 data to provide an analysis for the biologist
data4_2 <- readr::read_csv("4_2_filter_results.csv", 
                        show_col_types = FALSE)

data_4_2t = pivot_longer(data4_2, cols = -probesets, names_to = 'subject_id') %>%
  pivot_wider(names_from = probesets, values_from=value)

# first make a tibble with the cluster fit results
cluster_fit_df <- as.data.frame(cluster_fit) %>%
  as_tibble(rownames='subject_id') %>%
  dplyr::rename('subtype'=cluster_fit)

# add this data to data_t
data_with_subtypes <- dplyr::inner_join(data_t, 
                                        cluster_fit_df,
                                        by='subject_id')

# define new matrices for separate patient groups
patient_group1 <- data_with_subtypes %>% 
  dplyr::filter(subtype == 1) %>%
  dplyr::select(-'subtype')

patient_group2 <- data_with_subtypes %>% 
  dplyr::filter(subtype == 2) %>%
  dplyr::select(-'subtype')

# make dfs for t test
p1_df = as.data.frame(patient_group1)
p2_df = as.data.frame(patient_group2)
rownames(p1_df) <- p1_df$subject_id
rownames(p2_df) <- p2_df$subject_id

# get p values
p_val <- mapply(function(df1, df2) t.test(df1, df2)$p.value, 
                 df1=p1_df[-1], 
                 df2=p2_df[-1])%>%
  as_tibble(rownames='probesets') %>%
  dplyr::rename('p_val'=value)

# get t statistic
t_stat <- mapply(function(df1, df2) t.test(df1,df2)$statistic, 
                 df1=p1_df[-1], 
                 df2=p2_df[-1])%>%
  as_tibble(rownames='probesets') %>%
  dplyr::rename('t_stat'=value)

t_stat$probesets <- t_stat$probesets %>% stringr::str_sub(1,-3)
  
# combine the tibbles
data_stats <- dplyr::inner_join(p_val, t_stat, by='probesets') %>%
  dplyr::mutate(p_adj = p.adjust(p_val, method='fdr')) %>%
  dplyr::arrange(p_adj)

View(data_stats %>% dplyr::filter(p_adj<0.05))

# 1,251 gene probes meet the criteria using fully filtered data
# 23,781 gene probes meet the criteria using data from 4.5

# writing a file with the analysis results
readr::write_csv(data_stats, '5_6gene_stats.csv')
