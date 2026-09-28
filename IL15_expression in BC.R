#loading the data using tidyr
library(readr)
data <- read_csv("GSE183947_fpkm.csv.gz")
#load libraries
library(dbplyr)
library(tidyverse)
#############################################################################################
#get metadata
library(GEOquery)
gse <- getGEO(filename = "GSE183947_series_matrix.txt.gz", getGPL = FALSE)
meta_df <- pData(gse) #to get the phenodata
head(meta_df)
###########################################################################################
#select columns of interest in metadata
colnames(meta_df)[c(1, 10, 11, 17)] #select columns of interests
meta_df.modified <- meta_df%>%
  select(1,10,11,17)%>%
  rename(tissue=characteristics_ch1)%>%
  rename(metastasis=characteristics_ch1.1)%>%
  mutate(tissue=gsub("tissue: ", "", tissue))%>%
  mutate(metastasis=gsub("metastasis: ", "", metastasis))
head(data)
#########################################################################################
#reshaping the data(from wide, to long format)
data.long <- data%>%
  rename(gene=...1)%>%
  gather(key = "samples", value = "FPKM", -gene)
##############################################################################################
#Join the data frames, data.long and meta_df.modified
data.long <- data.long%>%
  left_join(., meta_df.modified, by = c("samples"="description"))
##########################################################################################
#Explore the data
data.long%>%
  filter(gene == "IL15"|gene == "IL15RA")%>%
  group_by(gene, tissue)%>%
  summarise(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM))%>%
  arrange(mean_FPKM)%>%
  head()
#####################################################################################################
###Visualizing the data using ggplot
#load libraries
library(tidyverse)
library(ggplot2)
#1. barplot
data.long%>%
  filter(gene == "IL15")%>%
  ggplot(., aes(x= samples, y= FPKM, fill = tissue))+
  geom_col()
#2. Density plot
data.long%>%
  filter(gene == "IL15")%>%
  ggplot(., aes(x= FPKM, fill = tissue))+
  geom_density(alpha = 0.3)
#3. violin plot
data.long%>%
  filter(gene == "IL15")%>%
  ggplot(., aes(x= metastasis, y= FPKM))+
  geom_violin()#can also use geom_violin
#4. Scatter plot to compare two genes
data.long%>%
  filter(gene == "IL15"|gene == "IL2")%>%
  spread(key = gene, value = FPKM)%>%
  ggplot(., aes(x = IL15, y = IL2, colour = tissue))+
  geom_point()+ geom_smooth(method = "lm", se = FALSE)
#5. Heatmap, to compare expression of multiple genes
genes.of.interest <- c(
  # Receptors
  "IL2RA", "IL2RB", "IL2RG", "IL15RA",
  # Cytokines
  "IL2", "IL15",
  # Signaling
  "JAK1", "JAK3", "STAT5A", "STAT5B", "STAT3",
  # Regulated / functional genes
  "FOXP3", "BCL2", "BCL2L1", "MCL1", "MYC", "CCND2", "CCND3",
  "PRF1", "GZMB", "IFNG", "TBX21", "EOMES",
  "SOCS1", "SOCS3", "CISH", "PRDM1"
)
heatmap <- data.long%>%
  filter(gene %in% genes.of.interest)%>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")
#6. Example on saving the plots
ggsave(heatmap, filename = "heatmap_save.pdf", width = 10, height = 8)
#############################################################################################################
