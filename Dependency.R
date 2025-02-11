# devtools::install_github("caijun/ggcorrplot2")
# devtools::install_github("noriakis/CBNplot")
# install.packages("D:/Desktop/CBNplot-main/CBNplot-main", repos = NULL, type = "source")
# devtools::install_github("BioSenior/ggvolcano")
# BiocManager::install('ropls')
# BiocManager::install('clusterProfiler')
# BiocManager::install('org.Hs.eg.db')
# BiocManager::install("ReactomePA")
# BiocManager::install(c("STRINGdb","igraph"),ask = F,update = F)
# install.packages("Vennerable", repos="http://R-Forge.R-project.org")
# BiocManager::install('mixOmics')

# devtools::install_github("taowenmicro/ggClusterNet")
# devtools::install_github("taowenmicro/EasyMicrobiome")
# devtools::install_github("WangLab-MSSM/DreamAI/Code")

library('ggClusterNet')
library('mixOmics')
library('ggsci')
library('rjson')

library('openxlsx')
library('ggplot2')
library('dcurves')
library('ggbreak')
library('ggpubr')
library('tidyverse')
library('patchwork')
library('gghalves')
library('reshape2')
library('stringr')
library('ggrepel')
library('ggplotify')
library('nnet')
library('VennDiagram')
library('cowplot')
library('Vennerable')
library('GGally')
library('coin')
library('DreamAI')
# library('ggforce')
library('clinfun')
library('biomaRt')

library('ropls')
library('corrplot')
library('ComplexHeatmap')
library('ggcorrplot2')
library('jjAnno')
library("FactoMineR")
library("factoextra") 
library('pheatmap')
library('org.Hs.eg.db')
library('clusterProfiler')
# library('enrichplot')


library('rms')
library('ReactomePA')
library('Mfuzz')
library('cluster')
library('WGCNA')
library('ggExtra')
library('circlize')
library('ggVolcano')
# library('CBNplot')
library('ROCR')
library("dunn.test")

library('table1')
# dependencies
library('caret')
library('ordinal')
library('kernlab')
library('kknn')
library('ordinalForest')
library('brms')
library('glmpathcr')
library('rpartScore')
library('ordinalNet')

origin_path = 'E:/04_irAE_blood&protein_2024_07_31/04_irAE_blood&protein/'

#########install
Secreted_class = data.frame()
Secreted_fun = function(file, name){
  result = fromJSON(file = paste("sa_location_", file, ".json", sep = ''))
  temp <- as.data.frame(t(sapply(result, "[", i = 1:max(sapply(result, length)))))
  Secreted_class = rbind(Secreted_class, data.frame(Class = name, Gene = unlist(temp$Gene)))
  return(Secreted_class)
}
Secreted_class = Secreted_fun('Immunoglobulin', 'Immunoglobulin genes')
Secreted_class = Secreted_fun('Intracellular', 'Intracellular and membrane')
Secreted_class = Secreted_fun('Secreted_blood', 'Secreted to blood')
Secreted_class = Secreted_fun('Secreted_brain', 'Secreted in brain')
Secreted_class = Secreted_fun('Secreted_digestive_system', 'Secreted to digestive system')
Secreted_class = Secreted_fun('Secreted_extracellular matrix', 'Secreted to extracellular matrix')
Secreted_class = Secreted_fun('Secreted_female_reproductive_system', 'Secreted in female reproductive system')
Secreted_class = Secreted_fun('Secreted_male_reproductive_system', 'Secreted in male reproductive system')
Secreted_class = Secreted_fun('Secreted_other_tissues', 'Secreted in other tissues')
Secreted_class = Secreted_fun('Secreted_unknown_location', 'Secreted - unknown location')

Tissue_input = fromJSON(file = paste("tissue_category_rna_any_Tissue.json", sep = ''))
Tissue_input = as.data.frame(t(sapply(Tissue_input, "[", i = 1:max(sapply(Tissue_input, length)))))
Tissue_enriched = data.frame(Class = unlist(lapply(Tissue_input$`RNA tissue specific nTPM`, function(x) names(x)[1])),
                             Gene = unlist(Tissue_input$Gene))
