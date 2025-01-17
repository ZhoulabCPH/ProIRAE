########################################################################################
# Figure 3a
# WGCNA clustering 

Module_info = read.xlsx('./Output/Module_info.xlsx')
ME = moduleEigengenes(WGCNA_data, Module_info$Module)$eigengenes

biomarker = c('CRP', 'CD177', 'CD74', 'LCP1', 'ITM2B', intersect(DEP_list$C_vs_M$Gene, DEP_list$M_vs_S$Gene))

cor_data = cor(cbind(ME, WGCNA_data[, colnames(WGCNA_data) %in% biomarker]), method = 'pearson')
cor_value = cor_data[(length(unique(Module_info$Module))+1):(dim(cor_data)[1]), 1:length(unique(Module_info$Module))]

cor_multi_test = function(out, scoreMN, cal_name_1, cal_name_2){
  num_1 = length(unique(scoreMN[,colnames(scoreMN) == cal_name_1]))
  num_2 = length(unique(scoreMN[,colnames(scoreMN) == cal_name_2]))
  
  if (num_1<=3 & num_2<=3){
    p_value = data.frame(val_name_1 = cal_name_1,
                         val_name_2 = cal_name_2,
                         test_type = "Fisher' Test",
                         p_value = fisher.test(table(scoreMN[,colnames(scoreMN) %in% c(cal_name_1, cal_name_2)]))$p.value)
    out = rbind(out, p_value)
  }else{
    if (num_1 > 3 & num_2 > 3){
      p_value = data.frame(val_name_1 = cal_name_1,
                           val_name_2 = cal_name_2,
                           test_type = "Fisher' Test",
                           p_value = cor.test(scoreMN[,colnames(scoreMN) == cal_name_1], scoreMN[,colnames(scoreMN) == cal_name_2])$p.value)
      out = rbind(out, p_value)
    }else{
      V1 = cal_name_1
      V2 = cal_name_2
      if (num_1 > num_2){
        V1 = cal_name_2
        V2 = cal_name_1
      }
      scoreMN$temp1 = as.vector(scoreMN[,colnames(scoreMN) == V1])
      scoreMN$temp2 = as.vector(scoreMN[,colnames(scoreMN) == V2])
      if (min(num_1, num_2) == 2){
        p_value = data.frame(val_name_1 = cal_name_1,
                             val_name_2 = cal_name_2,
                             test_type = "Fisher' Test",
                             p_value = wilcox.test(temp2 ~ temp1, data = scoreMN)$p.value)
        out = rbind(out, p_value)
      }
      if (min(num_1, num_2) == 3){
        p_value = data.frame(val_name_1 = cal_name_1,
                             val_name_2 = cal_name_2,
                             test_type = "Fisher' Test",
                             p_value = kruskal.test(temp2 ~ temp1, data = scoreMN)$p.value)
        out = rbind(out, p_value)
      }
    }
  }
  return(out)
}

cor_out = data.frame()
for (name1 in colnames(ME)){
  for (name2 in biomarker[biomarker %in% rownames(cor_value)]){
    if (name1 != name2){
      print(paste(name1, name2))
      cor_out = cor_multi_test(cor_out, cor_data, name1, name2)
    }
  }
}

cor_out$p_adjust = p.adjust(cor_out$p_value, method = 'fdr')
cor_out = dcast(cor_out, val_name_2 ~ val_name_1, value.var = "p_adjust")
cor_p = cor_out[match(cor_out$val_name_2, rownames(cor_value)), paste('MEM', 0:(length(unique(Module_info$Module))-1), sep='')]
colnames(cor_p) = paste('M', 0:(length(unique(Module_info$Module))-1), sep='')

col_anno = unique(Module_info[,c('Module_col_remap', 'Module')])$Module_col_remap
names(col_anno) = unique(Module_info[,c('Module_col_remap', 'Module')])$Module

top_anno = HeatmapAnnotation(Module = unique(Module_info$Module),
                             col = list(Module = col_anno),
                             border = T,
                             show_annotation_name = F,
                             show_legend = F
)

colnames(cor_value) = str_split_fixed(colnames(cor_value), 'ME', 2)[,2]

cor_value = rbind(cor_value[rownames(cor_value) %in% c('CRP', 'CD177', 'CD74', 'LCP1', 'ITM2B'),],
                  cor_value[!rownames(cor_value) %in% c('CRP', 'CD177', 'CD74', 'LCP1', 'ITM2B'),])

p_cluster_1 = Heatmap(cor_value,
                      cluster_rows = F,
                      column_dend_height = unit(20, 'mm'),
                      col = colorRamp2(c(-1, 0, 1), c("#a698d6", '#f2e6e6', "#ffbe48")),
                      # col = colorRamp2(c(-1, 0, 1), c("#21a8a3", '#f2e6e6', "#e2001a")),
                      rect_gp = gpar(col = 'black', lwd = 1),
                      row_names_side = 'left',
                      column_names_side = 'top',
                      top_annotation = top_anno,
                      
                      cell_fun = function(j, i, x, y, width, height, fill){
                        if (cor_p[i,j] >= 0.05){
                          grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
                        }else{
                          if (cor_p[i,j] < 0.01){
                            grid.text(paste(sprintf("**"), round(cor_value[i, j], 3)), x, y, gp = gpar(fontsize = 10))
                          }else{
                            grid.text(paste(sprintf("*"), round(cor_value[i, j], 3)), x, y, gp = gpar(fontsize = 10))
                          }
                        }
                      },
                      heatmap_legend_param = list(title = "Pearson's r")
)






########################################################################################
# Figure 3b
ME = moduleEigengenes(WGCNA_data, Module_info$Module)$eigengenes
signedKME = signedKME(WGCNA_data, ME)

kME = as.data.frame(cor(WGCNA_data, ME, method = 'pearson', use = "p"))
kME$Gene = rownames(kME)
kME_data = reshape::melt(kME, id = 'Gene')
colnames(kME_data) = c('Gene', 'Module', 'cor')
kME_data$Module = str_split_fixed(kME_data$Module, 'ME', 2)[,2]

tSNE_gene = merge(Module_info[, c('Gene', 'Module', 'Module_col_remap')], kME_data, by = c('Gene', 'Module'))
tSNE_gene_q1 = tSNE_gene %>% group_by(Module, Module_col_remap) %>% summarise(abs_Q1 = quantile(abs(cor), 0.75))
tSNE_gene = merge(tSNE_gene, tSNE_gene_q1, by = c('Module', 'Module_col_remap'))
tSNE_gene = tSNE_gene[abs(tSNE_gene$cor) > tSNE_gene$abs_Q1, ]

library('Rtsne')

tsne_seed = function(i){
  set.seed(i)
  tSNE_out = cbind(tSNE_gene,
                   Rtsne(t(WGCNA_data[, tSNE_gene$Gene]),
                         dims = 2,
                         pca = F,
                         perplexity = 10,
                         beta = 0)$Y)
  colnames(tSNE_out) = c(colnames(tSNE_gene), 'Dim1', 'Dim2')
  
  p = ggplot()+
    geom_point(data = tSNE_out[tSNE_out$Module %in% c('M0', 'M1', 'M2', 'M4', 'M16'),], 
               aes(x = Dim1, y = Dim2, color = Module, size = abs(cor)), alpha = 0.9)+
    scale_color_manual(breaks = unique(Module_info[,c('Module', 'Module_col_remap')])$Module,
                       values = unique(Module_info[,c('Module', 'Module_col_remap')])$Module_col_remap)+
    theme_bw()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1, "cm"),
      axis.ticks.length.y.left = unit(0.1, "cm"),
      axis.text = element_text(color = 'black'),
      legend.position = 'none',
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.1, "inches")
    )+
    labs(x = 't-SNE Dim 1',
         y = 't-SNE Dim 2')
  
  return(p)
}
pdf('./Output/WGCNA_tSNE_irAE.pdf', height = 4, width = 4)
print(tsne_seed(7))
dev.off()








########################################################################################
# Figure 3c

go_annotation = function(Gene_list, name, out){
  print(length(Gene_list))
  EntrezID = bitr(Gene_list,fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb="org.Hs.eg.db")
  
  go = enrichGO(EntrezID$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont='ALL',
                pAdjustMethod = 'fdr',
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                keyType = 'ENTREZID')@result
  
  kegg = enrichKEGG(gene = EntrezID$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff =0.05,
                    pAdjustMethod = 'fdr',
                    qvalueCutoff = 0.05)@result
  
  reactome = enrichPathway(gene = EntrezID$ENTREZID,
                           organism='human',
                           pAdjustMethod = "fdr",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)@result
  
  if (dim(go)[1] == 0){
    return(rbind(out, 
                 cbind(name = name, ONTOLOGY = 'KEGG', kegg),
                 cbind(name = name, ONTOLOGY = 'Reactome', reactome)
    ))
  }
  else{
    return(rbind(out, 
                 cbind(name = name, go), 
                 cbind(name = name, ONTOLOGY = 'KEGG', kegg),
                 cbind(name = name, ONTOLOGY = 'Reactome', reactome)
    ))
  }
}

# Co-expression analysis (WGCNA)

Module_info = read.xlsx('./Output/Module_info.xlsx')
Module_count = as.data.frame(table(Module_info$Module))
colnames(Module_count) = c('Module', 'Module_count')

Trend_count = data.frame(Group = c('Strong_Up', 'Strong_Down', 'Weak_Up', 'Weak_Down', 'Up', 'Down'),
                         Group_count = c(184, 33, 85, 21, 269, 54))

Strong_Up = as.data.frame(table(Module_info[Module_info$Gene %in% intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), ]$Module))
Strong_Down = as.data.frame(table(Module_info[Module_info$Gene %in% intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), ]$Module))
Weak_Up = as.data.frame(table(Module_info[Module_info$Gene %in% setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), ]$Module))
Weak_Down = as.data.frame(table(Module_info[Module_info$Gene %in% setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), ]$Module))

WGCNA_Trend = full_join(full_join(full_join(Strong_Up, Strong_Down, by = 'Var1'), Weak_Up, by = 'Var1'), Weak_Down, by = 'Var1')
WGCNA_Trend[, 2:5] = replace(WGCNA_Trend[, 2:5], is.na(WGCNA_Trend[, 2:5]), 0)

colnames(WGCNA_Trend) = c('Module', 'Strong_Up', 'Strong_Down', 'Weak_Up', 'Weak_Down')

WGCNA_Trend = WGCNA_Trend[, 1:3]

WGCNA_Trend = reshape2::melt(WGCNA_Trend, id = 'Module')
colnames(WGCNA_Trend) = c('Module', 'Group', 'a')
WGCNA_Trend = left_join(left_join(WGCNA_Trend, Module_count, by = 'Module'), Trend_count, by = 'Group')

hyper = function(x){
  phyper(x[1]-1, x[3], 1459-x[3], x[2], lower.tail = FALSE)
}
WGCNA_Trend$p_value = apply(WGCNA_Trend[, 3:5], 1, hyper)
WGCNA_Trend[WGCNA_Trend$a !=0 , 'p_value_adj'] = p.adjust(WGCNA_Trend[WGCNA_Trend$a !=0 , ]$p_value, method = 'fdr')

p_value_adj = WGCNA_Trend[WGCNA_Trend$a !=0 , ]$p_value_adj
p_value_adj_str = ifelse(p_value_adj<0.0001, "****", 
                         ifelse(p_value_adj<0.001, "***",
                                ifelse(p_value_adj<0.01, "**", 
                                       ifelse(p_value_adj<0.05, "*", 
                                              ""))))
WGCNA_Trend[WGCNA_Trend$a !=0 , 'p_value_adj_str'] = p_value_adj_str

WGCNA_Trend$Module = factor(WGCNA_Trend$Module, levels = paste('M', 0:(length(unique(Module_info$Module))-1), sep = ''))
WGCNA_Trend$Group = factor(WGCNA_Trend$Group, levels = c('Strong_Up', 'Weak_Up','Strong_Down', 'Weak_Down'))

p_WGCNA = ggplot()+
  geom_bar(data = WGCNA_Trend[WGCNA_Trend$Group %in% c('Strong_Down', 'Strong_Up'),], 
           aes(x = Module, y = a, color = Group), width = 0.8, position = position_dodge(), stat = 'identity', fill = NA)+
  scale_color_manual(breaks = c('Weak_Down', 'Strong_Down', 'Weak_Up', 'Strong_Up'),
                     values = c('#c4e0ff', '#1f8fff', "#f8bdc6", '#e7233e'))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black')
    # axis.text.x = element_text(angle = 90, hjust = 1)
    # legend.justification = c(1,1), 
    # legend.position = c(1,1),
    # legend.text = element_text(size = 8),
    # legend.key.size = unit(0.1, "inches")
  )+
  labs(x = '', y = "No. of proteins")

pdf('./Output/WGCNA_DEPs_Trend_ratio.pdf', height =3, width = 8)
print(p_WGCNA)
dev.off()






########################################################################################
# Figure 3d
# WGCNA Trend rotein expressions heatmap cluster 2
annotation_col = Clinical_info[,c('Patient_stats', "BMI", 'Age', "Duration", "Sex")]

annotation_row = select_Module_info[, c('Gene', 'Module')]
annotation_row[annotation_row$Gene %in% Secreted_class$Gene, 'Secreted'] = "Secreted annotated"
annotation_row[annotation_row$Gene %in% Secreted_class[Secreted_class$Class == "Secreted to blood", ]$Gene, 'Secreted'] = "Secreted to blood"
annotation_row[annotation_row$Gene %in% Secreted_class[Secreted_class$Class == "Intracellular and membrane", ]$Gene, 'Secreted'] = "Intracellular and membrane"
annotation_row[annotation_row$Gene %in% Secreted_class[Secreted_class$Class == "Immunoglobulin genes", ]$Gene, 'Secreted'] = "Immunoglobulin genes"

annotation_row[annotation_row$Gene %in% Tissue_enriched$Gene, 'Tissue'] = "Tissue annotated"
annotation_row[annotation_row$Gene %in% Tissue_enriched[Tissue_enriched$Class == 'liver', ]$Gene, 'Tissue'] = "liver annotated"

rownames(annotation_row) = annotation_row$Gene

rownames(annotation_col) = Clinical_info$Id
annotation_colors = list(
  Patient_stats = c('Control'='#969696', 'Mild'= '#f6d377', 'Severe'='#e23d24')
)


data = protein
data$Gene = protein$Gene.name
data = merge(data, select_Module_info, by = 'Gene')

data$Module = factor(data$Module, levels = c('M1', 'M5', 'M7', 'M16', 'M2'))
data = data[order(data$Module), ]

matrix_out = as.matrix(data[,colnames(data) %in% Clinical_info[Clinical_info$Hospital_stats == 'Outpatient', ]$Id])
rownames(matrix_out) = data$Gene.name

z_score = function(x){
  (x - mean(x))/sd(x)
}
matrix_out = as.data.frame(t(apply(matrix_out, 1, z_score)))
matrix_out = as.data.frame(t(matrix_out))
matrix_out[matrix_out > 3] = 3
matrix_out[matrix_out < -3] = -3

pdf('./Output/WGCNA_DEPs_cluster_2.pdf', height = 8, width = 8)
add.flag(pheatmap(t(matrix_out),
                  # clustering_method = 'average',
                  annotation_col = annotation_col,
                  annotation_row = annotation_row,
                  annotation_colors = annotation_colors,
                  cluster_rows = F,
                  cluster_cols = T,
                  show_rownames = T,
                  annotation_legend = TRUE,
                  color = colorRampPalette(c('#2196f3', '#f2e6e6', '#f44336'))(n = 100)),         
         kept.labels = DEP_list$C_vs_M$Gene,
         repel.degree = 0.2
)
dev.off()





########################################################################################
# Figure 3e & Supplementary Figure 4

gene_ratio_str = c()
for (i in WGCNA_Trend$GeneRatio){
  gene_ratio_str = c(gene_ratio_str, eval(parse(text = i)))
}
WGCNA_Trend$GeneRatio_num = gene_ratio_str
WGCNA_Trend_backup = WGCNA_Trend

wgcna_profile_cluster = function(Ontology, limits = c(-12, 12)){
  
  WGCNA_Trend = WGCNA_Trend_backup[WGCNA_Trend_backup$ONTOLOGY %in% Ontology,]
  
  item = union(WGCNA_Trend[WGCNA_Trend$p.adjust < 0.05 & WGCNA_Trend$GeneRatio_num > 0.1, ]$Description,
               WGCNA_Trend[WGCNA_Trend$p.adjust < 0.05 & WGCNA_Trend$GeneRatio_num > 0.1, ]$Description)
  
  Data = WGCNA_Trend[WGCNA_Trend$p.adjust < 0.05 & WGCNA_Trend$Description %in% item, ]
  
  Data$name = factor(Data$name, levels =  rev(c('M1', 'M2', 'M16', 'M4')))
  Data = Data[order(Data$name), ]
  
  Data$Description = factor(Data$Description, levels = unique(rev(Data$Description)))
  Data$log10_p_adj = -log10(Data$p.adjust)
  limits = c(-12, 12)
  Data[Data$log10_p_adj > max(limits), 'log10_p_adj'] = max(limits)
  
  Data[Data$name %in% c('M4'), 'log10_p_adj'] = -Data[Data$name %in% c('M4'), ]$log10_p_adj
  
  p = ggplot()+
    geom_point(data = Data, aes(x=Description, y=name, fill = log10_p_adj, size = GeneRatio_num), shape = 22)+
    
    scale_fill_gradient2(low = '#1f8fff', mid = '#f2e6e6', high = '#e7233e', limits = limits)+
    scale_size(limits = c(0, 0.5), breaks = seq(0,0.5,0.1))+
    theme_bw()+
    ggtitle(Ontology)+
    theme(
      # panel.grid = element_blank(),
      axis.text = element_text(color = 'black')
      # axis.text.x = element_text(angle = 90, hjust = 1),
    )+
    labs(x = "", y= '')
  p
}

# annotation 1
WGCNA_Trend = data.frame()
WGCNA_Trend = go_annotation(select_Module_info[select_Module_info$Module == 'M1', ]$Gene, 'M1', WGCNA_Trend) # 65
WGCNA_Trend = go_annotation(select_Module_info[select_Module_info$Module == 'M4', ]$Gene, 'M4', WGCNA_Trend) # 18
gene_ratio_str = c()
for (i in WGCNA_Trend$GeneRatio){
  gene_ratio_str = c(gene_ratio_str, eval(parse(text = i)))
}
WGCNA_Trend$GeneRatio_num = gene_ratio_str
WGCNA_Trend_backup = WGCNA_Trend

p1 = wgcna_profile_cluster('KEGG')
p2 = wgcna_profile_cluster('Reactome')

pdf('./Output/WGCNA_DEPs_kegg_0.pdf', height =2, width = 8)
print(p1+p2+ plot_layout(widths = c(3,6.5), guides = "collect")+ theme(legend.position="right"))
dev.off()

pdf('./Output/WGCNA_DEPs_GO.pdf', height =5, width = 16)
print((p1+p3+p4 + plot_layout(widths = c(3,3,4)))/(p5) + plot_layout(guides = "collect")+ theme(legend.position="right"))
dev.off()