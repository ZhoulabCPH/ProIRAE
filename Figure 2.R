########################################################################################
# Figure 2a & b

# Differential expression proteins: Kruskal Wills (Control vs. Mild vs. Severe)

data = as.data.frame(t(protein[,4:dim(protein)[2]]))
colnames(data) = protein$Gene.name
data$Id = rownames(data)
data = left_join(data, Clinical_info, by='Id')

temp = data[data$Patient_stats  %in% c('Control', 'Mild', 'Severe'),]

# 01: Jonckheere-Terpstra

p_value = c()
orientation = c()
for (name in protein$Gene.name){
  ready = data.frame(value = as.vector(temp[,name]),
                     group = c(rep(1, sum(temp$Patient_stats == 'Control')), 
                               rep(2, sum(temp$Patient_stats == 'Mild')), 
                               rep(3, sum(temp$Patient_stats == 'Severe'))))
  ready = ready[!is.na(ready$value),]
  
  if (median(ready[ready$group == 1, ]$value) < median(ready[ready$group == 3, ]$value)){
    orientation = c(orientation, 'Up')
  } else {
    orientation = c(orientation, 'Down')
  }
  
  set.seed(1)
  if (any(table(ready$group) > 1)){
    p_value = c(p_value, jonckheere.test(ready$value, ready$group, nperm=1000)$p.value)
  } else {
    p_value = c(p_value, jonckheere.test(ready$value, ready$group)$p.value)
  }
}

JT_test = data.frame(Gene = protein$Gene.name, 
                     p_value = p_value,
                     orientation = orientation,
                     Compare_name = 'Control vs. Mild vs. Severe (JT test)')
JT_test$p_adjust_fdr= p.adjust(JT_test$p_value, method = 'fdr')

# 02: Kruskal waills

p_value = c()
for (name in protein$Gene.name){
  compare = temp[!is.na(temp[, name]), c(name, 'Patient_stats')]
  
  if (length(compare[compare$Patient_stats  == 'Control', 1] != 0))
    p_value = c(p_value,
                kruskal.test(list(compare[compare$Patient_stats  == 'Control', 1], 
                                  compare[compare$Patient_stats  == 'Mild', 1], 
                                  compare[compare$Patient_stats  == 'Severe', 1]))$p.value
    )else{
      p_value = c(p_value, 1)
    }
}
KW_test = data.frame(Gene = protein$Gene.name, p_value = p_value, Compare_name = 'Control vs. Mild vs. Severe', log2FC = NA)
KW_test$p_adjust_fdr = p.adjust(KW_test$p_value, method = 'fdr')

# Control vs.Mild
temp = data[data$Patient_stats %in% c('Control', 'Mild', 'Severe'),]
Control = temp[temp$Patient_stats %in% c('Control'),]
Mild = temp[temp$Patient_stats %in% c('Mild'), ]
Severe = temp[temp$Patient_stats %in% c('Severe'),]

Wilcox = function(left, right, compare){
  p_value = c()
  for (gene in intersect(KW_test[KW_test$p_adjust_fdr<0.05,]$Gene, JT_test[JT_test$p_adjust_fdr<0.05,]$Gene)){
    p_value = c(p_value, wilcox.test(na.omit(left[, colnames(left) == gene]),
                                     na.omit(right[, colnames(right) == gene]))$p.value)
  }
  return(data.frame(Gene = intersect(KW_test[KW_test$p_adjust_fdr<0.05,]$Gene, JT_test[JT_test$p_adjust_fdr<0.05,]$Gene), 
                    p_value = p_value, compare = compare))
}
C1 = Wilcox(Control, Mild, 'Control vs. Mild')
C2 = Wilcox(Control, Severe, 'Control vs. Severe')
C3 = Wilcox(Mild, Severe, 'Mild vs. Severe')

Two_compare = rbind(C1, C2, C3)
Two_compare$p_adjust_fdr = p.adjust(Two_compare$p_value, method = 'fdr')

non_severe = Wilcox(rbind(Control, Mild), Severe, 'No severe vs. Severe')
non_irae = Wilcox(Control, rbind(Mild, Severe), 'No irAE vs. irAE')

irAE_compare = rbind(non_severe, non_irae)
irAE_compare$p_adjust_fdr = p.adjust(irAE_compare$p_value, method = 'fdr')

DEP_list = list(
  C_vs_M = Two_compare[Two_compare$compare == 'Control vs. Mild' & Two_compare$p_adjust_fdr < 0.05, ],
  C_vs_S = Two_compare[Two_compare$compare == 'Control vs. Severe' & Two_compare$p_adjust_fdr < 0.05, ],
  M_vs_S = Two_compare[Two_compare$compare == 'Mild vs. Severe' & Two_compare$p_adjust_fdr < 0.05, ],
  JT_test =JT_test[JT_test$p_adjust_fdr < 0.05,],
  KW_test = KW_test[KW_test$p_adjust_fdr < 0.05,],
  non_severe = irAE_compare[irAE_compare$compare == 'No severe vs. Severe' & irAE_compare$p_adjust_fdr < 0.05, ],
  non_irae = irAE_compare[irAE_compare$compare == 'No irAE vs. irAE' & irAE_compare$p_adjust_fdr < 0.05, ]
)

Sig_diff_Gene = intersect(KW_test[KW_test$p_adjust_fdr<0.05,]$Gene, JT_test[JT_test$p_adjust_fdr<0.05,]$Gene)

Two_compare_dcast = reshape2::dcast(Two_compare, Gene~compare, value.var = 'p_adjust_fdr')
irAE_compare_dcast = reshape2::dcast(irAE_compare, Gene~compare, value.var = 'p_adjust_fdr')
JT_test_dcast = JT_test[, c(1, 3, 5)]
KW_test_dcast = KW_test[, c(1, 5)]
colnames(JT_test_dcast) = c("Gene", 'Orientation', 'JT_test')
colnames(KW_test_dcast) = c("Gene", 'KW_test')

Test_value = full_join(full_join(JT_test_dcast[JT_test$p_adjust_fdr<0.05,], KW_test_dcast[KW_test$p_adjust_fdr<0.05,], by = 'Gene'),
                       full_join(Two_compare_dcast, irAE_compare_dcast, by = 'Gene'),
                       by = 'Gene')
# Stastical test venn
venn_plot_fun = function(List){
  print(plot(Venn(List), doWeights = TRUE))
  venn.plot = venn.diagram(
    # category.names = c('Up', 'Down', 'Sig'),
    x = List,
    col=pal_lancet()(length(List)),
    scaled = TRUE,
    # alpha=0.8,
    fill = "white",  # 描边颜色
    cex=0.8,
    cat.pos = 7,
    cat.dist = 0.05,
    cat.fontface=3,
    fontfamily="serif",
    filename = NULL,
    height = 2500,
    width = 2500
  )
  vb = as_grob(venn.plot)
  gb = as.ggplot(vb)
  print(gb)
}

pdf('./Output/DEPs_Venn_Trend.pdf', height = 8, width =8)

venn_plot_fun(list(KW_test = DEP_list$KW_test$Gene,
                   JT_test_up = DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene,
                   JT_test_down = DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene))
venn_plot_fun(list(C_vs_M = DEP_list$C_vs_M$Gene,
                   C_vs_S = DEP_list$C_vs_S$Gene,
                   M_vs_S = DEP_list$M_vs_S$Gene))
venn_plot_fun(list(non_severe = DEP_list$non_severe$Gene,
                   non_irae = DEP_list$non_irae$Gene))

venn_plot_fun(list(C_vs_M = DEP_list$C_vs_M$Gene,
                   M_vs_S = DEP_list$M_vs_S$Gene,
                   non_irae = DEP_list$non_irae$Gene))
dev.off()







########################################################################################
# Figure 2c
# clustering check

gene_set = JT_test[JT_test$p_adjust_fdr < 0.05, ]$Gene
print(length(gene_set))
annotation_row = JT_test[JT_test$p_adjust_fdr < 0.05, c('Gene', 'orientation')]
rownames(annotation_row) = annotation_row$Gene
annotation_row[annotation_row$Gene %in% Sig_diff_Gene & annotation_row$orientation == 'Up', 'Group'] = 'Strong Up'
annotation_row[annotation_row$Gene %in% Sig_diff_Gene & annotation_row$orientation == 'Down', 'Group'] = 'Strong Down'
annotation_row[!annotation_row$Gene %in% Sig_diff_Gene & annotation_row$orientation == 'Up', 'Group'] = 'Weak Up'
annotation_row[!annotation_row$Gene %in% Sig_diff_Gene & annotation_row$orientation == 'Down', 'Group'] = 'Weak Down'


annotation_row[annotation_row$Gene %in% DEP_list$C_vs_M$Gene, 'Fur_C_vs_M'] = 'Sig. in Control vs. Mild'
annotation_row[annotation_row$Gene %in% DEP_list$M_vs_S$Gene, 'Fur_M_vs_S'] = 'Sig. in Mild vs. Severe'
annotation_row[annotation_row$Gene %in% DEP_list$C_vs_S$Gene, 'Fur_C_vs_S'] = 'Sig. in Control vs. Severe'

annotation_row$Group = factor(annotation_row$Group, levels = c('Weak Down', 'Strong Down', 'Weak Up', 'Strong Up'))
annotation_row$Fur_C_vs_M = factor(annotation_row$Fur_C_vs_M, levels = c('Sig. in Control vs. Mild', 'NA'))
annotation_row[, c("Group", "Fur_C_vs_M", "Fur_M_vs_S", "Fur_C_vs_S")]

annotation_colors = list(
  Group = c('Weak Down'='#a9deef', 'Strong Down'='#22a8d6', 'Weak Up' = "#f59891", 'Strong Up'='#ed3224'),
  Fur_C_vs_M = c('Sig. in Control vs. Mild' = '#26b3a9', 'NA' = 'white'),
  Fur_M_vs_S = c('Sig. in Mild vs. Severe' = '#dcc431', 'NA' = 'white'),
  Fur_C_vs_S = c('Sig. in Control vs. Severe' = '#df0784', 'NA' = 'white')
)

Input = protein
rownames(Input) = protein$Gene.name
Input = Input[rownames(Input) %in% gene_set,
              colnames(Input) %in% Clinical_info[Clinical_info$Patient_stats %in% c('Control', 'Mild', 'Severe') , ]$Id]
Input = Input[order(factor(rownames(Input), levels = gene_set)), ]

z_score = function(x){
  (x - mean(x))/sd(x)
}
norm_data = as.data.frame(t(apply(Input, 1, z_score)))

norm_data = data.frame(Control = apply(norm_data[,colnames(norm_data) %in% Clinical_info[Clinical_info$Patient_stats == 'Control', ]$Id] , 1, median),
                       Mild = apply(norm_data[,colnames(norm_data) %in% Clinical_info[Clinical_info$Patient_stats == 'Mild', ]$Id] , 1, median),
                       Severe  = apply(norm_data[,colnames(norm_data) %in% Clinical_info[Clinical_info$Patient_stats == 'Severe', ]$Id] , 1, median))

norm_data = rbind(
  norm_data[rownames(norm_data) %in% rownames(annotation_row[annotation_row$Group == 'Weak Down', ]), ],
  norm_data[rownames(norm_data) %in% rownames(annotation_row[annotation_row$Group == 'Strong Down', ]), ],
  norm_data[rownames(norm_data) %in% rownames(annotation_row[annotation_row$Group == 'Strong Up', ]), ],
  norm_data[rownames(norm_data) %in% rownames(annotation_row[annotation_row$Group == 'Weak Up', ]), ])

pdf('./Output/DEPs_Trend_cluster.pdf', height =8, width = 6)
add.flag(pheatmap(as.matrix(norm_data),
                  clustering_distance_rows = "euclidean", 
                  clustering_distance_cols = "euclidean", 
                  clustering_method = 'average',
                  cluster_rows = T,
                  cluster_cols = F,
                  show_rownames = T,
                  # main = column_title,
                  annotation_row = annotation_row, 
                  annotation_colors = annotation_colors,
                  annotation_legend = TRUE,
                  color = c(colorRampPalette(c('#1a94c8', "#e5e9ec"))(n = abs(as.integer(min(norm_data)*50))),
                            colorRampPalette(c("#e5e9ec", '#ef5045'))(n = abs(as.integer(max(norm_data)*50))))),
         kept.labels = DEP_list$C_vs_M$Gene,
         repel.degree = 0.2
)
dev.off()







########################################################################################
# Figure 2e & f

# volcano plot log2FC
data = as.data.frame(t(protein[,4:dim(protein)[2]]))
colnames(data) = protein$Gene.name
data$Id = rownames(data)
data = left_join(data, Clinical_info, by='Id')

Control = data[data$Patient_stats %in% c('Control'), colnames(data) %in% protein$Gene.name]
Mild = data[data$Patient_stats %in% c('Mild'), colnames(data) %in% protein$Gene.name]
Severe = data[data$Patient_stats %in% c('Severe'), colnames(data) %in% protein$Gene.name]

DE_analysis = data.frame(Gene = colnames(Severe),
                         log2_S_M = log2(2^as.vector(apply(Severe, 2, median))/2^as.vector(apply(Mild, 2, median))),
                         log2_M_C = log2(2^as.vector(apply(Mild, 2, median))/2^as.vector(apply(Control, 2, median))),
                         log2_S_C = log2(2^as.vector(apply(Severe, 2, median))/2^as.vector(apply(Control, 2, median)))
)

DE_analysis = merge(DE_analysis, KW_test, by ='Gene')
DE_analysis[DE_analysis$Gene %in% setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), 'Trend'] = 'Weak-Up'
DE_analysis[DE_analysis$Gene %in% setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), 'Trend'] = 'Weak-Down'
DE_analysis[DE_analysis$Gene %in% intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), 'Trend'] = 'Strong-Up'
DE_analysis[DE_analysis$Gene %in% intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), 'Trend'] = 'Strong-Down'
DE_analysis[is.na(DE_analysis$Trend), 'Trend'] = 'Normal'
DE_analysis$Trend= factor(DE_analysis$Trend, levels = c('Weak-Up', 'Weak-Down', 'Strong-Up', 'Strong-Down', 'Normal'))

print(DE_analysis[DE_analysis$Gene %in% c("CRP", "CD74", "LCP1", "CD177", "ITM2B"),])

# biomarker = c('CSF3', 'CSF2', 'CX3CL1', 'FGF2', 'IFNA2', 'IL12A', 'IL1A', 'IL1B', 'IL1RA',  'IL2', 'IL13',  'CD74','GNAL', 'ADPGK', 
#               'LCP1', 'CRP','IL6','CXCL10','CX3CL1','CD177','CD74','ITM2B')

# Test_value[Test_value$Gene %in% c("CRP", "CD74", "LCP1", "CD177", "ITM2B"),]

p_log2FC_1 = ggplot()+
  geom_point(data = DE_analysis[DE_analysis$Trend %in% c('Weak-Up', 'Weak-Down', 'Normal'), ], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=16, alpha = 0.4, color = 'lightgrey')+
  geom_point(data = DE_analysis[DE_analysis$Trend %in% c('Strong-Up', 'Strong-Down'), ], 
             aes(y = log2_M_C, x = log2_S_M, color = Trend, fill = Trend, size = -log10(p_adjust_fdr)),
             shape=16, alpha = 0.4)+
  geom_point(data = DE_analysis[DE_analysis$Gene %in% c("CRP", "CD74", "LCP1", "CD177", "ITM2B"),], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=21, alpha = 1, color = '#06b41f', fill = NA)+
  geom_point(data = DE_analysis[DE_analysis$Gene %in% intersect(DEP_list$C_vs_M$Gene, DEP_list$M_vs_S$Gene),], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=21, alpha = 1, color = 'black', fill = NA)+
  
  scale_color_manual(values = c('#ed3224', '#22a8d6', '#ed3224', '#22a8d6', 'lightgrey'),
                     breaks = c('Strong-Up', 'Strong-Down', 'Weak-Up', 'Weak-Down', 'Normal'))+
  geom_abline(intercept = .5, slope = -1, color = '#333333', linetype = 'dashed')+
  geom_abline(intercept = -.5, slope = -1, color = '#333333', linetype = 'dashed')+
  geom_abline(intercept = 0, slope = 1, color = '#333333', linetype = 'dashed')+
  scale_x_continuous(limits = c(-2.53, 4), breaks = seq(-4, 30, 2))+
  scale_y_continuous(limits = c(-2.53, 4), breaks = seq(-4, 30, 2))+
  scale_size(limits = c(0, 3.5))+
  
  geom_text_repel(data = DE_analysis[DE_analysis$Gene %in% c("CRP", "CD74", "LCP1", "CD177", "ITM2B"),],
                  aes(y = log2_M_C, x = log2_S_M, label = Gene),size = 3, vjust = 1,hjust = -0.1, color= '#06b41f')+

  theme_bw()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    axis.text = element_text(color = 'black'),
    legend.position = 'right',
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.1, "inches")
  )+
  labs(x = 'log2(Severe vs. Mild)',
       y = 'log2(Mild vs. Control)')


p_log2FC_2 = ggplot()+
  geom_point(data = DE_analysis[!DE_analysis$Gene %in% DEP_list$C_vs_M$Gene, ], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=16, alpha = 0.6, color = 'lightgrey')+
  geom_point(data = DE_analysis[DE_analysis$Gene %in% DEP_list$C_vs_M$Gene, ], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=1, alpha = 0.7, color = '#ef9332', fill = '#ef9332')+
  geom_point(data = DE_analysis[DE_analysis$Gene %in% intersect(DEP_list$C_vs_M$Gene, DEP_list$M_vs_S$Gene), ], 
             aes(y = log2_M_C, x = log2_S_M, size = -log10(p_adjust_fdr)),
             shape=16, alpha = 0.7, color = '#ef9332', fill = '#ef9332')+
  geom_abline(intercept = .5, slope = -1, color = '#333333', linetype = 'dashed')+
  geom_abline(intercept = -.5, slope = -1, color = '#333333', linetype = 'dashed')+
  geom_abline(intercept = 0, slope = 1, color = '#333333', linetype = 'dashed')+
  scale_x_continuous(limits = c(-2.53, 4), breaks = seq(-4, 30, 2))+
  scale_y_continuous(limits = c(-2.53, 4), breaks = seq(-4, 30, 2))+
  
  geom_text_repel(data = DE_analysis[DE_analysis$Gene %in% intersect(DEP_list$C_vs_M$Gene, DEP_list$M_vs_S$Gene),],
                  aes(y = log2_M_C, x = log2_S_M, label = Gene),size = 3, vjust = 1,hjust = -0.1, color= 'black')+
  # theme_classic()+
  theme_bw()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    axis.text = element_text(color = 'black'),
    legend.position = 'right',
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.1, "inches")
  )+
  labs(x = 'log2(Severe vs. Mild)',
       y = 'log2(Mild vs. Control)')

pdf('./Output/DEPs_scatter_Trend.pdf', height = 5.3, width =14)
print(p_log2FC_1 + p_log2FC_2 )
dev.off()






########################################################################################
# Figure 2d
# Kegg_Go_for JT & KW

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

Up_Down_Trend = data.frame()
Up_Down_Trend = go_annotation(setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), 'Weak_Down_Trend', Up_Down_Trend) # 21
Up_Down_Trend = go_annotation(setdiff(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), 'Weak_Up_Trend', Up_Down_Trend) #85
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$KW_test$Gene), 'Strong_Down_Trend', Up_Down_Trend) #33
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$KW_test$Gene), 'Strong_Up_Trend', Up_Down_Trend) #184
Up_Down_Trend = go_annotation(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, 'Down_Trend', Up_Down_Trend) # 54 
Up_Down_Trend = go_annotation(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, 'Up_Trend', Up_Down_Trend) # 269 

Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$C_vs_M$Gene), 'M_vs_S_Up', Up_Down_Trend) # 17
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$M_vs_S$Gene), 'M_vs_S_Down', Up_Down_Trend) # 24
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$M_vs_S$Gene), 'M_vs_S_Up', Up_Down_Trend) # 159
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$C_vs_S$Gene), 'C_vs_S_Down', Up_Down_Trend) # 33
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$C_vs_S$Gene), 'C_vs_S_Up', Up_Down_Trend) # 184

Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$non_severe$Gene), 'non_severe_Up', Up_Down_Trend) # 184
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$non_severe$Gene), 'non_severe_Down', Up_Down_Trend) # 33
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Up', ]$Gene, DEP_list$non_irae$Gene), 'non_irae_Up', Up_Down_Trend) # 145
Up_Down_Trend = go_annotation(intersect(DEP_list$JT_test[DEP_list$JT_test$orientation == 'Down', ]$Gene, DEP_list$non_irae$Gene), 'non_irae_Down', Up_Down_Trend) # 27

gene_ratio_str = c()
for (i in Up_Down_Trend$GeneRatio){
  gene_ratio_str = c(gene_ratio_str, eval(parse(text = i)))
}

Up_Down_Trend$GeneRatio_num = gene_ratio_str

convert_entrez_to_symbol <- function(entrez_ids) {
  symbol_genes <- bitr(unlist(strsplit(entrez_ids, "/")), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  return(paste(symbol_genes, collapse = "/"))
}

Up_Down_Trend$geneID = sapply(Up_Down_Trend$geneID, convert_entrez_to_symbol)


Up_Down_Trend_backup = Up_Down_Trend
write.xlsx(Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05, c(1:4, 10:12, 7, 8)], 
           './Output/GO_KEGG_Reactome.xlsx', row.names = F)

Up_Down_Trend = Up_Down_Trend[Up_Down_Trend$ONTOLOGY %in% c('KEGG', 'Reactome'),]

# Strong_trend_infomation

Trend_anno = function(Input_left, Input_right, Title, ORI = NULL){
  Input = full_join(Input_left,Input_right, by = c('ONTOLOGY', 'ID', 'Description'))
  Input$ord1 = factor(Input$ONTOLOGY, levels = c("MF", "CC", "BP", 'Reactome', "KEGG"))
  
  Input[is.na(Input$p.adjust.x) & !is.na(Input$p.adjust.y), 'ord2'] = 1
  Input[!is.na(Input$p.adjust.x) & is.na(Input$p.adjust.y), 'ord2'] = 3
  Input[is.na(Input$ord2), 'ord2'] = 2
  Input$GeneRatio_num.x = replace(Input$GeneRatio_num.x, is.na(Input$GeneRatio_num.x), 0.00001)
  Input$GeneRatio_num.y = replace(Input$GeneRatio_num.y, is.na(Input$GeneRatio_num.y), 0.00001)
  Input$ord3 = Input$GeneRatio_num.x-Input$GeneRatio_num.y
  Input = Input[order(Input$ord1, Input$ord2, Input$ord3),]
  
  Input$Description = factor(Input$Description, levels = rev(unique(Input$Description)))
  Input = Input[order(Input$Description),]
  
  num = Input %>% 
    group_by(ONTOLOGY) %>%
    summarise(n=n())
  num$ONTOLOGY = factor(num$ONTOLOGY, levels = rev(c("MF", "CC", "BP", 'Reactome', "KEGG")))
  num = num[order(num$ONTOLOGY),]
  num$num = rownames(num)
  num_sum = c(0, cumsum(num$n))
  num$min = num_sum[1:sum(c("MF", "CC", "BP", 'Reactome', "KEGG") %in%  num$ONTOLOGY)]
  num$max = num_sum[-1]
  
  p_part = ggplot()+
    geom_bar(data = Input,
             aes(x = Description, y = -GeneRatio_num.x, fill = log10(Input$p.adjust.x)>0),
             color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_bar(data = Input,
             aes(x = Description, y = GeneRatio_num.y, fill = -log10(Input$p.adjust.y)>0),
             color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_vline(xintercept = unique(cumsum(num$n)) +0.5, linetype ="dashed", color = 'lightgrey') +
    geom_hline(yintercept = c(-0.1, 0.1), linetype ="dashed", color = 'darkgrey') +
    geom_rect(data = num, aes(xmin = min+0.5, 
                              xmax = max+0.5, ymin=max(Input$GeneRatio_num.y)+0.01, 
                              ymax=max(Input$GeneRatio_num.y)+0.07, color=ONTOLOGY), fill='white')+
    geom_text(data = num, aes(x = (min+max)/2 +0.5, y=0.205, label = ONTOLOGY), color='black')+
    scale_x_discrete(breaks = unique(Input$Description))+
    ggtitle(Title)+
    labs(x = '', y = 'Gene Ratio')+
    coord_flip() 
  
  theme_self_define = function(){
    theme_bw()+
      theme_classic()+
      theme(
        axis.ticks.length.x.bottom = unit(0.1, "cm"),
        axis.ticks.length.y.left = unit(0.1, "cm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black')
        # axis.text.y = element_text(angle = 90, hjust = 1)
        # legend.justification = c(1,1), 
        # legend.position = c(1,1),
        # legend.text = element_text(size = 8),
        # legend.key.size = unit(0.1, "inches")
      )
  }
  
  if (is.null(ORI)){
    p = p_part + 
      scale_color_manual(breaks = c("MF", "CC", "BP", 'Reactome', "KEGG"),
                         values = c('#fad9a6', '#bfdaa0', '#d7edf1', '#fad9a6', '#bfdaa0'))+
      scale_fill_manual(breaks = c(TRUE, FALSE), values = c('#e7233e', '#1f8fff'))+
      guides(fill = guide_legend(title.hjust = 0, title.theme = element_text(size = 8)))+
      theme_self_define()
    
  } else if (ORI == 'Up'){
    p = p_part +
      scale_color_manual(breaks = c("MF", "CC", "BP", 'Reactome', "KEGG"),
                         values = c('#fad9a6', '#bfdaa0', '#d7edf1', '#fad9a6', '#bfdaa0'))+
      scale_fill_manual(breaks = c(TRUE, FALSE), values = c('#e7233e', '#1f8fff'))+
      guides(fill = guide_legend(title.hjust = 0, title.theme = element_text(size = 8)))+
      theme_self_define()
  } else {
    p = p_part +
      scale_color_manual(breaks = c("MF", "CC", "BP", "KEGG"),
                         values = c('#fad9a6', '#bfdaa0', '#d7edf1', '#fad9a6', '#bfdaa0'))+
      scale_fill_manual(breaks = c(TRUE, FALSE), values = c('#1f8fff', '#1f8fff'))+
      theme_self_define()
  }
  return(p)
}
theme_self_define_2 = function(){
  theme_bw()+
    theme_classic()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1, "cm"),
      axis.ticks.length.y.left = unit(0.1, "cm"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(color = 'black'),
      axis.text.y = element_blank()
    )
}


profile_cluster = function(Ontology, limits = c(-12, 12), item = 0){
  Up_Down_Trend = Up_Down_Trend_backup[Up_Down_Trend_backup$ONTOLOGY %in% Ontology, ]
  if (item == 0){
    item = union(Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 & Up_Down_Trend$GeneRatio_num > 0.1 & Up_Down_Trend$name %in% c("Down_Trend"), ]$Description,
                 Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 & Up_Down_Trend$GeneRatio_num > 0.1 & Up_Down_Trend$name %in% c("Up_Trend"), ]$Description)
  } else{
    item = Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 & Up_Down_Trend$GeneRatio_num > 0.1, ]$Description
  }
  
  Input_left=Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 & Up_Down_Trend$Description %in% item & Up_Down_Trend$name %in% c("Down_Trend"), ]
  Input_right=Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 & Up_Down_Trend$Description %in% item & Up_Down_Trend$name %in% c("Up_Trend"), ]
  
  Input = full_join(Input_left,Input_right, by = c('ONTOLOGY', 'ID', 'Description'))
  Input$ord1 = factor(Input$ONTOLOGY, levels = c("MF", "CC", "BP", 'Reactome', "KEGG"))
  
  Input[is.na(Input$p.adjust.x) & !is.na(Input$p.adjust.y), 'ord2'] = 1
  Input[!is.na(Input$p.adjust.x) & is.na(Input$p.adjust.y), 'ord2'] = 3
  Input[is.na(Input$ord2), 'ord2'] = 2
  Input$GeneRatio_num.x = replace(Input$GeneRatio_num.x, is.na(Input$GeneRatio_num.x), 0.00001)
  Input$GeneRatio_num.y = replace(Input$GeneRatio_num.y, is.na(Input$GeneRatio_num.y), 0.00001)
  Input$ord3 = Input$GeneRatio_num.x-Input$GeneRatio_num.y
  Input = Input[order(Input$ord1, Input$ord2, Input$ord3),]
  
  Input$Description = factor(Input$Description, levels = rev(unique(Input$Description)))
  
  Input = Input[order(Input$Description),]
  
  num = Input %>%
    group_by(ONTOLOGY) %>%
    summarise(n=n())
  num$ONTOLOGY = factor(num$ONTOLOGY, levels = rev(c("MF", "CC", "BP", 'Reactome', "KEGG")))
  num = num[order(num$ONTOLOGY),]
  num$num = rownames(num)
  num_sum = c(0, cumsum(num$n))
  num$min = num_sum[1:sum(c("MF", "CC", "BP", 'Reactome', "KEGG") %in%  num$ONTOLOGY)]
  num$max = num_sum[-1]
  
  Input[-log10(Input$p.adjust.x) > max(limits) & !is.na(Input$p.adjust.x), 'p.adjust.x'] = 1/(10^max(limits))
  Input[-log10(Input$p.adjust.y) > max(limits) & !is.na(Input$p.adjust.y), 'p.adjust.y'] = 1/(10^max(limits))
  
  p_part = ggplot()+
    geom_bar(data = Input,
             aes(x = Description, y = -GeneRatio_num.x, fill = log10(Input$p.adjust.x)),
             color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_bar(data = Input,
             aes(x = Description, y = GeneRatio_num.y, fill = -log10(Input$p.adjust.y)),
             color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_vline(xintercept = unique(cumsum(num$n)) +0.5, linetype ="dashed", color = 'lightgrey') +
    geom_hline(yintercept = c(-0.1, 0.1), linetype ="dashed", color = 'darkgrey') +
    geom_rect(data = num, aes(xmin = min+0.5,
                              xmax = max+0.5, ymin=max(Input$GeneRatio_num.y)+0.01,
                              ymax=max(Input$GeneRatio_num.y)+0.07, color=ONTOLOGY), fill='white')+
    geom_text(data = num, aes(x = (min+max)/2 +0.5, y=0.205, label = ONTOLOGY), color='black')+
    scale_x_discrete(breaks = unique(Input$Description))+
    # ggtitle(Title)+
    labs(x = '', y = 'Gene Ratio')
  # coord_flip()
  
  theme_self_define = function(){
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
      )
  }
  
  p_1 = p_part +
    scale_color_manual(breaks = c("MF", "CC", "BP", 'Reactome', "KEGG"),
                       values = c('#fad9a6', '#bfdaa0', '#d7edf1', '#fad9a6', '#bfdaa0'))+
    scale_fill_gradient2(low = '#1f8fff', mid = '#f2e6e6', high = '#e7233e', limits = limits)+
    # scale_fill_manual(breaks = c(TRUE, FALSE), values = c('#e7233e', '#1f8fff'))+
    # guides(fill = guide_legend(title.hjust = 0, title.theme = element_text(size = 8)))+
    theme_self_define()
  
  Data = Up_Down_Trend[Up_Down_Trend$p.adjust < 0.05 &
                         Up_Down_Trend$Description %in% item &
                         Up_Down_Trend$name %in% c("Weak_Down_Trend", "Strong_Down_Trend",
                                                   'M_vs_S_Down', 'non_irae_Down',
                                                   "Weak_Up_Trend", "Strong_Up_Trend",
                                                   'M_vs_S_Up', 'non_irae_Up'), ]
  Data$Description = factor(Data$Description, levels = levels(Input$Description))
  Data$name = factor(Data$name, levels =  c('M_vs_S_Down', 'non_irae_Down',
                                            "Strong_Down_Trend","Weak_Down_Trend",
                                            "Weak_Up_Trend", "Strong_Up_Trend",
                                            'non_irae_Up', 'M_vs_S_Up'))
  Data$log10_p_adj = -log10(Data$p.adjust)
  Data[Data$log10_p_adj > max(limits), 'log10_p_adj'] = max(limits)
  
  Data[Data$name %in% c("Weak_Down_Trend", "Strong_Down_Trend",
                        'M_vs_S_Down', 'non_irae_Down'), 'log10_p_adj'] =
    -Data[Data$name %in% c("Weak_Down_Trend", "Strong_Down_Trend",
                           'M_vs_S_Down', 'non_irae_Down'), ]$log10_p_adj
  
  p_2 = ggplot()+
    geom_point(data = Data, aes(x=Description, y=name, fill = log10_p_adj, size = GeneRatio_num), shape = 22)+
    
    scale_fill_gradient2(low = '#1f8fff', mid = '#f2e6e6', high = '#e7233e', limits = limits)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = 'black')
          # axis.text.x = element_text(angle = 90, hjust = 1)
    )+
    labs(x = "", y= '')
  return(p_1 / p_2)
}

pdf('./Output/DEPs_kegg_anno_trend.pdf', height = 5, width = 8)
print(profile_cluster(c('Reactome', 'KEGG')) +
        plot_layout(height = c(4 ,3), guides = "collect")+ theme(legend.position="right"))
dev.off()






########################################################################################
# Supplementary Figure 2

pdf('./Output/DEPs_GO_anno_trend.pdf', height = 5, width = 16)
print(profile_cluster(c("MF", "CC", "BP"), c(-10,10)) +
        plot_layout(height = c(4 ,3), guides = "collect")+ theme(legend.position="right"))
dev.off()
