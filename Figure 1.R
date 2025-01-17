
########################################################################################
# Figure 1b
data = data.frame(x = c('Identified peptides', 'Identified proteins', 'Comparable proteins'),
                  y = c(15745, 2304, 1459))
data$x = factor(data$x, levels = c('Identified peptides', 'Identified proteins', 'Comparable proteins'))

p1 = ggplot()+
  geom_bar(data = data, aes(x = x, y = log10(y), fill = x), stat = 'identity')+
  scale_fill_manual(values = c('#8b8b8b', '#bebebe', '#d2b4de'))+
  geom_text_repel(data = data, aes(x = x, y = log10(y), label = y),  color = 'black')+
  geom_text_repel(data = data, aes(x = x, y = log10(y)/2, label = x),  color = 'black', angle = 90)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position = 'none',
    axis.text.x = element_blank()
  )+
  labs(x = '', y = 'Number')






########################################################################################
# Figure 1c

p_count = as.data.frame(colSums(!is.na(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Hospital_stats == 'Outpatient', ]$Id])))
colnames(p_count) = 'P_count'
p_count$Id = rownames(p_count)
p_count = merge(p_count, Clinical_info[,c('Id', 'Patient_stats')])

kruskal.test(P_count~Patient_stats, p_count)

p4 = ggplot(p_count, aes(x=Patient_stats, y=P_count, color=Patient_stats,fill=Patient_stats)) +
  geom_half_boxplot(color='black', width = 0.7, position = position_dodge(0.8), alpha=1, outlier.shape = NA)+
  geom_half_point(position = position_dodge(0.7), shape = 21, alpha=0.6, size = 2, transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
    # axis.text.x = element_blank()
  )+
  labs(x = '', y = 'Number of proteins')






########################################################################################
# Figure 1d

abundance = apply(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Control', ]$Id], 1, median, na.rm=T)
mean = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Control', ]$Id]), 1, mean, na.rm=T)
sd = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Control', ]$Id]), 1, sd, na.rm=T)
names(abundance) = protein_raw$Gene.name
abundance_1 = data.frame(value = as.vector(abundance), cv = as.vector(sd/mean*100),
                         Gene = names(abundance),
                         rank = as.vector(rank(-abundance)),
                         group = 'Control')
print(sum(abundance_1[abundance_1$rank <= 10, ]$value)/sum(abundance_1$value, na.rm = T))
a_list_1 = abundance_1[abundance_1$rank <= 10, ]$Gene

abundance = apply(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Mild', ]$Id], 1, median, na.rm=T)
mean = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Mild', ]$Id]), 1, mean, na.rm=T)
sd = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Mild', ]$Id]), 1, sd, na.rm=T)
names(abundance) = protein_raw$Gene.name
abundance_2 = data.frame(value = as.vector(abundance), cv = as.vector(sd/mean*100),
                         Gene = names(abundance),
                         rank = as.vector(rank(-abundance)),
                         group = 'Mild')
print(sum(abundance_2[abundance_2$rank <= 10, ]$value)/sum(abundance_2$value, na.rm = T))
a_list_2 = abundance_2[abundance_2$rank <= 10, ]$Gene

abundance = apply(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Severe', ]$Id], 1, median, na.rm=T)
mean = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Severe', ]$Id]), 1, mean, na.rm=T)
sd = apply(log2(protein_raw[, colnames(protein_raw) %in% Clinical_info[Clinical_info$Patient_stats == 'Severe', ]$Id]), 1, sd, na.rm=T)
names(abundance) = protein_raw$Gene.name
abundance_3 = data.frame(value = as.vector(abundance), cv = as.vector(sd/mean*100),
                         Gene = names(abundance),
                         rank = as.vector(rank(-abundance)),
                         group = 'Severe')
print(sum(abundance_3[abundance_3$rank <= 10, ]$value)/sum(abundance_3$value, na.rm = T))
a_list_3 = abundance_3[abundance_3$rank <= 10, ]$Gene

abundance = rbind(abundance_1, abundance_2, abundance_3)


p2 = ggplot()+
  geom_point(data = abundance, aes(x = rank, y = log10(value), fill = group, color = group), alpha = 0.3, shape = 1)+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  # geom_text_repel(data = data, aes(x = x, y = log10(y), label = y),  color = 'black')+
  # geom_text_repel(data = data, aes(x = x, y = log10(y)/2, label = x),  color = 'black', angle = 90)+
  scale_y_continuous(breaks = seq(0, 10, 1))+
  scale_x_continuous(breaks = seq(0, 1600, 300))+
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
    # axis.text.x = element_blank()
  )+
  labs(x = 'Abundance rank', y = 'Log10(mean protein intensity)')







########################################################################################
# Figure 1e

print(abundance %>% group_by(group) %>% summarise(median = median(cv, na.rm = T), ratio_30 = sum(cv < 30, na.rm = T)))

p3 = ggplot()+
  geom_histogram(data = abundance, 
                 aes(x = cv, fill = group),
                 binwidth = 1, color = "black", alpha = 1)+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  geom_vline(xintercept = 30, linetype ="dashed", color = 'red') +
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position = 'none'
    # axis.text.x = element_blank()
  )+
  labs(x = 'Cofficient of variantion (%)', y = 'Number of proteins')





########################################################################################
# Figure 1f

data = t(protein[, colnames(protein) %in% Clinical_info[Clinical_info$Hospital_stats == 'Outpatient',]$Id])

res.pca = opls(x = data)
scoreMN = as.data.frame(res.pca@scoreMN)
scoreMN$Id = rownames(scoreMN)
scoreMN = left_join(scoreMN, Clinical_info, by='Id')

p_3_PCA = ggplot(scoreMN, aes(x = p1, y = p2, color = Patient_stats, fill = Patient_stats, shape = Patient_stats))+
  geom_point(alpha=0.8)+
  stat_ellipse()+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'),
                     breaks = c('Control', 'Mild', 'Severe')) +
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'),
                    breaks = c('Control', 'Mild', 'Severe')) +
  geom_vline(xintercept = 0, linetype ="dashed") +
  geom_hline(yintercept = 0, linetype ="dashed") +
  # theme(title = element_text(size = 15), text = element_text(size = 15)) +
  # geom_text_repel(aes(label = Gene),size = 3, vjust = 1,hjust = -0.1)+
  # ggtitle('Outpatient irAEs stats: PCA')+
  scale_x_continuous(breaks = seq(-60, 60, 20))+
  scale_y_continuous(breaks = seq(-60, 60, 20))+
  theme_bw()+
  guides(fill = guide_legend(title = 'Outpatient stats'))+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position="none"
    # legend.justification = c(1,1),
    # legend.position = c(1,1),
    # legend.text = element_text(size = 8),
    # legend.key.size = unit(0.1, "inches")
  )+
  labs(x=paste('PCA dimension 1 (', res.pca@modelDF$R2X[1]*100, '%)', sep=''),
       y=paste('PCA dimension 2 (', res.pca@modelDF$R2X[2]*100, '%)', sep=''))






########################################################################################
# Figure 1g

plsda.datatm <-unclass(mixOmics::plsda(data, Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)], ncomp = 2))

scoreMN = as.data.frame(plsda.datatm$variates$X)
scoreMN$SampleType = Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)]
scoreMN$Id = rownames(scoreMN)
scoreMN = left_join(scoreMN, Clinical_info, by='Id')
scoreMN$p1 = scoreMN$comp1
scoreMN$p2 = scoreMN$comp2

p_3_PLS_DA = ggplot(scoreMN, aes(x = p1, y = p2, color = Patient_stats, fill = Patient_stats, shape = Patient_stats))+
  geom_point(alpha=0.8)+
  stat_ellipse()+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'),
                     breaks = c('Control', 'Mild', 'Severe')) +
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'),
                    breaks = c('Control', 'Mild', 'Severe')) +
  geom_vline(xintercept = 0, linetype ="dashed") +
  geom_hline(yintercept = 0, linetype ="dashed") +
  scale_x_continuous(breaks = seq(-60, 60, 20))+
  scale_y_continuous(breaks = seq(-60, 60, 20))+
  # theme(title = element_text(size = 15), text = element_text(size = 15)) +
  # geom_text_repel(aes(label = Gene),size = 3, vjust = 1,hjust = -0.1)+
  # ggtitle('Outpatient irAEs stats: PLS-DA')+
  guides(fill = 'none')+
  theme_bw()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black'),
    legend.position="right"
    # legend.justification = c(1,1),
    # legend.position = c(1,1),
    # legend.text = element_text(size = 8),
    # legend.key.size = unit(0.1, "inches")
  )+
  labs(x=paste('PLS-DA dimension 1 (', round(plsda.datatm$prop_expl_var$X[1]*100, 1), '%)', sep=''),
       y=paste('PLS-DA dimension 2 (', round(plsda.datatm$prop_expl_var$X[2]*100, 1), '%)', sep=''))






########################################################################################
# Figure 1h

cor_multi_test = function(out, cal_name, val_type = 'C'){
  if (val_type == 'C'){
    p_value = data.frame(val_type = val_type,
                         val_name = cal_name,
                         test_type = 'Spearman',
                         P1_p_value = cor.test(scoreMN$p1, scoreMN[,colnames(scoreMN) == cal_name])$p.value,
                         P2_p_value = cor.test(scoreMN$p2, scoreMN[,colnames(scoreMN) == cal_name])$p.value)
    out = rbind(out, p_value)
  }else{
    scoreMN$temp = as.vector(scoreMN[,colnames(scoreMN) == cal_name])
    if (length(unique(scoreMN$temp)) == 3){
      p_value = data.frame(val_type = val_type,
                           val_name = cal_name,
                           test_type = 'Kruskal Wallis Test',
                           P1_p_value = kruskal.test(p1 ~ temp, data = scoreMN)$p.value,
                           P2_p_value = kruskal.test(p2 ~ temp, data = scoreMN)$p.value)
      out = rbind(out, p_value)
    }else{
      p_value = data.frame(val_type = val_type,
                           val_name = cal_name,
                           test_type = 'Two-side Wilcox Test',
                           P1_p_value = wilcox.test(p1 ~ temp, data = scoreMN)$p.value,
                           P2_p_value = wilcox.test(p2 ~ temp, data = scoreMN)$p.value)
      out = rbind(out, p_value)
    }
  }
}

data = t(protein[, colnames(protein) %in% Clinical_info[Clinical_info$Hospital_stats == 'Outpatient',]$Id])
y = class.ind(Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)])
y = y[, colSums(y) !=0]

res.pca = opls(x = data)
scoreMN = as.data.frame(res.pca@scoreMN)
scoreMN$Id = rownames(scoreMN)
scoreMN = left_join(scoreMN, Clinical_info, by='Id')

scoreMN$Stage = scoreMN$Stage == 'IV'
scoreMN$Lung = scoreMN$Lung == 'None'

cor_out = data.frame()
cor_out = cor_multi_test(cor_out, 'Age')
cor_out = cor_multi_test(cor_out, 'BMI')
cor_out = cor_multi_test(cor_out, 'Duration')

cor_out = cor_multi_test(cor_out, 'Patient_stats', 'D')
cor_out = cor_multi_test(cor_out, 'Sex', 'D')
cor_out = cor_multi_test(cor_out, 'Age_category', 'D')
cor_out = cor_multi_test(cor_out, 'BMI_category', 'D')
cor_out = cor_multi_test(cor_out, 'Lung_cancer', 'D')
cor_out = cor_multi_test(cor_out, 'Stage', 'D')
cor_out = cor_multi_test(cor_out, 'Cardiovascular', 'D')
cor_out = cor_multi_test(cor_out, 'Lung', 'D')
cor_out = cor_multi_test(cor_out, 'Liver', 'D')

cor_data = reshape2::melt(cor_out[, c(2, 4:5)], id='val_name')
colnames(cor_data) = c('val_name', 'PCA', 'PCA_p')
cor_data = merge(cor_data, cor_out[, c('val_name', 'test_type')])
cor_data$PCA_p_adjust = p.adjust(cor_data$PCA_p, method = 'fdr')

cor_data$val_name = factor(cor_data$val_name,
                           levels = c("Age", "BMI", "Duration",
                                      "Sex", "Age_category", "BMI_category", "Lung_cancer",
                                      'Stage', "Lung","Cardiovascular", "Liver", "Patient_stats"))

p_pca_overview_1 = ggplot(cor_data, aes(x=val_name, y=PCA, fill = -log10(PCA_p)))+
  geom_point(data = cor_data[cor_data$PCA_p >= 0.05,], shape = 4, size = 2)+
  geom_point(data = cor_data[cor_data$PCA_p < 0.05,], shape = 22, size = 3)+
  geom_point(data = cor_data[cor_data$PCA_p < 0.01,], shape = 22, size = 6)+
  scale_x_discrete(breaks = levels(cor_data$val_name),
                   labels = c("Age", "BMI", "Duration",
                              "Sex", "Age category", "BMI category", "Lung cancer stats",
                              'Stage', "Lung","Cardiovascular", "Liver", "Patient group"),
                   limits = levels(cor_data$val_name)
  )+
  scale_y_discrete(breaks = unique(cor_data$PCA),
                   labels = paste(paste('PC', 1:2, sep = ''), ' (', res.pca@modelDF$R2X[1:2] * 100, '%)', sep=''))+
  
  scale_fill_continuous(low = "#a0bec7", high = "#da3b46", limits=c(0,8))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "", y= '')


data = t(protein[, colnames(protein) %in% Clinical_info[Clinical_info$Hospital_stats == 'Outpatient',]$Id])
y = class.ind(Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)])
y = y[, colSums(y) !=0]

plsda.datatm <-unclass(mixOmics::plsda(data, Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)], ncomp = 2))

scoreMN = as.data.frame(plsda.datatm$variates$X)
scoreMN$SampleType = Clinical_info$Patient_stats[match(rownames(data), Clinical_info$Id)]
scoreMN$Id = rownames(scoreMN)
scoreMN = left_join(scoreMN, Clinical_info, by='Id')
scoreMN$p1 = scoreMN$comp1
scoreMN$p2 = scoreMN$comp2

scoreMN$Stage = scoreMN$Stage == 'IV'
scoreMN$Lung = scoreMN$Lung == 'None'

cor_out = data.frame()
cor_out = cor_multi_test(cor_out, 'Age')
cor_out = cor_multi_test(cor_out, 'BMI')
cor_out = cor_multi_test(cor_out, 'Duration')

cor_out = cor_multi_test(cor_out, 'Patient_stats', 'D')
cor_out = cor_multi_test(cor_out, 'Sex', 'D')
cor_out = cor_multi_test(cor_out, 'Age_category', 'D')
cor_out = cor_multi_test(cor_out, 'BMI_category', 'D')
cor_out = cor_multi_test(cor_out, 'Lung_cancer', 'D')
cor_out = cor_multi_test(cor_out, 'Stage', 'D')
cor_out = cor_multi_test(cor_out, 'Cardiovascular', 'D')
cor_out = cor_multi_test(cor_out, 'Lung', 'D')
cor_out = cor_multi_test(cor_out, 'Liver', 'D')

cor_data = reshape2::melt(cor_out[, c(2, 4:5)], id='val_name')
colnames(cor_data) = c('val_name', 'PCA', 'PCA_p')
cor_data = merge(cor_data, cor_out[, c('val_name', 'test_type')])
# cor_data$PCA_p_adjust = p.adjust(cor_data$PCA_p, method = 'fdr')

cor_data$val_name = factor(cor_data$val_name,
                           levels = c("Age", "BMI", "Duration",
                                      "Sex", "Age_category", "BMI_category", "Lung_cancer",
                                      'Stage', "Lung","Cardiovascular", "Liver", "Patient_stats"))

p_pca_overview_1_1 = ggplot(cor_data, aes(x=val_name, y=PCA, fill = -log10(PCA_p)))+
  geom_point(data = cor_data[cor_data$PCA_p >= 0.05,], shape = 4, size = 2)+
  geom_point(data = cor_data[cor_data$PCA_p < 0.05,], shape = 22, size = 3)+
  geom_point(data = cor_data[cor_data$PCA_p < 0.01,], shape = 22, size = 6)+
  scale_x_discrete(breaks = levels(cor_data$val_name),
                   labels = c("Age", "BMI", "Duration",
                              "Sex", "Age category", "BMI category", "Lung cancer stats",
                              'Stage', "Lung","Cardiovascular", "Liver", "Patient group"),
                   limits = levels(cor_data$val_name)
  )+
  scale_y_discrete(breaks = unique(cor_data$PCA),
                   labels = paste(paste('PC', 1:2, sep = ''), ' (',   round(plsda.datatm$prop_expl_var$X[1:2]*100, 1), '%)', sep=''))+
  
  scale_fill_continuous(low = "#a0bec7", high = "#da3b46", limits=c(0,8))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "", y= '')

#########

cor_multi_test_2 = function(out, cal_name_1, cal_name_2){
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
for (name1 in unique(cor_data$val_name)){
  for (name2 in unique(cor_data$val_name)){
    if (name1 != name2){
      print(paste(name1, name2))
      cor_out = cor_multi_test_2(cor_out, name1, name2)
    }
  }
}

cor_out$p_adjust = p.adjust(cor_out$p_value, method = 'fdr')

cor_out$val_name_1 = factor(cor_out$val_name_1,
                            levels = c("Age", "BMI", "Duration",
                                       "Sex", "Age_category", "BMI_category", "Lung_cancer",
                                       'Stage', "Lung","Cardiovascular", "Liver", "Patient_stats"))
cor_out$val_name_2 = factor(cor_out$val_name_2,
                            levels = c("Age", "BMI", "Duration",
                                       "Sex", "Age_category", "BMI_category", "Lung_cancer",
                                       'Stage', "Lung","Cardiovascular", "Liver", "Patient_stats"))

p_pca_overview_2 = ggplot()+
  geom_point(data = cor_out[cor_out$p_value >= 0.05,],
             aes(x=val_name_1, y=val_name_2, fill = -log10(p_value)),
             shape = 4, size = 2)+
  geom_point(data = cor_out[cor_out$p_value < 0.05,],
             aes(x=val_name_1, y=val_name_2, fill = -log10(p_value)),
             shape = 22, size = 3)+
  geom_point(data = cor_out[cor_out$p_value < 0.01,],
             aes(x=val_name_1, y=val_name_2, fill = -log10(p_value)),
             shape = 22, size = 6)+
  scale_x_discrete(breaks = levels(cor_out$val_name_1),
                   labels = c("Age", "BMI", "Duration",
                              "Sex", "Age category", "BMI category", "Lung cancer stats",
                              'Stage', "Lung","Cardiovascular", "Liver", "Patient group")
  )+
  scale_y_discrete(breaks = levels(cor_out$val_name_2),
                   labels = c("Age", "BMI", "Duration",
                              "Sex", "Age category", "BMI category", "Lung cancer stats",
                              'Stage', "Lung","Cardiovascular", "Liver", "Patient group")
  )+
  
  scale_fill_continuous(low = "#a0bec7", high = "#da3b46", limits=c(0,12))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "", y= '')








########################################################################################
# Figure 1i

data_HPA = Tissue_enriched %>% group_by(Class) %>% summarise(Count = n())
data_HMU = Tissue_enriched[Tissue_enriched$Gene %in% protein$Gene.name,] %>% group_by(Class) %>% summarise(Count = n())

data_HPA$Class = factor(data_HPA$Class, levels = c(data_HMU[order(-data_HMU$Count), ]$Class, setdiff(unique(data_HPA$Class), unique(data_HMU$Class))))
data_HMU$Class = factor(data_HMU$Class, levels = c(data_HMU[order(-data_HMU$Count), ]$Class, setdiff(unique(data_HPA$Class), unique(data_HMU$Class))))

Tissue_plot = ggplot()+
  geom_bar(data = data_HPA, aes(x = Class, y = Count), stat = 'identity', fill = NA, color = 'darkgrey')+
  geom_bar(data = data_HMU, aes(x = Class, y = Count), stat = 'identity', fill = '#8dd3c7', color = 'black')+
  
  geom_text_repel(data = data_HPA, aes(x = Class, y = Count+5, label = Count),  color = 'darkgrey')+
  geom_text_repel(data = data_HMU, aes(x = Class, y = Count+5, label = Count),  color = 'black')+
  
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 90),
    legend.position = 'none'
  )+
  labs(x = '',
       y = 'No. of proteins')

pdf('./Output/Overview_Tissue_class.pdf', height = 6, width = 16)
print(Tissue_plot)
List = list(HMU = protein$Gene.name,
            Tissue_enriched = Tissue_enriched$Gene)
print(plot(Venn(List), doWeights = TRUE))
dev.off()








########################################################################################
# Figure 1j

data_HPA = Secreted_class %>% group_by(Class) %>% summarise(Count = n())
data_HMU = Secreted_class[Secreted_class$Gene %in% protein$Gene.name,] %>% group_by(Class) %>% summarise(Count = n())

data_HPA$Class = factor(data_HPA$Class, levels = c(data_HMU[order(-data_HMU$Count), ]$Class, setdiff(unique(data_HPA$Class), unique(data_HMU$Class))))
data_HMU$Class = factor(data_HMU$Class, levels = c(data_HMU[order(-data_HMU$Count), ]$Class, setdiff(unique(data_HPA$Class), unique(data_HMU$Class))))

Secreted_plot = ggplot()+
  geom_bar(data = data_HPA, aes(x = Class, y = Count), stat = 'identity', fill = NA, color = 'darkgrey')+
  geom_bar(data = data_HMU, aes(x = Class, y = Count), stat = 'identity', fill = '#8dd3c7', color = 'black')+
  
  geom_text_repel(data = data_HPA, aes(x = Class, y = Count+5, label = Count),  color = 'darkgrey')+
  geom_text_repel(data = data_HMU, aes(x = Class, y = Count+5, label = Count),  color = 'black')+
  
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', angle = 90),
    legend.position = 'none'
  )+
  labs(x = '',
       y = 'No. of proteins')

pdf('./Output/Overview_Secreted_class.pdf', height = 6, width = 8)
print(Secreted_plot)
List = list(HMU = protein$Gene.name,
            HPA_Secreted = Secreted_class$Gene)
print(plot(Venn(List), doWeights = TRUE))
dev.off()
