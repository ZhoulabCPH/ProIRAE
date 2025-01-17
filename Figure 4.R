
########################################################################################
# Figure 4a

pdf('./Output/DEPs_ELISA_compare.pdf', height = 3.2, width = 6)

for (name in intersect(DEP_list$C_vs_M$Gene, DEP_list$M_vs_S$Gene)){
  temp = ELISA_temp[, grepl(name, colnames(ELISA_temp)) | colnames(ELISA_temp) == 'Patient_stats']
  colnames(temp) = c('Original', 'Patient_stats', 'ELISA')
  
  p1 = ggplot()+
    geom_point(data = temp,
               aes(x = Original, y = ELISA, fill = Patient_stats, color= Patient_stats),
               shape=16, alpha = 0.7, size = 3)+
    scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'),
                       breaks = c('Control', 'Mild', 'Severe')) +
    scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'),
                      breaks = c('Control', 'Mild', 'Severe')) +
    scale_y_continuous(breaks = seq(0, 30, 1))+
    # theme_classic()+
    geom_smooth(data = temp,
                aes(x = Original, y = ELISA, size = 1.5),
                method = 'lm',se = T,size = 1.5, formula = y~x) +
    stat_cor(data = temp,
             aes(x = Original, y = ELISA, size = 1.5),
             method = "pearson",digits = 3,size=3)+
    ggtitle(name)+
    theme_bw()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1, "cm"),
      axis.ticks.length.y.left = unit(0.1, "cm"),
      axis.text = element_text(color = 'black'),
      legend.position = 'none',
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.1, "inches")
    )+
    labs(x = 'log2(MS)',
         y = 'log2(ELISA)')
  
  p2 = ggplot(temp, aes(x=Patient_stats, y= ELISA, color=Patient_stats, fill=Patient_stats)) +
    geom_half_boxplot(color='black', width = 0.7, position = position_dodge(0.8), alpha=1, outlier.shape = NA)+
    geom_half_point(position = position_dodge(0.7), shape = 21, alpha=0.6, size = 2, transformation = position_jitter(height = 0))+
    scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
    scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
    scale_y_continuous(breaks = seq(0, 30, 1))+
    annotate("text", x = 2, y = max(temp$ELISA, na.rm = T) *0.99, 
             label = paste('Kruskal Wallis\nP = ',  kruskal.test(ELISA~Patient_stats, temp)$p.value, sep=''))+
    geom_signif(comparisons = list(c('Control', 'Mild'), c('Mild', 'Severe'), c('Control', 'Severe')),
                map_signif_level = F, test = "wilcox.test", 
                y_position = c(0.83, 0.84, 0.85) * max(temp$ELISA, na.rm = T),
                tip_length = 0)+
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
    labs(x = '', y = 'log2(ELISA)')
  
  p3 = ggplot(temp, aes(x=Patient_stats, y= Original, color=Patient_stats, fill=Patient_stats)) +
    geom_half_boxplot(color='black', width = 0.7, position = position_dodge(0.8), alpha=1, outlier.shape = NA)+
    geom_half_point(position = position_dodge(0.7), shape = 21, alpha=0.6, size = 2, transformation = position_jitter(height = 0))+
    scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
    scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
    scale_y_continuous(breaks = seq(0, 30, 1))+
    annotate("text", x = 2, y = max(temp$Original, na.rm = T) *0.99, 
             label = paste('Kruskal Wallis\nP = ',  kruskal.test(Original~Patient_stats, temp)$p.value, sep=''))+
    geom_signif(comparisons = list(c('Control', 'Mild'), c('Mild', 'Severe'), c('Control', 'Severe')),
                map_signif_level = F, test = "wilcox.test", 
                y_position = c(0.83, 0.84, 0.85) * max(temp$Original, na.rm = T),
                tip_length = 0)+
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
    labs(x = '', y = 'log2(MS)')
  
  print(p1+p2+p3+plot_layout(widths = c(3, 1.5, 1.5)))
}
dev.off()


########################################################################################
# Figure 4b

feature_selection = function(selected_DIA){
  # VIP: Predict variable importance in projection
  vip_imp = as.data.frame(mixOmics::vip(mixOmics::plsda(selected_DIA, Clinical_info$Patient_stats[match(rownames(selected_DIA), Clinical_info$Id)]))[,1])
  colnames(vip_imp) = c('VIP')
  vip_imp$Protein = rownames(vip_imp)
  importance = vip_imp
  
  selected_DIA$Group = factor(str_split_fixed(rownames(selected_DIA), '-', 2)[,1], levels = c('C', 'M', 'S'))
  
  # Random Forest
  library(randomForest)
  set.seed(123)
  rf_model <- randomForest(Group ~ ., data = selected_DIA)
  rf_imp = as.data.frame(importance(rf_model))
  colnames(rf_imp) = c('RF')
  rf_imp$Protein = rownames(rf_imp)
  importance <- merge(importance, rf_imp, by = "Protein", all = TRUE)
  
  # XGBoost
  library(xgboost)
  set.seed(123)
  label <- as.numeric(selected_DIA$Group) - 1
  data_matrix <- xgb.DMatrix(data.matrix(selected_DIA[, -9]), label = label)
  model <- xgboost(data = data_matrix, max.depth = 3, eta = 0.1, nrounds = 10, objective = "multi:softmax", num_class = 3)
  xgb_imp <- xgb.importance(feature_names = colnames(selected_DIA[, -9]), model = model)
  names = 
    xgb_imp = data.frame(Protein = xgb_imp$Feature, XGBoost = xgb_imp$Gain)
  importance <- merge(importance, xgb_imp, by = "Protein", all = TRUE)
  
  # RFE: Recursive Feature Elimination
  library(caret)
  set.seed(123)
  control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
  rfe_model <- rfe(selected_DIA[, -9], selected_DIA$Group, sizes = c(1:3), rfeControl = control)
  print(rfe_model)
  # IL1RL1, FABP3, TNFAIP6
  
  # # Boruta
  # library(Boruta)
  # boruta_model <- Boruta(Group ~ ., data = selected_DIA, doTrace = 2)
  # print(boruta_model)
  # # No attributes deemed unimportant.
  return(importance)
}

importance_ELISA = feature_selection(na.omit(ELISA))
importance_MS = feature_selection(ELISA_DIA)

feature_selection_plot = function(model, hline){
  temp1 = importance_MS[c('Protein', model)]
  colnames(temp1) = c('Protein', 'model')
  p1 = ggplot(temp1,
              aes(x = factor(Protein, levels = temp1[order(temp1$model), ]$Protein), label = round(model, 3), y = model))+
    geom_hline(yintercept = hline, linetype ="dashed", color = 'darkgrey') +
    geom_bar(fill = '#f7d796', color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_text(color = 'black')+
    coord_flip()+
    ggtitle(paste(model,' for feature selection with'))+
    theme_bw()+
    theme_classic()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1, "cm"),
      axis.ticks.length.y.left = unit(0.1, "cm"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(color = 'black')
    )+
    labs(x = '', y = 'MS')
  
  temp2 = importance_ELISA[c('Protein', model)]
  colnames(temp2) = c('Protein', 'model')
  p2 = ggplot(temp2,
              aes(x = factor(Protein, levels = temp2[order(temp2$model), ]$Protein), label = round(model, 3), y = model))+
    geom_hline(yintercept = hline, linetype ="dashed", color = 'darkgrey') +
    geom_bar(fill = '#f7d796', color='black', stat  = 'identity', width = 0.7, position = "dodge")+
    geom_text(color = 'black')+
    coord_flip()+
    # ggtitle(model)+
    theme_bw()+
    theme_classic()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1, "cm"),
      axis.ticks.length.y.left = unit(0.1, "cm"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(color = 'black')
    )+
    labs(x = '', y = 'ELISA')
  return(p1+p2)
}

VIP_plot = feature_selection_plot('VIP', 1)
RF_plot = feature_selection_plot('RF', 3.5)
XGBoost_plot = feature_selection_plot('XGBoost', 0.1)

temp = data.frame(Protein = rep(c('IL1RL1', 'FABP3', 'PLIN3', 'VWCE', 'ADM', 'TNFAIP6', 'TNC', 'VSIG4'),2),
                  value = c(1,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0),
                  class = c(rep('MS', 8), rep('ELISA', 8)))
temp$Protein = factor(temp$Protein, levels = rev(c('IL1RL1', 'FABP3', 'PLIN3', 'VWCE', 'ADM', 'TNFAIP6', 'TNC', 'VSIG4')))

RFE_plot = ggplot()+
  geom_point(data = temp, aes(x=Protein, y=class, size = value), shape = 22, fill = '#f7d796')+
  coord_flip()+
  # ggtitle(model)+
  theme_bw()+
  theme_classic()+
  theme(
    axis.ticks.length.x.bottom = unit(0.1, "cm"),
    axis.ticks.length.y.left = unit(0.1, "cm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black')
  )+
  labs(x = '', y = '')

pdf('./Output/Performance_feature_selection_plot.pdf', height = 4, width =16)
print(VIP_plot | RF_plot | XGBoost_plot | RFE_plot + plot_layout(widths = c(3,3,3,2)))
dev.off()


########################################################################################
# Figure 4c-f

################################################################################################
# Model predicted value for discovery and external cohort

Gene_analysis = Gene_analysis = c('IL1RL1', 'FABP3')
family = "acat"
link = "cauchit"

set.seed(123)
model <- ordinalNet(x = as.matrix(input[, c(Gene_analysis)]), y = input[, ]$Group, family=family,link=link)
print(coef(model,matrix=TRUE))
result_value = as.data.frame(predict(model, newx = as.matrix(input[, Gene_analysis])))
rownames(result_value) = rownames(input)
colnames(result_value) = c('C', 'M', 'S')

result_value_val = as.data.frame(predict(model, newx = as.matrix(protein_external_Haerbin[, Gene_analysis])))
rownames(result_value_val) = rownames(protein_external_Haerbin)
colnames(result_value_val) = c('C', 'M', 'S')

result_value_3d = as.data.frame(predict(model, newx = as.matrix(input_3d[, Gene_analysis])))
rownames(result_value_3d) = rownames(input_3d)
colnames(result_value_3d) = c('C', 'M', 'S')


################################################################################################
# Model predicted value for discovery and external cohort multi+compare

Discovery_aupr = rbind(data_generate_aupr(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe'),
                       data_generate_aupr(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE'),
                       data_generate_aupr(prediction(input$IL1RL1, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe IL1RL1'),
                       data_generate_aupr(prediction(input$IL1RL1, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE IL1RL1'),
                       data_generate_aupr(prediction(input$FABP3, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe FABP3'),
                       data_generate_aupr(prediction(input$FABP3, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE FABP3')
)
Discovery_auc = rbind(data_generate_auc(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe'),
                      data_generate_auc(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE'),
                      data_generate_auc(prediction(input$IL1RL1, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe IL1RL1'),
                      data_generate_auc(prediction(input$IL1RL1, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE IL1RL1'),
                      data_generate_auc(prediction(input$FABP3, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe FABP3'),
                      data_generate_auc(prediction(input$FABP3, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE FABP3')
)

External_aupr = rbind(data_generate_aupr(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe'),
                      data_generate_aupr(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE'),
                      data_generate_aupr(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe IL1RL1'),
                      data_generate_aupr(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE IL1RL1'),
                      data_generate_aupr(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe FABP3'),
                      data_generate_aupr(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE FABP3')
)
External_auc = rbind(data_generate_auc(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe'),
                     data_generate_auc(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE'),
                     data_generate_auc(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe IL1RL1'),
                     data_generate_auc(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE IL1RL1'),
                     data_generate_auc(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe FABP3'),
                     data_generate_auc(prediction(protein_external_Haerbin$IL1RL1, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE FABP3')
)
pdf('./Output/Performance_model_performance_external_auc_aupr.pdf', height = 4, width =8)
print(auc_plot(Discovery_auc[grepl("Severe", Discovery_auc$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')) |
        aupr_plot(Discovery_aupr[grepl("Severe", Discovery_aupr$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')))
print(auc_plot(Discovery_auc[grepl("irAE", Discovery_auc$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')) |
        aupr_plot(Discovery_aupr[grepl("irAE", Discovery_aupr$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')))

print(auc_plot(External_auc[grepl("Severe", External_auc$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')) |
        aupr_plot(External_aupr[grepl("Severe", External_aupr$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')))
print(auc_plot(External_auc[grepl("irAE", External_auc$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')) |
        aupr_plot(External_aupr[grepl("irAE", External_aupr$group), ], c('#b81b25', '#eac4ab', '#cdd5d4')))
dev.off()

pdf('./Output/Patent_Performance_model_performance_external_auc_aupr.pdf', height = 4, width =8)
print(auc_plot(Discovery_auc[grepl("IL1RL1", Discovery_auc$group), ], c('#b81b25', '#f6d377')) |
        aupr_plot(Discovery_aupr[grepl("IL1RL1", Discovery_aupr$group), ], c('#b81b25',  '#f6d377')))
print(auc_plot(External_auc[grepl("IL1RL1", External_auc$group), ], c('#b81b25',  '#f6d377')) |
        aupr_plot(External_aupr[grepl("IL1RL1", External_aupr$group), ], c('#b81b25', '#f6d377')))

print(auc_plot(Discovery_auc[grepl("FABP3", Discovery_auc$group), ], c('#b81b25', '#f6d377')) |
        aupr_plot(Discovery_aupr[grepl("FABP3", Discovery_aupr$group), ], c('#b81b25',  '#f6d377')))
print(auc_plot(External_auc[grepl("FABP3", External_auc$group), ], c('#b81b25',  '#f6d377')) |
        aupr_plot(External_aupr[grepl("FABP3", External_aupr$group), ], c('#b81b25', '#f6d377')))

dev.off()




temp = protein_external_Haerbin[, c('IL1RL1', 'FABP3')]
temp$Patient_stats = substr(rownames(temp), 1,2)

p1 = ggplot(temp, aes(x=Patient_stats, y= IL1RL1, color=Patient_stats, fill=Patient_stats)) +
  geom_half_boxplot(color='black', width = 0.7, position = position_dodge(0.8), alpha=1, outlier.shape = NA)+
  geom_half_point(position = position_dodge(0.7), shape = 21, alpha=0.6, size = 2, transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  scale_y_continuous(breaks = seq(0, 30, 1))+
  annotate("text", x = 2, y = max(temp$IL1RL1, na.rm = T) *0.99, 
           label = paste('Kruskal Wallis\nP = ',  kruskal.test(IL1RL1~Patient_stats, temp)$p.value, sep=''))+
  geom_signif(comparisons = list(c('CC', 'CS'), c('CS', 'SS'), c('CC', 'SS')),
              map_signif_level = F, test = "wilcox.test", 
              y_position = c(0.83, 0.84, 0.85) * max(temp$IL1RL1, na.rm = T),
              tip_length = 0)+
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
  labs(x = '', y = 'log2(MS)')

p2 = ggplot(temp, aes(x=Patient_stats, y= FABP3, color=Patient_stats, fill=Patient_stats)) +
  geom_half_boxplot(color='black', width = 0.7, position = position_dodge(0.8), alpha=1, outlier.shape = NA)+
  geom_half_point(position = position_dodge(0.7), shape = 21, alpha=0.6, size = 2, transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  scale_color_manual(values = c('#969696', '#f6d377', '#e23d24'))+
  scale_y_continuous(breaks = seq(0, 30, 1))+
  annotate("text", x = 2, y = max(temp$IL1RL1, na.rm = T) *0.99, 
           label = paste('Kruskal Wallis\nP = ',  kruskal.test(FABP3~Patient_stats, temp)$p.value, sep=''))+
  geom_signif(comparisons = list(c('CC', 'CS'), c('CS', 'SS'), c('CC', 'SS')),
              map_signif_level = F, test = "wilcox.test", 
              y_position = c(0.83, 0.84, 0.85) * max(temp$FABP3, na.rm = T),
              tip_length = 0)+
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
  labs(x = '', y = 'log2(MS)')

pdf('./Output/Patent_validation_compare.pdf', height = 3.2, width = 3)
print(p1+p2)
dev.off()


# Model predicted value for discovery and external cohort

Discovery_aupr = rbind(data_generate_aupr(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe'),
                       data_generate_aupr(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE'))
Discovery_auc = rbind(data_generate_auc(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) == 'S'), 'Severe'),
                      data_generate_auc(prediction(result_value$S, unname(sapply(rownames(result_value), substr, 1, 1)) != 'C'), 'irAE'))

External_aupr = rbind(data_generate_aupr(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe'),
                      data_generate_aupr(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE'))
External_auc = rbind(data_generate_auc(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) == 'SS'), 'Severe'),
                     data_generate_auc(prediction(result_value_val$S, unname(sapply(rownames(result_value_val), substr, 1, 2)) != 'CC'), 'irAE'))


pdf('./Output/Performance_model_performance_external_auc_aupr_simple.pdf', height = 4, width =8)
print(auc_plot(Discovery_auc, c('#eac4ab', '#cdd5d4')) | aupr_plot(Discovery_aupr, c('#eac4ab', '#cdd5d4')))
print(auc_plot(External_auc, c('#eac4ab', '#cdd5d4')) | aupr_plot(External_aupr, c('#eac4ab', '#cdd5d4')))
dev.off()

# Predicted risk value density

colors <- colorRamp2(c(0, unique(Discovery_auc$cut_off)[2], unique(Discovery_auc$cut_off)[1], 1),
                     c('#798bc3', '#b0d8f5', 'white', "#bf5663"))


data = result_value[order(result_value$S), ]
data$Id = factor(rownames(data), levels = rownames(data))
data$group = unname(sapply(rownames(result_value), substr, 1, 1))

p1 = ggplot(data, aes(x=Id, y=S, fill=S), color = 'black')+
  geom_bar(stat = 'identity')+
  scale_y_continuous(breaks = c(0, unique(Discovery_auc$cut_off)[1], 1))+
  scale_fill_gradientn(colors = colors(seq(0, 1, length.out = 256))) +
  theme_classic()+
  theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        axis.ticks.length.y.left = unit(-0.1, "cm"),
        # axis.text.y = element_blank(),
        legend.position = "right")+
  labs(x=NULL, y='Predicted risk scores')

p2 = ggplot(data, aes(x = S, group = group, fill = group))+
  geom_density(alpha = 0.8)+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'),
                    breaks = c('C', 'M', 'S')) +
  theme_classic()+
  theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        axis.ticks.length.y.left = unit(-0.1, "cm"),
        # axis.text.y = element_blank(),
        legend.position = "right")+
  labs(x='Predicted risk scores', 'Density')

data = result_value_val[order(result_value_val$S), ]
data$Id = factor(rownames(data), levels = rownames(data))
data$group = unname(sapply(rownames(result_value_val), substr, 1, 2))

p3 = ggplot(data, aes(x=Id, y=S, fill=S), color = 'black')+
  geom_bar(stat = 'identity')+
  scale_fill_gradientn(colors = colors(seq(0, 1, length.out = 256))) +
  scale_y_continuous(breaks = c(0, unique(Discovery_auc$cut_off)[1], 1))+
  theme_classic()+
  theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        axis.ticks.length.y.left = unit(-0.1, "cm"),
        # axis.text.y = element_blank(),
        legend.position = "right")+
  labs(x=NULL, y='Predicted risk scores')

p4 = ggplot(data, aes(x = S, group = group, fill = group))+
  geom_density(alpha = 0.8)+
  scale_fill_manual(values = c('#969696', '#f6d377', '#e23d24'),
                    breaks = c('CC', 'CS', 'SS')) +
  theme_classic()+
  theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        axis.ticks.length.y.left = unit(-0.1, "cm"),
        # axis.text.y = element_blank(),
        legend.position = "right")+
  labs(x='Predicted risk scores', 'Density')

pdf('./Output/Performance_model_performance_density.pdf', height = 5, width = 10)
(p1/p2)|(p3/p4)
dev.off()

