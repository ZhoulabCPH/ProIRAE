
#########################################################################################################
# Table 1: clinical characteristic Discovery

Clinical_info[Clinical_info$Outcome %in% 'Death', 'Outcome'] = 'Death / Withdrawn'
Clinical_info[Clinical_info$Age >= 65 , 'Age_category'] = '65+'
Clinical_info[Clinical_info$Age < 65 , 'Age_category'] = '30-64'

Clinical_info[Clinical_info$BMI < 18.5 , 'BMI_category'] = 'Underweight (<18.5)'
Clinical_info[Clinical_info$BMI > 24.9 , 'BMI_category'] = 'Overweight (>24.9)'
Clinical_info[is.na(Clinical_info$BMI_category) , 'BMI_category'] = 'Healthyweight (18.5-24.9)'

Clinical_info[Clinical_info$Cancer_types %in% '01-Lung', 'Lung_cancer'] = 'With'
Clinical_info[!Clinical_info$Cancer_types %in% '01-Lung', 'Lung_cancer'] = 'Without'

Clinical_info$BMI_category = factor(Clinical_info$BMI_category, 
                                    levels = c('Underweight (<18.5)', 'Healthyweight (18.5-24.9)','Overweight (>24.9)'))
Clinical_info$Stage = factor(Clinical_info$Stage, levels =c('II', 'III', 'IV'))

label(Clinical_info$Sex) <- "Sex"
label(Clinical_info$Age) <- "Age (years)"
label(Clinical_info$Age_category) <- "Age (years)"
label(Clinical_info$BMI) <- "BMI (kg/m2)"
label(Clinical_info$BMI_category) <- "BMI (kg/m2)"
label(Clinical_info$Duration) <- "Duration (weeks)"
label(Clinical_info$Cancer_types) <- "Cancer types"
label(Clinical_info$Lung_cancer) <- "Lung cancer"
label(Clinical_info$Stage) <- "Stage"
label(Clinical_info$Outcome) <- "Outcome"
label(Clinical_info$Cardiovascular) <- "Cardiovascular grade"
label(Clinical_info$Liver) <- "Liver grade"
label(Clinical_info$Lung) <- "Lung grade"
label(Clinical_info$Skin) <- "Skin grade"
label(Clinical_info$diagnose) <- "Cancer"


Outpatient = Clinical_info[Clinical_info$Hospital_stats == 'Outpatient',]
Outpatient$Patient_stats = factor(Outpatient$Patient_stats, levels = c('Control', 'Mild', 'Severe'))

strata <- c(list(`Total` = Clinical_info[Clinical_info$Hospital_stats == 'Outpatient',]),
            split(Outpatient, Outpatient$Patient_stats))
)

labels <- list(variables = list(Sex = render.varlabel(Clinical_info$Sex), 
                                Age = render.varlabel(Clinical_info$Age),
                                Age_category = render.varlabel(Clinical_info$Age_category),
                                BMI = render.varlabel(Clinical_info$BMI),
                                BMI_category = render.varlabel(Clinical_info$BMI_category),
                                Duration = render.varlabel(Clinical_info$Duration),
                                Cancer_types = render.varlabel(Clinical_info$Cancer_types),
                                Lung_cancer = render.varlabel(Clinical_info$Lung_cancer),
                                Stage = render.varlabel(Clinical_info$Stage),
                                Outcome = render.varlabel(Clinical_info$Outcome),
                                Lung = render.varlabel(Clinical_info$Lung),
                                Cardiovascular = render.varlabel(Clinical_info$Cardiovascular),
                                Liver = render.varlabel(Clinical_info$Liver),
                                Skin = render.varlabel(Clinical_info$Skin),
                                diagnose = render.varlabel(Clinical_info$diagnose)
),
groups = list('', "Discovery")
)

caption = "Clinical characteristic of xxxxxxxxxx"  
footnote = "*All irAEs were classifed according to the United States Health and Human Services Common Terminology Criteria for Adverse Events (CTCAE) v.5.0 with grade ≥ 3 considered to be severe (High grade)"
table1(strata, 
       labels, 
       render.continuous = c(. = "Mean (SD)",
                             . = "Median [Min, Max]"),
       groupspan = c(1,3,1),
       caption = caption, 
       footnote = footnote)


#########################################################################################################
# Table 1: clinical characteristic validation

Clinical_info_external[Clinical_info_external$Age >= 65 , 'Age_category'] = '65+'
Clinical_info_external[Clinical_info_external$Age < 65 , 'Age_category'] = '30-64'

Clinical_info_external[Clinical_info_external$BMI < 18.5 , 'BMI_category'] = 'Underweight (<18.5)'
Clinical_info_external[Clinical_info_external$BMI > 24.9 , 'BMI_category'] = 'Overweight (>24.9)'
Clinical_info_external[is.na(Clinical_info_external$BMI_category) , 'BMI_category'] = 'Healthyweight (18.5-24.9)'


Clinical_info_external$BMI_category = factor(Clinical_info_external$BMI_category, 
                                             levels = c('Underweight (<18.5)', 'Healthyweight (18.5-24.9)','Overweight (>24.9)'))
Clinical_info_external$Stage = factor(Clinical_info_external$Stage, levels =c('I', 'II', 'III', 'IV'))

label(Clinical_info_external$Sex) <- "Sex"
label(Clinical_info_external$Age) <- "Age (years)"
label(Clinical_info_external$Age_category) <- "Age (years)"
label(Clinical_info_external$BMI) <- "BMI (kg/m2)"
label(Clinical_info_external$BMI_category) <- "BMI (kg/m2)"
label(Clinical_info_external$Duration) <- "Duration (weeks)"
label(Clinical_info_external$Stage) <- "Stage"
label(Clinical_info_external$Cancer) <- "Cancer"


Clinical_info_external$Patient_stats = factor(Clinical_info_external$Patient_stats, levels = c('Control', 'Mild', 'Severe'))

strata <- c(list(`Total` = Clinical_info_external),
            split(Clinical_info_external, Clinical_info_external$Patient_stats))

labels <- list(variables = list(Sex = render.varlabel(Clinical_info_external$Sex), 
                                Age = render.varlabel(Clinical_info_external$Age),
                                Age_category = render.varlabel(Clinical_info_external$Age_category),
                                BMI = render.varlabel(Clinical_info_external$BMI),
                                BMI_category = render.varlabel(Clinical_info_external$BMI_category),
                                Duration = render.varlabel(Clinical_info_external$Duration),
                                Stage = render.varlabel(Clinical_info_external$Stage),
                                Cancer = render.varlabel(Clinical_info_external$Cancer)
),
groups = list('', "External cohort")
)

caption = "Clinical characteristic of xxxxxxxxxx"  
footnote = "*All irAEs were classifed according to the United States Health and Human Services Common Terminology Criteria for Adverse Events (CTCAE) v.5.0 with grade ≥ 3 considered to be severe (High grade)"
table1(strata, 
       labels, 
       render.continuous = c(. = "Mean (SD)",
                             . = "Median [Min, Max]"),
       groupspan = c(1,3),
       caption = caption, 
       footnote = footnote)
