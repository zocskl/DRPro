#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# 蛋白——DR----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————

library(survival)
library(fdrtool)
library(readxl)
library(openxlsx)
library(dplyr)

subgroups <- c("main result"
               , "exlude CVD & CKD"
               , "2yrsExcluded"
               , "man"
               , "woman"
               , "PRS high"
               , "PRS low"
)

wb <- createWorkbook()
results_simple <- list()
for (subgroup in subgroups) {
  cox_result <- data.frame()
  
  for (protein in proteins) {
    formula <- paste0("Surv(followup_dr, incident_dr) ~ ", protein, " + a1 + as.factor(a2) + as.factor(a9)")
    
    if (subgroup == "man" | subgroup == "woman") {
      formula <- gsub("\\+ as\\.factor\\(a2\\)", "", formula)
    }
    formula <- as.formula(formula)
    
    cox_model <- coxph(formula, data = switch(subgroup,
                                              "main result" = data_pro
                                              , "man" = data_pro_man
                                              , "woman" = data_pro_woman
                                              , "PRS high" = data_pro_prs_high
                                              , "PRS low" = data_pro_prs_low
                                              , "2yrsExcluded" = data_pro_2yrsExcluded
    ))
    
    coef <- summary(cox_model)$coefficients[,1][1]
    HR <- summary(cox_model)$coefficients[,2][1]
    HR_low <- summary(cox_model)$conf.int[,3][1]
    HR_up <- summary(cox_model)$conf.int[,4][1]
    se <- summary(cox_model)$coefficients[,3][1]
    zvalue <- summary(cox_model)$coefficients[,4][1]
    pvalue <- summary(cox_model)$coefficients[,5][1]
    cox_result <- rbind(cox_result, c(protein, coef, HR, HR_low, HR_up, se, zvalue, pvalue))
  }
  
  colnames(cox_result) <- c("pro", "coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  columns_to_convert <- c("coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  cox_result[columns_to_convert] <- lapply(cox_result[columns_to_convert], as.numeric)
  cox_result$pvalue_sig <- ifelse(cox_result$pvalue <= 0.05, "significant", "-")
  
  cox_result$BH <- p.adjust(cox_result$pvalue, "BH")
  cox_result$BH_sig <- ifelse(cox_result$BH <= 0.05, "significant", "-")
  
  cox_result <- left_join(cox_result, proname, by = "pro")
  cox_result <- left_join(cox_result, panel, by = "pro")
  
  results_simple[[subgroup]] <- cox_result
  
  sheetname <- paste0("{model1} ", subgroup)
  addWorksheet(wb, sheetname)
  writeData(wb, sheetname, results_simple[[subgroup]])
}


results_full <- list()
for (subgroup in subgroups) {
  cox_result <- data.frame()
  
  for (protein in proteins) {
    formula <- paste0("Surv(followup_dr, incident_dr) ~", protein, " + a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm)")
    
    if (subgroup == "man" | subgroup == "woman") {
      formula <- gsub("\\+ as\\.factor\\(a2\\)", "", formula)
    }
    formula <- as.formula(formula)
    
    cox_model <- coxph(formula, data = switch(subgroup,
                                              "main result" = data_pro
                                              , "man" = data_pro_man
                                              , "woman" = data_pro_woman
                                              , "PRS high" = data_pro_prs_high
                                              , "PRS low" = data_pro_prs_low
                                              , "2yrsExcluded" = data_pro_2yrsExcluded
                                              ))
    
    coef <- summary(cox_model)$coefficients[,1][1]
    HR <- summary(cox_model)$coefficients[,2][1]
    HR_low <- summary(cox_model)$conf.int[,3][1]
    HR_up <- summary(cox_model)$conf.int[,4][1]
    se <- summary(cox_model)$coefficients[,3][1]
    zvalue <- summary(cox_model)$coefficients[,4][1]
    pvalue <- summary(cox_model)$coefficients[,5][1]
    
    cox_result <- rbind(cox_result, c(protein, coef, HR, HR_low, HR_up, se, zvalue, pvalue))
  }
  
  colnames(cox_result) <- c("pro", "coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  columns_to_convert <- c("coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  cox_result[columns_to_convert] <- lapply(cox_result[columns_to_convert], as.numeric)
  cox_result$pvalue_sig <- ifelse(cox_result$pvalue <= 0.05, "significant", "-")
  
  cox_result$BH <- p.adjust(cox_result$pvalue, "BH")
  cox_result$BH_sig <- ifelse(cox_result$BH <= 0.05, "significant", "-")
  
  cox_result <- left_join(cox_result, proname, by = "pro")
  cox_result <- left_join(cox_result, panel, by = "pro")
  
  results_full[[subgroup]] <- cox_result
  
  addWorksheet(wb, subgroup)
  writeData(wb, subgroup, results_full[[subgroup]])
}

covariates <- list(
  simple = "a1 + as.factor(a2) + as.factor(a9)",
  full = "a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm)"
)

wb <- createWorkbook()
result_med <- list()

for (covariate in names(covariates)) {
  proteins <- paste0("result_",covariate) %>% get() %>% pull(pro)
  
  cox_result <- data.frame()
    for (protein in proteins) {
    formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~", protein, " + ", covariates[[covariate]], " + as.factor(med_dm)"))
    cox_model <- coxph(formula, data = data_pro)
    
    coef <- summary(cox_model)$coefficients[,1][1]
    HR <- summary(cox_model)$coefficients[,2][1]
    HR_low <- summary(cox_model)$conf.int[,3][1]
    HR_up <- summary(cox_model)$conf.int[,4][1]
    se <- summary(cox_model)$coefficients[,3][1]
    zvalue <- summary(cox_model)$coefficients[,4][1]
    pvalue <- summary(cox_model)$coefficients[,5][1]
    
    cox_result <- rbind(cox_result, c(protein, coef, HR, HR_low, HR_up, se, zvalue, pvalue))
  }
  
  colnames(cox_result) <- c("pro", "coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  columns_to_convert <- c("coef", "HR", "HR_low", "HR_up", "se", "zvalue", "pvalue")
  cox_result[columns_to_convert] <- lapply(cox_result[columns_to_convert], as.numeric)
  cox_result$pvalue_sig <- ifelse(cox_result$pvalue <= 0.05, "significant", "-")
  
  cox_result$BH <- p.adjust(cox_result$pvalue, "BH")
  cox_result$BH_sig <- ifelse(cox_result$BH <= 0.05, "significant", "-")
  
  cox_result <- left_join(cox_result, proname, by = "pro")
  cox_result <- left_join(cox_result, panel, by = "pro")
  
  result_med[[covariate]] <- cox_result
  sheet_name <- paste0(covariate, " + med")
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, result_med[[covariate]])
}


r2_result <- data.frame()
for (pro in result_full$pro) {
  
  formula_base <- paste0(pro, " ~ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm)")
  model_base <- lm(formula_base, data = data_pro)
  
  formula_med <- paste0(pro, " ~ as.factor(med_dm)")
  model_med <- lm(formula_med, data = data_pro)
  
  formula_integrated <- paste0(formula_base, " + as.factor(med_dm)")
  model_integrated <- lm(formula_integrated, data = data_pro)
  
  r2_base = summary(model_base)$adj.r.squared
  r2_med = summary(model_med)$adj.r.squared
  r2_integrated = summary(model_integrated)$adj.r.squared
  r2_delta = r2_integrated - r2_base
  
  ci_base <- r2_ci(r2_base, n = nobs(model_base), k = length(coef(model_base)) - 1)
  ci_med <- r2_ci(r2_med, n = nobs(model_med), k = length(coef(model_med)) - 1)
  ci_integrated <- r2_ci(r2_integrated, n = nobs(model_integrated), k = length(coef(model_integrated)) - 1)
  
  se_delta <- sqrt((ci_base$se)^2 + (ci_integrated$se)^2)
  r2_delta_low <- r2_delta - 1.96 * se_delta
  r2_delta_up <- r2_delta + 1.96 * se_delta
  z_delta <- r2_delta / se_delta
  p_delta = 2 * (1 - pnorm(abs(z_delta)))
  
  r2_result <- rbind(r2_result,
                     data.frame(pro = pro,
                                r2_base = r2_base, r2_base_low = ci_base$ci_low, r2_base_up = ci_base$ci_up,
                                r2_med = r2_med, r2_med_low = ci_med$ci_low, r2_med_up = ci_med$ci_up,
                                r2_integrated = r2_integrated, r2_integrated_low = ci_integrated$ci_low, r2_integrated_up = ci_integrated$ci_up,
                                r2_delta = r2_delta, r2_delta_low = r2_delta_low, r2_delta_up = r2_delta_up,
                                p_delta = p_delta))
}

r2_result$BH_delta <- p.adjust(r2_result$p_delta, method = "BH")


wb <- createWorkbook()
# 基因
interaction <- data.frame()
for (pro in proteins) {
  
  formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm) + ", pro, " + as.factor(prs_category) + ", pro, " * as.factor(prs_category)"))
  cox_model <- coxph(formula, data = data_pro)
  
  pro_coef <- summary(cox_model)$coefficients[pro,"coef"]
  pro_HR <- summary(cox_model)$coefficients[pro,"exp(coef)"]
  pro_HR_low <- summary(cox_model)$conf.int[pro,"lower .95"]
  pro_HR_up <- summary(cox_model)$conf.int[pro,"upper .95"]
  pro_se <- summary(cox_model)$coefficients[pro,"se(coef)"]
  pro_pvalue <- summary(cox_model)$coefficients[pro,"Pr(>|z|)"]
  
  prs_coef <- summary(cox_model)$coefficients["as.factor(prs_category)1","coef"]
  prs_HR <- summary(cox_model)$coefficients["as.factor(prs_category)1","exp(coef)"]
  prs_HR_low <- summary(cox_model)$conf.int["as.factor(prs_category)1","lower .95"]
  prs_HR_up <- summary(cox_model)$conf.int["as.factor(prs_category)1","upper .95"]
  prs_se <- summary(cox_model)$coefficients["as.factor(prs_category)1","se(coef)"]
  prs_pvalue <- summary(cox_model)$coefficients["as.factor(prs_category)1","Pr(>|z|)"]
  
  interact_coef <- summary(cox_model)$coefficients[paste0(pro,":as.factor(prs_category)1"),"coef"]
  interact_HR <- summary(cox_model)$coefficients[paste0(pro,":as.factor(prs_category)1"),"exp(coef)"]
  interact_HR_low <- summary(cox_model)$conf.int[paste0(pro,":as.factor(prs_category)1"),"lower .95"]
  interact_HR_up <- summary(cox_model)$conf.int[paste0(pro,":as.factor(prs_category)1"),"upper .95"]
  interact_se <- summary(cox_model)$coefficients[paste0(pro,":as.factor(prs_category)1"),"se(coef)"]
  interact_pvalue <- summary(cox_model)$coefficients[paste0(pro,":as.factor(prs_category)1"),"Pr(>|z|)"]
  
  interaction <- rbind(interaction,
                       data.frame(pro = pro,
                                  pro_coef = pro_coef, pro_HR = pro_HR, pro_HR_low = pro_HR_low, pro_HR_up = pro_HR_up, pro_se = pro_se, pro_pvalue = pro_pvalue,
                                  prs_coef = prs_coef, prs_HR = prs_HR, prs_HR_low = prs_HR_low, prs_HR_up = prs_HR_up, prs_se = prs_se, prs_pvalue = prs_pvalue,
                                  interact_coef = interact_coef, interact_HR = interact_HR, interact_HR_low = interact_HR_low, interact_HR_up = interact_HR_up, interact_se = interact_se, interact_pvalue = interact_pvalue))
}

addWorksheet(wb, "prs-interaction")
writeData(wb, "prs-interaction", interaction)
# 性别
interaction <- data.frame()
for (pro in proteins) {
  
  formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm) + ", pro, " + as.factor(a2) + ", pro, " * as.factor(a2)"))
  cox_model <- coxph(formula, data = data_pro)
  
  pro_coef <- summary(cox_model)$coefficients[pro,"coef"]
  pro_HR <- summary(cox_model)$coefficients[pro,"exp(coef)"]
  pro_HR_low <- summary(cox_model)$conf.int[pro,"lower .95"]
  pro_HR_up <- summary(cox_model)$conf.int[pro,"upper .95"]
  pro_se <- summary(cox_model)$coefficients[pro,"se(coef)"]
  pro_pvalue <- summary(cox_model)$coefficients[pro,"Pr(>|z|)"]
  
  sex_coef <- summary(cox_model)$coefficients["as.factor(a2)1","coef"]
  sex_HR <- summary(cox_model)$coefficients["as.factor(a2)1","exp(coef)"]
  sex_HR_low <- summary(cox_model)$conf.int["as.factor(a2)1","lower .95"]
  sex_HR_up <- summary(cox_model)$conf.int["as.factor(a2)1","upper .95"]
  sex_se <- summary(cox_model)$coefficients["as.factor(a2)1","se(coef)"]
  sex_pvalue <- summary(cox_model)$coefficients["as.factor(a2)1","Pr(>|z|)"]
  
  interact_coef <- summary(cox_model)$coefficients[paste0(pro,":as.factor(a2)1"),"coef"]
  interact_HR <- summary(cox_model)$coefficients[paste0(pro,":as.factor(a2)1"),"exp(coef)"]
  interact_HR_low <- summary(cox_model)$conf.int[paste0(pro,":as.factor(a2)1"),"lower .95"]
  interact_HR_up <- summary(cox_model)$conf.int[paste0(pro,":as.factor(a2)1"),"upper .95"]
  interact_se <- summary(cox_model)$coefficients[paste0(pro,":as.factor(a2)1"),"se(coef)"]
  interact_pvalue <- summary(cox_model)$coefficients[paste0(pro,":as.factor(a2)1"),"Pr(>|z|)"]
  
  interaction <- rbind(interaction,
                       data.frame(pro = pro,
                                  pro_coef = pro_coef, pro_HR = pro_HR, pro_HR_low = pro_HR_low, pro_HR_up = pro_HR_up, pro_se = pro_se, pro_pvalue = pro_pvalue,
                                  sex_coef = sex_coef, sex_HR = sex_HR, sex_HR_low = sex_HR_low, sex_HR_up = sex_HR_up, sex_se = sex_se, sex_pvalue = sex_pvalue,
                                  interact_coef = interact_coef, interact_HR = interact_HR, interact_HR_low = interact_HR_low, interact_HR_up = interact_HR_up, interact_se = interact_se, interact_pvalue = interact_pvalue))
}


#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# 蛋白——基因----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————

lm_result <- data.frame()
for (protein in protein_sig) {
  formula <- paste0("prs ~", protein, "+ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm) + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
  lm_model <- lm(formula, data = data_pro)
  
  coef <- summary(lm_model)$coefficients[, 1][2]
  std <- summary(lm_model)$coefficients[, 2][2]
  low <- coef - 1.96*std
  high <- coef + 1.96*std
  tvalue <- summary(lm_model)$coefficients[, 3][2]
  pvalue <- summary(lm_model)$coefficients[, 4][2]
  
  lm_result <- rbind(lm_result, c(protein, coef, low, high, std, pvalue))
}

colnames(lm_result) <- c("pro", "coef", "low", "high", "std", "pvalue")
columns_to_convert <- c("coef", "low", "high", "std", "pvalue")
lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
lm_result$pvalue_sig <- ifelse(lm_result$pvalue <= 0.05, "significant", "-")

lm_result$BH <- p.adjust(lm_result$pvalue, "BH")
lm_result$BH_sig <- ifelse(lm_result$BH <= 0.05, "significant", "-")

lm_result <- left_join(lm_result, proname, by = "pro")
lm_result <- left_join(lm_result, panel, by = "pro")

lm_result <- data.frame()
for (protein in protein_sig) {
  formula <- paste0("prs_2 ~", protein, "+ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm) + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
  lm_model <- lm(formula, data = data_pro)
  
  coef <- summary(lm_model)$coefficients[, 1][2]
  std <- summary(lm_model)$coefficients[, 2][2]
  low <- coef - 1.96*std
  high <- coef + 1.96*std
  tvalue <- summary(lm_model)$coefficients[, 3][2]
  pvalue <- summary(lm_model)$coefficients[, 4][2]
  
  lm_result <- rbind(lm_result, c(protein, coef, low, high, std, pvalue))
}

colnames(lm_result) <- c("pro", "coef", "low", "high", "std", "pvalue")
columns_to_convert <- c("coef", "low", "high", "std", "pvalue")
lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
lm_result$pvalue_sig <- ifelse(lm_result$pvalue <= 0.05, "significant", "-")

lm_result$BH <- p.adjust(lm_result$pvalue, "BH")
lm_result$BH_sig <- ifelse(lm_result$BH <= 0.05, "significant", "-")

lm_result <- left_join(lm_result, proname, by = "pro")
lm_result <- left_join(lm_result, panel, by = "pro")


gene_predictors <- c('rs1801133_G',    
                     'rs12092121_G' ,   'rs2811893_T' ,   'rs6427247_A' ,   'rs4762_G'    ,   'rs699549_T'  ,   'rs763970_A'  ,   'rs1399634_T'   , 
                     'rs2380261_C'  ,  'rs1801282_C'  ,  'rs1197310_T'  ,  'rs4470583_A'  ,  'rs2910964_G'  ,  'rs1445754_T'  ,  'rs13163610_A' ,  'rs17376456_A'  , 
                     'rs2300782_C'  ,  'rs1800629_G'  ,  'rs184003_C'   ,  'rs2070600_C'  ,  'rs1800624_A'  ,  'rs1800625_A'  ,  'rs1224329_G'   , 
                     'rs1150790_A'  ,  'rs833061_C'   ,  'rs2010963_C'  ,  'rs2146323_C'  ,  'rs3025039_C'  ,  'rs713050_T'   ,  'rs487083_T'   ,  'rs4880_A'     ,  'rs39059_A'     , 
                     'rs759853_G'   ,  'rs2070744_C'  ,  'rs1571942_A'  ,  'rs12219125_G' ,  'rs4838605_C'  ,  'rs11101355_C'  , 
                     'rs11101357_G' ,  'rs4462262_T'  ,  'rs7903146_C'  ,  'rs899036_G'   ,  'rs10501943_T' ,
                     'rs9565164_T'  ,  'rs2031236_G'  ,  'rs7986566_C'  ,  'rs2038823_G'  ,  'rs3742872_G'  ,  'rs10519765_G' ,  'rs832882_G'   ,  'rs1024611_A'   , 
                     'rs599019_C'   ,  'rs13306430_G' ,  'rs5498_A'     ,  'rs1800469_A'  ,  'rs761207_T'   ,  'rs6031415_A')
snp <- snp[,c("eid_brain", "snp_sex", gene_predictors)]

lm_result <- data.frame()
for (protein in protein_prs_sig) {
  
  for (snp in gene_predictors) {
    formula <- paste0(protein, " ~ ", snp, "+ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm) + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    lm_model <- lm(formula, data = data_pro)
    
    coef <- summary(lm_model)$coefficients[snp,"Estimate"]
    se <- summary(lm_model)$coefficients[snp,"Std. Error"]
    low <- coef - 1.96*se
    high <- coef + 1.96*se
    pvalue <- 2 * pnorm(-abs(coef/se))
    
    lm_result <- rbind(lm_result, c(protein, snp, coef, low, high, pvalue))
    }
}

colnames(lm_result) <- c("pro", "snp", "coef", "low", "high", "pvalue")
columns_to_convert <- c("coef", "low", "high", "pvalue")
lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
lm_result$pvalue_sig <- ifelse(lm_result$pvalue <= 0.05, "significant", "-")
lm_result$BH <- p.adjust(lm_result$pvalue, "BH")
lm_result$BH_sig <- ifelse(lm_result$BH <= 0.05, "significant", "-")

lm_result <- left_join(lm_result, proname, by = "pro")
lm_result <- left_join(lm_result, panel, by = "pro")


#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# 蛋白——OCT----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————
oct_para <- c("overall", "rnfl", "gcipl", "isos", "rpe_layer")
wb <- createWorkbook()
lm_results_list <- list()
for (oct in oct_para) {
  
  lm_result <- data.frame()
  for (protein in protein_sig) {
    formula <- paste0(oct, " ~ ", protein, "+ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm)")
    lm_model <- lm(formula, data = data_pro_all)
    
    coef <- summary(lm_model)$coefficients[, 1][2]
    std <- summary(lm_model)$coefficients[, 2][2]
    low <- coef - 1.96*std
    high <- coef + 1.96*std
    tvalue <- summary(lm_model)$coefficients[, 3][2]
    pvalue <- summary(lm_model)$coefficients[, 4][2]
    
    lm_result <- rbind(lm_result, c(protein, coef, low, high, std, pvalue))
  }
  
  colnames(lm_result) <- c("pro", "coef", "low", "high", "std", "pvalue")
  columns_to_convert <- c("coef", "low", "high", "std", "pvalue")
  lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
  lm_result$pvalue_sig <- ifelse(lm_result$pvalue <= 0.05, "significant", "-")
  
  lm_result$BH <- p.adjust(lm_result$pvalue, "BH")
  lm_result$BH_sig <- ifelse(lm_result$BH <= 0.05, "significant", "-")
  
  lm_result <- left_join(lm_result, proname, by = "pro")
  lm_result <- left_join(lm_result, panel, by = "pro")
  lm_results_list[[oct]] <- lm_result
  addWorksheet(wb, oct)
  writeData(wb, oct, lm_results_list[[oct]])
}

# 亚组
data_list <- list(male = data_pro_all[data_pro_all$a2 == 1,],
                  female = data_pro_all[data_pro_all$a2 == 0,],
                  PRShigh = data_pro_all[data_pro_all$prs >= median(data_pro_all$prs, na.rm = T),],
                  PRSlow = data_pro_all[data_pro_all$prs < median(data_pro_all$prs, na.rm = T),])

for (subgroup in names(data_list)) {
  wb <- createWorkbook()
  lm_results_list <- list()
  
  for (oct in oct_para) {
    lm_result <- data.frame()
    
    for (protein in protein_sig) {
      formula <- paste0(oct, " ~ ", protein, "+ a1 + as.factor(a2) + as.factor(a9) + a18 + dm_duration + a10 + a16 + as.factor(predm)")
      lm_model <- lm(formula, data = data_list[[subgroup]])
      
      coef <- summary(lm_model)$coefficients[, 1][2]
      std <- summary(lm_model)$coefficients[, 2][2]
      low <- coef - 1.96*std
      high <- coef + 1.96*std
      tvalue <- summary(lm_model)$coefficients[, 3][2]
      pvalue <- summary(lm_model)$coefficients[, 4][2]
      
      lm_result <- rbind(lm_result, c(protein, coef, low, high, std, pvalue))
    }
    
    colnames(lm_result) <- c("pro", "coef", "low", "high", "std", "pvalue")
    columns_to_convert <- c("coef", "low", "high", "std", "pvalue")
    lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
    lm_result$pvalue_sig <- ifelse(lm_result$pvalue <= 0.05, "significant", "-")
    
    lm_result$BH <- p.adjust(lm_result$pvalue, "BH")
    lm_result$BH_sig <- ifelse(lm_result$BH <= 0.05, "significant", "-")
    
    lm_result <- left_join(lm_result, proname, by = "pro")
    lm_result <- left_join(lm_result, panel, by = "pro")
    
    lm_results_list[[oct]] <- lm_result
    addWorksheet(wb, oct)
    writeData(wb, oct, lm_results_list[[oct]])
  }
  
  path <- paste("F:/杨少鹏/000_2 文章/062 文章22-DR蛋白组/UKB_OCT——蛋白质_", subgroup, ".xlsx")
  saveWorkbook(wb, path, overwrite = TRUE)
}


#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# 蛋白——预测----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————
library(xgboost)
library(pc)
library(dplyr)
library(openxlsx)
library(haven)
library(cutoff)
library(caret)

outcome_list <- c("dr")

results_pro_all <- data.frame()
data_score <- data.frame(eid_ageing = data_pro$eid_ageing)
seeds <- c(678)

params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.1,
  max_depth = 3,
  min_child_weight = 3,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 2,
  alpha = 1
)

for (outcome in outcome_list) {
  test_results <- data.frame()
  data_train$dataset <- "train"
  data_test$dataset <- "test"
  
  data_train <- data_train[!is.na(data_train[[paste0("incident_", outcome)]]), ]
  data_test <- data_test[!is.na(data_test[[paste0("incident_", outcome)]]), ]
  
  train_y <- as.matrix(data_train[[paste0("incident_", outcome)]])
  test_y <- as.matrix(data_test[[paste0("incident_", outcome)]])
  
  train_x <- as.matrix(data_train[protein_predictor_all])
  test_x <- as.matrix(data_test[protein_predictor_all])
  train <- xgb.DMatrix(data = train_x, label = train_y)
  test <- xgb.DMatrix(data = test_x, label = test_y)
  
  val.cv <- xgb.cv(params = params, data = train, nfold = 20, nrounds = 1000, early_stopping_rounds = 20, metrics = "auc", verbose_eval = T)
  model <- xgboost(data = train, params = params, nrounds = val.cv$best_ntreelimit, verbose = 2)
  
  # 查看训练集
  pred_train <- predict(model, train)
  # 测试集结果
  pred <- predict(model, test)
  
  data_score_temp <- data.frame(eid_ageing = data_pro$eid_ageing)
  pred_score <- as.data.frame(cbind(data.frame(eid_ageing = data_test$eid_ageing), pred))
  data_score_temp <- left_join(data_score_temp, pred_score, by = "eid_ageing")
  train_score <- as.data.frame(cbind(data.frame(eid_ageing = data_train$eid_ageing), pred_train))
  data_score_temp <- left_join(data_score_temp, train_score, by = "eid_ageing")
  
  data_score_temp$dataset <- ifelse(!is.na(data_score_temp$pred), "test", ifelse(!is.na(data_score_temp$pred_train), "train", NA))
  data_score_temp$pred <- ifelse(!is.na(data_score_temp$pred), data_score_temp$pred, ifelse(!is.na(data_score_temp$pred_train), data_score_temp$pred_train, NA))
  data_score_temp$pred_train <- NULL
}

library(SHAPforxgboost)
shap_data <- shap.prep(model, X_train = train_x)

shap_summary <- shap_data %>%
  group_by(variable) %>%
  summarise(AverageMeanValue = mean(mean_value, na.rm = TRUE))

# 单个蛋白预测
c_results <- data.frame()
for (protein in proteins) {
  formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~", protein))
  cox_model <- coxph(formula, data = data_train)
  data_test$ypred <- predict(cox_model, type = "risk", newdata = data_test)
  
  c_obj <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred, data = data_test)
  c <- c_obj$concordance
  c_se <- sqrt(c_obj$var)
  c_low <- c - 1.96 * c_se
  c_up <- c + 1.96 * c_se
  
  c_results <- rbind(c_results, c(protein, c, c_low, c_up))
}
colnames(c_results) <- c("pro", "c", "c_low", "c_up")
cols_to_convert <- c("c", "c_low", "c_up")
c_results[cols_to_convert] <- lapply(c_results[cols_to_convert], function(x) as.numeric(as.character(x)))
c_results <- c_results[!duplicated(c_results$pro), ]


covariates_Aspelund <- c("as.factor(a2)","a16","dm_duration","a18")
covariates_Hippisley <- c("as.factor(a2)","a10","a16","TC_divide_HDL","a18")
covariates_Dagliati <- c("baselineage","as.factor(a2)","dm_duration","a10","a18","prev_hbp","as.factor(a6)")
covariates_all <- union(union(covariates_Aspelund, covariates_Hippisley), covariates_Dagliati)
covariates_all_formula <- paste(covariates_all, collapse = " + ")
formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ ", covariates_all_formula))
cox_model_conv <- coxph(formula, data = train_data)
test_data$ypred_conv <- predict(cox_model_conv, type = "risk", newdata = test_data)
train_data$ypred_conv <- predict(cox_model_conv, type = "risk", newdata = train_data)
c_conv <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv, data = test_data)

c_protein <- survConcordance(Surv(followup_dr, incident_dr) ~ score, data = test_data)

covariates_all_formula <- paste(covariates_all, collapse = " + ")
formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ score + ", covariates_all_formula))
cox_model_combined <- coxph(formula, data = train_data)
test_data$ypred_combined <- predict(cox_model_combined, type = "risk", newdata = test_data)
train_data$ypred_combined <- predict(cox_model_combined, type = "risk", newdata = train_data)
c_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_combined, data = test_data)

formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ prs"))
cox_model_prs <- coxph(formula, data = train_data)
test_data$ypred_prs <- predict(cox_model_prs, type = "risk", newdata = test_data)
train_data$ypred_prs <- predict(cox_model_prs, type = "risk", newdata = train_data)
c_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs, data = test_data)

formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ prs + score"))
cox_model_prs_combined <- coxph(formula, data = train_data)
test_data$ypred_prs_combined <- predict(cox_model_prs_combined, type = "risk", newdata = test_data)
train_data$ypred_prs_combined <- predict(cox_model_prs_combined, type = "risk", newdata = train_data)
c_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs_combined, data = test_data)

formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ prs + ", covariates_all_formula))
cox_model_conv_prs <- coxph(formula, data = train_data)
test_data$ypred_conv_prs <- predict(cox_model_conv_prs, type = "risk", newdata = test_data)
train_data$ypred_conv_prs <- predict(cox_model_conv_prs, type = "risk", newdata = train_data)
c_conv_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs, data = test_data)

formula <- as.formula(paste0("Surv(followup_dr, incident_dr) ~ score + prs + ", covariates_all_formula))
cox_model_conv_prs_combined <- coxph(formula, data = train_data)
test_data$ypred_conv_prs_combined <- predict(cox_model_conv_prs_combined, type = "risk", newdata = test_data)
train_data$ypred_conv_prs_combined <- predict(cox_model_conv_prs_combined, type = "risk", newdata = train_data)
c_conv_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs_combined, data = test_data)

test_data$sex <- ifelse(test_data$a2 == 1, "male", "female")
male_data_test <- test_data[test_data$a2 == 1,]
female_data_test <- test_data[test_data$a2 == 0,]

c_conv <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv, data = male_data_test)
c_protein <- survConcordance(Surv(followup_dr, incident_dr) ~ score, data = male_data_test)
c_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_combined, data = male_data_test)
c_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs, data = male_data_test)
c_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs_combined, data = male_data_test)
c_conv_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs, data = male_data_test)
c_conv_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs_combined, data = male_data_test)

c_conv <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv, data = female_data_test)
c_protein <- survConcordance(Surv(followup_dr, incident_dr) ~ score, data = female_data_test)
c_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_combined, data = female_data_test)
c_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs, data = female_data_test)
c_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs_combined, data = female_data_test)
c_conv_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs, data = female_data_test)
c_conv_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs_combined, data = female_data_test)

test_data$PRS <- ifelse(!is.na(test_data$prs) & (test_data$prs >= median(data_pro$prs, na.rm = T)), "High",
                        ifelse(!is.na(test_data$prs) & (test_data$prs < median(data_pro$prs, na.rm = T)), "Low", NA))
PRShigh_data_test <- test_data[!is.na(test_data$prs) & (test_data$prs >= median(data_pro$prs, na.rm = T)),]
PRSlow_data_test <- test_data[!is.na(test_data$prs) & (test_data$prs < median(data_pro$prs, na.rm = T)),]

c_conv <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv, data = PRShigh_data_test)
c_protein <- survConcordance(Surv(followup_dr, incident_dr) ~ score, data = PRShigh_data_test)
c_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_combined, data = PRShigh_data_test)
c_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs, data = PRShigh_data_test)
c_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs_combined, data = PRShigh_data_test)
c_conv_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs, data = PRShigh_data_test)
c_conv_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs_combined, data = PRShigh_data_test)

c_conv <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv, data = PRSlow_data_test)
c_protein <- survConcordance(Surv(followup_dr, incident_dr) ~ score, data = PRSlow_data_test)
c_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_combined, data = PRSlow_data_test)
c_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs, data = PRSlow_data_test)
c_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_prs_combined, data = PRSlow_data_test)
c_conv_prs <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs, data = PRSlow_data_test)
c_conv_prs_combined <- survConcordance(Surv(followup_dr, incident_dr) ~ ypred_conv_prs_combined, data = PRSlow_data_test)

model_list <- c("Clinical", "Polygenic", "Clinical + Polygenic", "Clinical + Proteomic", "Polygenic + Proteomic", "Clinical + Polygenic + Proteomic")
dca_overall <- dca(cox_model_conv, cox_model_prs, cox_model_conv_prs, cox_model_combined, cox_model_prs_combined, c_conv_prs_combined,
                   model.names = model_list,
                   new.data = test_data)
dca_male <- dca(cox_model_conv, cox_model_prs, cox_model_conv_prs, cox_model_combined, cox_model_prs_combined, c_conv_prs_combined,
                   model.names = model_list,
                   new.data = male_data_test)
dca_female <- dca(cox_model_conv, cox_model_prs, cox_model_conv_prs, cox_model_combined, cox_model_prs_combined, c_conv_prs_combined,
                   model.names = model_list,
                   new.data = female_data_test)
dca_PRShigh <- dca(cox_model_conv, cox_model_prs, cox_model_conv_prs, cox_model_combined, cox_model_prs_combined, c_conv_prs_combined,
                   model.names = model_list,
                   new.data = PRShigh_data_test)
dca_PRSlow <- dca(cox_model_conv, cox_model_prs, cox_model_conv_prs, cox_model_combined, cox_model_prs_combined, c_conv_prs_combined,
                   model.names = model_list,
                   new.data = PRSlow_data_test)


library(survival)
library(survIDINRI)
library(PredictABEL)
covariates_all_formula
logistic.model.list <- list(Basic_model = glm(incident_dr ~ as.factor(a2) + a16 + dm_duration + a18 + a10 + TC_divide_HDL + baselineage + prev_hbp + as.factor(a6), data = train_data, family = binomial),
                            New_model = glm(incident_dr ~ as.factor(a2) + a16 + dm_duration + a18 + a10 + TC_divide_HDL + baselineage + prev_hbp + as.factor(a6), data = train_data, family = binomial))
reclassification(data = test_data, cOutcome = which(names(test_data) == "incident_dr"),
                 predrisk1 = fitted(logistic.model.list[["Basic_model"]]),
                 predrisk2 = fitted(logistic.model.list[["New_model"]]),
                 cutoff = c(0, 0.3, 1))


test_data$dr_quantile_risk <- as.vector(ifelse(test_data$score >= quantile(test_data$score, 0.8), "HIGH",
                                               ifelse(test_data$score < quantile(test_data$score, 0.2),"LOW",
                                                      ifelse(test_data$score >= quantile(test_data$score, 0.4) & test_data$score < quantile(test_data$score, 0.6),"MIDDLE",
                                                             NA))))
library(survival)
library(survminer)
library(ggthemes)
diff <- survdiff(Surv(followup_dr, incident_dr) ~ dr_quantile_risk, data = test_data)
pValue <- 1 - pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 3)
pValue <- format(pValue, scientific = TRUE)
test_data$dr_quantile_risk <- factor(test_data$dr_quantile_risk, levels=c("HIGH","MIDDLE","LOW"))
fit <- survfit(Surv(followup_dr, incident_dr) ~ dr_quantile_risk, data = test_data)




#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# DR——OCTA----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————

library(haven)
library(readr)
library(tidyr)
library(dplyr)
library(openxlsx)

octa_varlist <- setdiff(names(octa_varlist), c("study", "no", "year", "eye", "id")) # 获取血流的变量列表
bloodflow_outcomes <- paste(octa_varlist, "_change", sep = "") # 血流结局清单
proteins <- c("plxnb2","gdf15","ren")

library(openxlsx)
wb <- createWorkbook()

lm_results_list <- list()

for (protein in proteins) {
  lm_result <- data.frame()
  
  # 循环每个代谢物并存储结果
  for (bloodflow_outcome in bloodflow_outcomes) {
    if (sum(!is.na(data_all_6mm[[bloodflow_outcome]])) > 50) {
      formula <- paste0(bloodflow_outcome, "~", protein, " + age + as.factor(sex) + hba1c + dm_duration + bmi + sbp")
      lm_model <- lm(formula, data = data_all_6mm)
      coef <- summary(lm_model)$coefficients[, 1][2]
      std <- summary(lm_model)$coefficients[, 2][2]
      tvalue <- summary(lm_model)$coefficients[, 3][2]
      low_limit <- coef - 1.96*std
      up_limit <- coef + 1.96*std
      pvalue <- summary(lm_model)$coefficients[, 4][2]
      sig <- ifelse(pvalue <= 0.05, "Significant", "-")
      
      formula <- paste0(bloodflow_outcome, "_top ~", protein, " + age + as.factor(sex) + hba1c + dm_duration + bmi + sbp")
      logit_model <- glm(formula, data = data_all_6mm, family = binomial(link = "logit"))
      coef_c <- summary(logit_model)$coefficients[, 1][2]
      std_c <- summary(logit_model)$coefficients[, 2][2]
      zvalue_c <- summary(logit_model)$coefficients[, 3][2]
      HR_c <- exp(coef_c)
      HR_low_c <- exp(coef_c - 1.96*std_c)
      HR_up_c <- exp(coef_c + 1.96*std_c)
      pvalue_c <- summary(logit_model)$coefficients[, 4][2]
      sig_c <- ifelse(pvalue_c <= 0.05, "Significant", "-")
      
      lm_result <- rbind(lm_result, c(bloodflow_outcome, coef, low_limit, up_limit, std, pvalue, sig,
                                      bloodflow_outcome, coef_c, HR_c, HR_low_c, HR_up_c, std_c, pvalue_c, sig_c))
    } else {
      # 如果所有行都为空值，则直接输出空值
      lm_result <- rbind(lm_result, c(bloodflow_outcome, NA, NA, NA, NA, NA, NA,
                                      bloodflow_outcome, NA, NA, NA, NA, NA, NA, NA))
    }
  }
  
  colnames(lm_result) <- c("bloodflow", "coef", "low_limit", "up_limit", "std", "pvalue", "sig",
                           "bloodflow_c", "coef_c", "HR_c", "HR_low_c", "HR_up_c", "std_c", "pvalue_c", "sig_c")
  
  columns_to_convert <- c("coef", "low_limit", "up_limit", "std", "pvalue", "coef_c", "HR_c", "HR_low_c", "HR_up_c", "std_c", "pvalue_c") #转换数值变量
  lm_result[columns_to_convert] <- lapply(lm_result[columns_to_convert], as.numeric)
  
  lm_results_list[[protein]] <- lm_result
  addWorksheet(wb, protein)
  writeData(wb, protein, lm_results_list[[protein]])
}



#——————————————————————————————————————————————————————————————————————————————————————————————————————————####
# 孟德尔随机化----------------------------------------------------------------------------------------------------
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————

library(data.table)    
library(dplyr)         
library(TwoSampleMR)   

exposure_iv <- fread("PREPARED_CIS_PQTL_FILE.txt") %>% as.data.frame()

clump_in <- exposure_iv %>%
  transmute(
    rsid          = SNP,
    chr           = as.integer(CHR),
    pos           = POS,
    pval          = P,
    beta          = BETA,
    se            = SE,
    effect_allele = EFFECT_ALLELE,
    other_allele  = OTHER_ALLELE,
    eaf           = EAF,
    samplesize    = SAMPLE_SIZE
  )

dat_clumped <- ld_clump(
  d         = clump_in,
  bfile     = "PATH_TO_REFERENCE_PANEL",  
  plink_bin = get_plink_exe(),
  clump_r2  = 0.001
)

exposure_dat <- dat_clumped %>%
  transmute(
    SNP                    = rsid,
    chr                    = chr,
    pos                    = pos,
    effect_allele.exposure = effect_allele,
    other_allele.exposure  = other_allele,
    eaf.exposure           = eaf,
    beta.exposure          = beta,
    se.exposure            = se,
    pval.exposure          = pval,
    samplesize.exposure    = samplesize,
    id.exposure            = "EXPOSURE_ID",
    exposure               = "EXPOSURE"
  )

outcome_raw <- fread("OUTCOME_GWAS_FILE.txt") %>% as.data.frame()
outcome_dat <- outcome_raw %>%
  transmute(
    SNP                   = SNP,
    effect_allele.outcome = EFFECT_ALLELE,
    other_allele.outcome  = OTHER_ALLELE,
    eaf.outcome           = EAF,
    beta.outcome          = BETA,
    se.outcome            = SE,
    pval.outcome          = P,
    samplesize.outcome    = SAMPLE_SIZE,
    id.outcome            = "OUTCOME_ID",
    outcome               = "OUTCOME"
  )

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)

mr_res <- mr(
  harmonised_dat,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
)

mr_presso <- TwoSampleMR::run_mr_presso(harmonised_dat)

mr_het <- mr_heterogeneity(harmonised_dat)

mr_pleio <- mr_pleiotropy_test(harmonised_dat)

exposure_rev_raw <- fread("OUTCOME_GWAS_FILE.txt") %>% as.data.frame()
clump_in_rev <- exposure_rev_raw %>%
  transmute(
    rsid          = SNP,
    chr           = as.integer(CHR),
    pos           = POS,
    pval          = P,
    beta          = BETA,
    se            = SE,
    effect_allele = EFFECT_ALLELE,
    other_allele  = OTHER_ALLELE,
    eaf           = EAF,
    samplesize    = SAMPLE_SIZE
  )

dat_clumped_rev <- ld_clump(
  d         = clump_in_rev,
  bfile     = "PATH_TO_REFERENCE_PANEL",  
  plink_bin = get_plink_exe(),
  clump_r2  = 0.001
)

exposure_dat_rev <- dat_clumped_rev %>%
  transmute(
    SNP                    = rsid,
    chr                    = chr,
    pos                    = pos,
    effect_allele.exposure = effect_allele,
    other_allele.exposure  = other_allele,
    eaf.exposure           = eaf,
    beta.exposure          = beta,
    se.exposure            = se,
    pval.exposure          = pval,
    samplesize.exposure    = samplesize,
    id.exposure            = "OUTCOME_ID",
    exposure               = "OUTCOME"
  )

outcome_rev_raw <- fread("CIS_PQTL_FILE.txt") %>% as.data.frame()
outcome_dat_rev <- outcome_rev_raw %>%
  transmute(
    SNP                   = SNP,
    effect_allele.outcome = EFFECT_ALLELE,
    other_allele.outcome  = OTHER_ALLELE,
    eaf.outcome           = EAF,
    beta.outcome          = BETA,
    se.outcome            = SE,
    pval.outcome          = P,
    samplesize.outcome    = SAMPLE_SIZE,
    id.outcome            = "EXPOSURE_ID",
    outcome               = "EXPOSURE"
  )

harmonised_dat_rev <- harmonise_data(exposure_dat_rev, outcome_dat_rev)

mr_res_rev <- mr(
  harmonised_dat_rev,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
)

mr_presso_rev <- TwoSampleMR::run_mr_presso(harmonised_dat_rev)

mr_het_rev <- mr_heterogeneity(harmonised_dat_rev)

mr_pleio_rev <- mr_pleiotropy_test(harmonised_dat_rev)
















