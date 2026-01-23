library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(lubridate)
library(survival)
library(broom)

Proteomics_Olink_Instance0 <- read_csv("/mnt/data2/epidemiological/ukb_original/Proteomics_Olink_Instance0.csv")
colnames(Proteomics_Olink_Instance0)[1] = 'eid'
cumulative_cox <- data.frame()

for (i in 2:ncol(Proteomics_Olink_Instance0)) {
  P1 <- Proteomics_Olink_Instance0[, c(1,i)]
  P1$eid<-as.character(P1$eid)
  P1 <-inner_join(P1,X53,by="eid")
  colnames(P1)[2] = 'protein'
  
  All<-inner_join(MASLD4cov_unique,P1)
  
  time<-All %>% mutate(time = as.numeric(
    difftime(endTime, beginTime, units = "days")) / 365.25)
  
  alldata<-filter(time,time>0)
  
  catvars<-c("sex","Ethnicity","education","income","diabetes","Hypertension","healthy_diet","smoking","alcohol")
  alldata[catvars] <- lapply(alldata[catvars], factor)
  
  #model 1
  tryCatch({
    mycox <- coxph(Surv(time, MASLD== 1) ~ protein+age+sex+Ethnicity+education+TDindex+income+BMI+diabetes+Hypertension+healthy_diet+activity+smoking+alcohol, data = alldata)
        test.ph <- cox.zph(mycox)
    test.phtable <- data.frame(test.ph$table) 
    test.phtable <- test.phtable[1,3]
    summary_df <- data.frame(summary(mycox)$coefficients)
    confint_df <- data.frame(exp(confint(mycox, level = 0.95)))
    summary_df$exp_coef <- exp(summary_df[,1])
    summary_df$exp_lower<- confint_df[,1]
    summary_df$exp_upper<- confint_df[,2]
    summary_df<-summary_df[1,5:8]
    colnames(summary_df)[1] = 'Pvalue'
    summary_df$allevent<-mycox$nevent
    summary_df$total<-mycox$n
    summary_df$ph_Test<-test.phtable
    summary_df$exposure<-names(Proteomics_Olink_Instance0)[i]
    summary_df$outcome <- "MASLD"
    
    cumulative_cox <- rbind(cumulative_cox, summary_df)
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
  })
} 
cumulative_cox$Pvalue_bonferroni<-p.adjust(cumulative_cox$Pvalue, method = "bonferroni")
cumulative_cox$Pvalue_fdr<-p.adjust(cumulative_cox$Pvalue, method = "fdr")
write_csv(full_dataframe,file = "~/MASLD_COX_results.csv")

