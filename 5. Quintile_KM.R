library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(readr)

Big_dataframe <- read_csv("~/Full_dataframe.csv")

KRT <- Big_dataframe[,c("eid", "KRT18")]
KRT <- KRT %>%
  mutate(group = ntile(KRT18, 5))  

CDHR <- Big_dataframe[,c("eid", "CDHR2")]
CDHR <- CDHR %>%
  mutate(group = ntile(CDHR2, 5))  

GGT <- Big_dataframe[,c("eid", "GGT1")]
GGT <- GGT %>%
  mutate(group = ntile(GGT1, 5))  

FUOM <- Big_dataframe[,c("eid", "FUOM")]
FUOM <- FUOM %>%
  mutate(group = ntile(FUOM, 5))  

ACY1 <- Big_dataframe[,c("eid", "ACY1")]
ACY1 <- ACY1 %>%
  mutate(group = ntile(ACY1, 5))  

FigPath = "~/KM/plot5/"
DataPath = "~/KM/plot5/"

TargetPro <- "KRT18"
ProData_name <- paste(DataPath, TargetPro, ".csv", sep = "")

ProData <- read_csv(ProData_name)
MASLDtime <- Big_dataframe[, c(1, 47, 51)]

tmpgroup <- inner_join(ProData, MASLDtime, by = "eid")
tmpgroup$group <- factor(tmpgroup$group, levels = c(1,2,3,4,5), labels = c("Q1", "Q2","Q3","Q4","Q5"))

sfit <- survfit(Surv(time, MASLD) ~ group, data = tmpgroup)
cox_model <- coxph(Surv(time, MASLD) ~ group, data = tmpgroup)

summary_df <- data.frame(summary(cox_model)$coefficients)
confint_df <- data.frame(exp(confint(cox_model, level = 0.95)))
summary_df$exp_coef <- exp(summary_df[,1])
summary_df$exp_lower<- confint_df[,1]
summary_df$exp_upper<- confint_df[,2]

summary_df$allevent<-cox_model$nevent
summary_df$total<-cox_model$n
colnames(summary_df)
HR<- round(summary_df$exp_coef[4],2)
CI_lower<-round(summary_df$exp_lower[4],2)
CI_upper<-round(summary_df$exp_upper[4],2)
p_value1<-summary_df$Pr...z..[4]
p_value <- format(p_value1, scientific = TRUE, digits = 3)

custom_theme <- function() {
  theme(
    text = element_text(family = "Candara"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    risktable.line = element_line(color = "white"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, family = "Candara"),
    axis.title = element_text(size = 12, family = "Candara"),
    axis.text = element_text(size = 10, family = "Candara"),
    legend.text = element_text(size = 8, family = "Candara"),
    legend.title = element_text(size = 12, family = "Candara")
  )
}

p <- ggsurvplot(sfit,
                pval = FALSE,
                conf.int = TRUE,
                conf.int.alpha = 0.3,
                fun = "cumhaz",
                xlab = "Follow-up time (years)",
                ylab = "Cumulative incidence",
                title = TargetPro,
                title.y = 1.5,
                palette = c("#729EC8", "#0ED1D4", "#A5F390","#F2CC45" , "#F27945"),
                size = 0.6,
                legend.title = ggplot2::element_blank(),
                ylim = c(0, 0.06),
                xlim = c(0, 15),
                break.y.by = 0.02,
                break.x.by = 2.5,
                risk.table = TRUE,
                risk.table.y.text = TRUE,
                risk.table.fontsize = 3.5,
                risk.table.height = 0.35,
                ggtheme = custom_theme()) 

# Add annotation for HR and P value
p$plot <- p$plot +
  annotate("text", x = 0.4, y = 0.060, family = "Candara",
           label = paste("HR(Q5 vs Q1) =", HR, "[", CI_lower, "-", CI_upper, "]"),
           hjust = 0, vjust = 1, size = 4) +
  annotate("text", x = 0.4, y = 0.054, family = "Candara",
           label = paste("P value(Q5 vs Q1) =", p_value),
           hjust = 0, vjust = 1, size = 4)


# Adjust legend position and add custom legend
p$plot <- p$plot +
  theme(legend.position = c(0.04, 0.85),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 11))

# Print the plot
print(p)
