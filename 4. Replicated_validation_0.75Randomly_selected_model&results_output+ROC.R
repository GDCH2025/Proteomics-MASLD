#paralleled model&results output:75:25 randomly selected
library(dplyr)
library(survival)
library(pROC)
library(readr)
library(parallel)

# --------- roots ----------
data_file <- "~/All_data.csv"
region_file <- "~/54.csv"
outcome_list <- c("MASLD","MASLD5","MASLD10","MASLD16")

ukb_root <- "~/ukb_randomly0.75"
model_root_dir   <- file.path(ukb_root, "models")   
results_root_dir <- file.path(ukb_root, "results")  
plot_root_dir    <- file.path(ukb_root, "plots")    

R_boot <- 1000
min_n_for_boot <- 50
min_events_per_boot <- 5
n_cores <- detectCores() - 1
set.seed(123)
# ------------------------

dir.create(model_root_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root_dir, recursive = TRUE, showWarnings = FALSE)

All <- read_csv(data_file)
X54 <- read_csv(region_file)[, c(1,2)]
names(X54)[2] <- "Region"
All <- left_join(All, X54, by = "eid")

set.seed(123)  
train_idx <- sample(seq_len(nrow(All)), size = 0.75 * nrow(All))
train_data <- All[train_idx, ]
test_data  <- All[-train_idx, ]

clinical_factors <- "age + sex + Ethnicity + education + TDindex + BMI + diabetes + Hypertension + activity + smoking + alcohol + income + healthy_diet"
proteins_5 <- c("FUOM","ACY1","GGT1","CDHR2","KRT18")
proteins_30 <- unique(c(proteins_5, "BST2","SCLY","PCBD1","HAO1","AGXT","PALM2","IGSF3","ACE2","FTCD","MME","ACAA1","GSTA1","GRPEL1","KYNU","THBS2","CEACAM1","PDZK1","CBS","UPB1","ADH4","PLA2G15","PTS","ADAMTSL2","ITGBL1","GOT1"))

formulas_for_index <- function(outcome_var) {
  formulas <- list()
  labels <- c()
  for (p in proteins_5) {
    formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", p))
    labels <- c(labels, paste0("Single protein: ", p))
  }
  for (p in proteins_5) {
    formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", p, " + ", clinical_factors))
    labels <- c(labels, paste0("Single protein + Clinical: ", p))
  }
  formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", paste(proteins_5, collapse = " + ")))
  labels <- c(labels, "5 proteins")
  formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", paste(proteins_5, collapse = " + "), " + ", clinical_factors))
  labels <- c(labels, "5 proteins + Clinical")
  formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", paste(proteins_30, collapse = " + ")))
  labels <- c(labels, "30 proteins")
  formulas[[length(formulas)+1]] <- as.formula(paste0("Surv(time, ", outcome_var, " == 1) ~ ", paste(proteins_30, collapse = " + "), " + ", clinical_factors))
  labels <- c(labels, "30 proteins + Clinical")
  return(list(formulas = formulas, labels = labels))
}

bootstrap_auc_on_data <- function(model, data_df, outcome_var, R = 1000, min_n = 50, min_events = 5, cores = 1) {
  if (nrow(data_df) < min_n || sum(data_df[[outcome_var]]) < min_events) 
    return(list(mean = NA, sd = NA, ci = c(NA,NA), nsamp = nrow(data_df)))
  
  aucs <- mclapply(1:R, function(r) {
    idx <- sample(seq_len(nrow(data_df)), replace = TRUE)
    d <- data_df[idx, , drop=FALSE]
    if (length(unique(d[[outcome_var]])) < 2 || sum(d[[outcome_var]]) < min_events) return(NA)
    pred <- tryCatch(predict(model, newdata = d, type = "lp"), error = function(e) rep(NA, nrow(d)))
    if (all(is.na(pred)) || length(unique(pred)) < 2) return(NA)
    roc_obj <- tryCatch(pROC::roc(d[[outcome_var]], pred, quiet = TRUE), error = function(e) NULL)
    if (is.null(roc_obj)) return(NA)
    return(as.numeric(pROC::auc(roc_obj)))
  }, mc.cores = cores)
  
  aucs2 <- unlist(aucs)
  aucs2 <- aucs2[!is.na(aucs2)]
  if (length(aucs2) == 0) return(list(mean=NA, sd=NA, ci=c(NA,NA), nsamp=nrow(data_df)))
  return(list(
    mean = mean(aucs2),
    sd = sd(aucs2),
    ci = quantile(aucs2, c(0.025, 0.975)),
    nsamp = nrow(data_df)
  ))
}

# --------- Main ----------
for (outcome in outcome_list) {
  cat("\n=== Processing outcome:", outcome, "===\n")
  ff <- formulas_for_index(outcome)
  formulas <- ff$formulas
  labels <- ff$labels
  
  out_mod_dir <- file.path(model_root_dir, outcome)
  dir.create(out_mod_dir, recursive = TRUE, showWarnings = FALSE)
  out_result_file <- file.path(results_root_dir, paste0("UKB_test_bootstrap_results_", outcome, ".csv"))
  
  results_list <- list()
  for (i in seq_along(formulas)) {
    cat(sprintf(" Fitting model %02d / %02d: %s\n", i, length(formulas), labels[i]))
    f <- formulas[[i]]
    
    formula_text <- paste(deparse(f), collapse="")
    label_text <- gsub("\n"," ", labels[i])
    
    model <- tryCatch(coxph(f, data = train_data), error = function(e) NULL)
    if (is.null(model)) {
      results_list[[i]] <- data.frame(Model=paste0("model_", sprintf("%02d", i)), Formula=formula_text, Label=label_text,
                                      Train_AUC_mean=NA, Train_AUC_sd=NA, Train_CI_low=NA, Train_CI_high=NA, Train_N=nrow(train_data),
                                      Test_AUC_mean=NA, Test_AUC_sd=NA, Test_CI_low=NA, Test_CI_high=NA, Test_N=nrow(test_data),
                                      stringsAsFactors = FALSE)
      next
    }
    
    vars_used <- setdiff(all.vars(f), c("Surv","time", outcome))
    model$call <- NULL; model$y <- NULL; model$linear.predictors <- NULL; model$weights <- NULL
    
    save_file <- file.path(out_mod_dir, sprintf("%s_model_%02d.RData", outcome, i))
    tryCatch(save(model, vars_used, formula_text, file = save_file), error=function(e) message("Save failed: ", e$message))
    
    # Bootstrap AUC
    bs_train <- bootstrap_auc_on_data(model, train_data %>% select(all_of(c("time", outcome, vars_used))), outcome,
                                      R = R_boot, min_n = min_n_for_boot, min_events = min_events_per_boot, cores = n_cores)
    bs_test  <- bootstrap_auc_on_data(model, test_data  %>% select(all_of(c("time", outcome, vars_used))), outcome,
                                      R = R_boot, min_n = min_n_for_boot, min_events = min_events_per_boot, cores = n_cores)
    
    results_list[[i]] <- data.frame(Model=paste0("model_", sprintf("%02d", i)), Formula=formula_text, Label=label_text,
                                    Train_AUC_mean=bs_train$mean, Train_AUC_sd=bs_train$sd, Train_CI_low=bs_train$ci[1], Train_CI_high=bs_train$ci[2], Train_N=bs_train$nsamp,
                                    Test_AUC_mean=bs_test$mean, Test_AUC_sd=bs_test$sd, Test_CI_low=bs_test$ci[1], Test_CI_high=bs_test$ci[2], Test_N=bs_test$nsamp,
                                    stringsAsFactors = FALSE)
  }
  
  # CSV
  final_df <- do.call(rbind, results_list)
  write.csv(final_df, out_result_file, row.names = FALSE)
  cat("Saved results for", outcome, "->", out_result_file, "\n")
  
  # ---------- ROC ----------
  year_text <- gsub("MASLD", "", outcome)
  if (year_text=="") year_text <- "baseline"
  
  colors_train <- c("#F3D030","#F19464","#B46034","#DD7C82","#be71b1","#736baf","#0082bf","#5ECFE1","#00C9AC","#6A8E50","#A1D99B","#9E9AC8","#FD8D3C","#6BAED6")
  colors_test  <- c(
    "#F79C82",  
    "#9E0000",  
    "#FEE08B",  
    "#A67F4F",  
    "#4169E1",  
    "#313695",  
    "#74ADD1",  
    "#A6CEE3",  
    "#E41A1C",  
    "#66C2A5",  
    "#F46D43",  
    "#3288BD",  
    "#FDAE61",  
    "#5E4FA2"   
  )
  
  for (dataset in c("train","test")) {
    if (dataset=="train") {
      ddata <- train_data
      title_text <- paste0("train_", year_text, "_years_ROC")
      colors_use <- colors_train
    } else {
      ddata <- test_data
      title_text <- paste0("test_", year_text, "_years_ROC")
      colors_use <- colors_test
    }
    
    pdf(file.path(plot_root_dir, paste0(title_text,".pdf")), width=7.5, height=7.5)  # 900 px ≈ 7.5 inch * 120 dpi
    tryCatch({
      par(family="sans", font.main=2, font.lab=2, font.axis=2, cex.main=1.5, cex.lab=1.2, cex.axis=1.2)
      plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", main=title_text)
      abline(0,1,col="red",lty=2,lwd=1)
      
      auc_values <- c()
      color_index <- 1
      
      for (i in seq_along(formulas)) {
        model_file <- file.path(out_mod_dir, sprintf("%s_model_%02d.RData", outcome, i))
        if (!file.exists(model_file)) next
        
        e <- new.env()
        load(model_file, envir = e)
        if (!exists("model", envir = e)) next
        model_obj <- get("model", envir = e)
        
        pred <- tryCatch(predict(model_obj, newdata=ddata, type="lp"), error=function(e) NULL)
        if (is.null(pred)) next
        roc_obj <- roc(ddata[[outcome]], pred, quiet=TRUE)
        auc_val <- round(auc(roc_obj),3)
        
        auc_values[labels[i]] <- auc_val
        
        sp <- smooth.spline(roc_obj$specificities, roc_obj$sensitivities, spar=0.5)
        lines(1-sp$x, sp$y, col=colors_use[color_index], lwd=1.5)
        color_index <- color_index + 1
      }
      
      legend_labels <- sapply(seq_along(auc_values), function(i) {
        x <- names(auc_values)[i]
        auc_val <- auc_values[i]
        
        x_new <- sub("^Single protein: ", "", x)
        x_new <- sub("^Single protein \\+ Clinical: (.*)", "\\1 + Clinical", x_new)
        
        if (i == 1) {
          paste0(x_new, " (AUC=", auc_val, ")")
        } else {
          paste0(x_new, " (", auc_val, ")")
        }
      })
      
      legend("bottomright", legend = legend_labels,
             col = colors_use[1:length(auc_values)], lty=1, lwd=1.5, cex=1)
      
    }, error=function(e){message("Error in plotting: ", e$message)},
    finally={dev.off()})
  }
}