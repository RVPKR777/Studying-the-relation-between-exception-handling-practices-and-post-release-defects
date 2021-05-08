suppressMessages(require(rms, quietly = TRUE, warn.conflicts = FALSE))
require(splines, quietly = TRUE)
require(plotly, quietly = TRUE, warn.conflicts = FALSE)
require(Hmisc, quietly = TRUE)
require(e1071, quietly = TRUE)
require(caret, quietly = TRUE)
require(BiodiversityR, quietly = TRUE)
require(logistf, quietly = TRUE)

myCorrelationsByProject <- function(x, projects) {
  temp_list = vector("list", 0)
  
  for (i in 1:length(projects)) {
    temp_project_data <- x[x[,"Project"] == projects[i],]
    print(paste("Project:", projects[i]))
    
    temp_list <- c(temp_list, list(myCorrelations(temp_project_data)))
  }
  return (temp_list)
}

myCorrelations <- function(x) {
  
  # x = as.data.frame(data_list_m0[i])
  
  # Remove identification columns and response variable.
  temp_independent = x[,!(names(x) %in% dropToPredict )]
  ncol(temp_independent)
  
  # Remove predictors that are constant values.
  temp_independent_var = temp_independent[,apply(temp_independent, 2, var, na.rm=TRUE) >= 0.01]
  ncol(temp_independent_var)
  
  # Plot varclus for visuals
  temp_project_varclus = varclus(~ ., data=temp_independent_var, trans="abs")
  plot(temp_project_varclus, main = projects[i])
  abline(h = 1 - CORR_THRESHOLD, col="grey", lty=2)
  
  # Remove highly correlated
  temp_correlations <- cor(temp_independent_var, method="spearman")
  temp_high_corr_names <- findCorrelation(temp_correlations, cutoff= CORR_THRESHOLD, names = TRUE)
  temp_low_cor = temp_independent_var[!(names(temp_independent_var) %in% temp_high_corr_names)]
  ncol(temp_low_cor)
  
  # Put back identification columns and independant variables.
  output = x[(names(x) %in% c(names(temp_low_cor), dropToPredict) )]
  ncol(output)
  
  # Print all that is being kept
  print(paste("NumberOfMetricsInitial:", ncol(temp_independent)))
  print(paste("NumberOfMetricsKept:", ncol(temp_low_cor)))
  print(paste(names(temp_low_cor),collapse=" + "))
  
  return (output)
}

myCorrelationsReRun <- function(x) {
  
  #x = data
  
  # Remove identification columns and response variable.
  temp_independent = x[,!(names(x) %in% dropToPredict )]
  ncol(temp_independent)
  nrow(temp_independent)
  
  budget = floor(nrow(temp_independent) / 15)
  nvar_initial = ncol(temp_independent)
  nvar = nvar_initial
  
  if(nvar > budget) {
    # Plot varclus for visuals
    temp_project_varclus = varclus(~ ., data=temp_independent, trans="abs")
    plot(temp_project_varclus, main = projects[i])
    abline(h = 1 - CORR_THRESHOLD, col="grey", lty=2)
    
    # Remove until match budget
    temp_low_cor_info <- myfindCorrelation(temp_independent)
    temp_low_cor = temp_independent[(names(temp_independent) %in% temp_low_cor_info$vars)]
    nvar = ncol(temp_low_cor)
    
    final_cutoff = temp_low_cor_info$cutoff
    abline(h = 1 - final_cutoff, col="red", lty=2)
    
    # Put back identification columns and independant variables.
    output = x[(names(x) %in% c(names(temp_low_cor), dropToPredict) )]
    ncol(output)
    
    print(paste(names(temp_low_cor),collapse=" + "))
    
  } else {
    output = x
    final_cutoff = CORR_THRESHOLD
  }
  
  list(data=output,nvar_initial=nvar_initial,budget=budget, over_budget=(nvar_initial > budget), cutoff=final_cutoff,
       nvar_final=nvar)
}

myfindCorrelation <- function(x, opvExp = 15, names = TRUE){
  
  # x = temp_independent
  
  budget = floor(nrow(x) / opvExp)
  cors <- cor(x, method="spearman")
  
  dyn_cutoff = CORR_THRESHOLD # initial
  vars = names(x)
  nvar = ncol(x)
  while( nvar > budget){
    dyn_cutoff = dyn_cutoff - 0.001
    
    temp_high_corr_names <- findCorrelation(cors, cutoff=dyn_cutoff , names = TRUE)
    vars = vars[!(vars %in% temp_high_corr_names)]
    nvar = length(vars)
  }
  list(vars=vars,cutoff=dyn_cutoff)
}

myMissingData <- function(x) {  
  
  temp_naclus = naclus(x)
  plot(temp_naclus)
  
  if (sum(temp_naclus$sim < NA_THRESHOLD)) {
    temp_project_data_omitted <- na.omit(x)
    print(paste("MissingData: Small fraction found: Single imputation or discard. Discard is chosen. Discarded rows:",nrow(x)-nrow(temp_project_data_omitted), "/", nrow(x)))
  } else {
    # TODO: implement multiple imputation
    temp_project_data_omitted = x
    print(paste("MissingData: High fraction found: Multiple imputation needed! Discarded rows:",nrow(x)-nrow(temp_project_data_omitted), "/", nrow(x)))
  }
  return(as.data.frame(temp_project_data_omitted))
}

modelSelectionAndNormalityAdjustment <- function(data_list) { 
  
  ### Estimate budget for degrees of freedom (MC-2)
  # As we have the missing data taken care of, we can evaluate our budget for degrees of freedom. The modeling procedure can afford a certain amount of degrees of freedom based on its size. A degree of freedom will be then allocated according to the model formulation and number of variables.
  for (i in 1:length(projects)) {
    budget = floor(nrow(data_list[[i]]$data) / 15)
    print(paste("Project:", projects[i],"D.F. Budget:", budget))
  }
  
  ### Selecting Appropriate Model and Normality adjustment
  # Since the number of observations are defined, we evaluate the dependent variable, the number of post-release bugs to understand its distribution and choose the type of model to be used.
  
  #### Skewness and Kurtosis (MC-2)
  for (i in 1:length(projects)) {
    post_bugs_skewness = skewness(data_list[[i]]$data$Distinct.count.of.Issue.Key.POST)
    print(paste("Project:", projects[i],"skewness,", post_bugs_skewness))
    post_bugs_kurtosis = kurtosis(data_list[[i]]$data$Distinct.count.of.Issue.Key.POST)
    print(paste("Project:", projects[i],"kurtosis,", post_bugs_kurtosis))
  }
  
  #### Cumulative distribution
  # The cumulative distribution also represents the data and how much of the data is explained with a same value.
  # As the Ecdf plots show, more than 70% of the files in all projects have 0 bugs. This suggests a poisson distribution, which can be modeled using Ordinary Least Squares (OLS).
  # 
  for (i in 1:length(projects)) {
    plot(Ecdf(~ Distinct.count.of.Issue.Key.POST, data=data_list[[i]]$data))
  }
  
  # #### Normality adjustment and cumulative distribution (MC3)
  # For the poisson distribution using OLS we can adjust the dependant variable using a logarithmic transformation. We apply the logarithm in the base 10 to the values of post-release bugs plus one. Since most of the data is zero, we need the adjustment to one, to avoid infinite values.
  # We can then plot a new cumulative distribution and include the data points from our data set.
  for (i in 1:length(projects)) {
    project_log = log(data_list[[i]]$data$Distinct.count.of.Issue.Key.POST+1, base=10)
    project_log_max = max(project_log)
    project_log_min = min(project_log)
    project_log_mean = mean(project_log)
    project_log_sd = sd(project_log)
    
    # plot(Ecdf(~ project_log, data=data.frame(project_log)))
    xs = seq(floor(project_log_min), ceiling(project_log_max), length=100)
    cdf = pnorm(xs, project_log_mean, project_log_sd)
    plot(xs, cdf, type='l')
    lines(ecdf(project_log))
  }
}

dataApplyReduction <- function(data_list) { 

  ### Correlation analysis (MC4)
  data_list_low_corr = vector("list", 0)
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_data = myCorrelations(as.data.frame(data_list[i]))
    data_list_low_corr <- c(data_list_low_corr, list(temp_data))
  }
  
  ### Redundancy Analysis (MC5)
  data_list_after_redun = vector("list", 0)
  
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_corr = as.data.frame(data_list_low_corr[i])
    temp_corr_ind = temp_corr[,!(names(temp_corr) %in% dropToPredict )]
    redun_obj = redun(~., data = temp_corr_ind, nk=0)
    
    print(paste("Redudant variables:", paste(redun_obj$Out, collapse=",␣")))
    
    after_redun = temp_corr[,!(names(temp_corr) %in% redun_obj$Out)]
    print(paste(names(after_redun),collapse=" + "))
    
    data_list_after_redun <- c(data_list_after_redun, list(as.data.frame(after_redun)))
  }
  
  ### Budget based correlation analysis (MC6)
  data_list_model = vector("list", 0)
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_corr_rerun_info = myCorrelationsReRun(as.data.frame(data_list_after_redun[i]))
    
    print(paste("NumberOfMetricsInitial:", temp_corr_rerun_info$nvar_initial, "Budget:", temp_corr_rerun_info$budget,"Over Budget:", temp_corr_rerun_info$over_budget, "NumberOfMetricsKept:", temp_corr_rerun_info$nvar_final, "CorrelationCutoff:", temp_corr_rerun_info$cutoff))
    data_list_model <- c(data_list_model, list(temp_corr_rerun_info))
  }
  
  data_list_model
}

dataApplyReductionByModel <- function(data_list,model_type) { 
  
  data_list = all_list_omitted_m2
  ### Correlation analysis (MC4)
  data_list_low_corr = vector("list", 0)
  for (i in 1:length(data_list)) {
    name = data_list[[i]]$name
    project = data_list[[i]]$project
    data = data_list[[i]]$data
    sig = data_list[[i]]$sig
    print(paste("Project:", project, "Name:", name))
    
    temp_data = myCorrelations(data)
    print(names(temp_data))
    data_list_low_corr <- c(data_list_low_corr, list(list(name=name, project=project, data=temp_data, sig=sig)))
  }
  
  ### Redundancy Analysis (MC5)
  data_list_after_redun = vector("list", 0)
  for (i in 1:length(data_list_low_corr)) {
    name = data_list_low_corr[[i]]$name
    project = data_list_low_corr[[i]]$project
    data = data_list_low_corr[[i]]$data
    sig = data_list_low_corr[[i]]$sig
    print(paste("Project:", project, "Name:", name))
    
    temp_corr = data
    temp_corr_ind = temp_corr[,!(names(temp_corr) %in% dropToPredict )]
    redun_obj = redun(~., data = temp_corr_ind, nk=0)
    
    print(paste("Redudant variables:", paste(redun_obj$Out, collapse=",␣")))
    
    after_redun = temp_corr[,!(names(temp_corr) %in% redun_obj$Out)]
    print(paste("Final variables:",paste(names(after_redun),collapse=" + ")))
    
    print(paste("NumberOfMetricsInitial:", ncol(temp_corr)))
    print(paste("NumberOfMetricsKept:", ncol(after_redun)))
    
    data_list_after_redun <- c(data_list_after_redun, list(list(name=name, project=project, data=after_redun, sig=sig)))
  }
  
  ### Budget based correlation analysis (MC6)
  data_list_model = vector("list", 0)
  for (i in 1:length(data_list_after_redun)) {
    name = data_list_after_redun[[i]]$name
    project = data_list_after_redun[[i]]$project
    data = data_list_after_redun[[i]]$data
    sig = data_list_after_redun[[i]]$sig
    
    print(paste("------------------------------------Project:", project, "Name:", name))
    temp_corr_rerun_info = myCorrelationsReRun(data)
    print(paste("NumberOfMetricsInitial:", temp_corr_rerun_info$nvar_initial, "Budget:", temp_corr_rerun_info$budget,"Over Budget:", temp_corr_rerun_info$over_budget, "NumberOfMetricsKept:", temp_corr_rerun_info$nvar_final, "CorrelationCutoff:", temp_corr_rerun_info$cutoff))
    
    if( !is.null(sig))
    for(j in 1:length(sig))
      if( !(sig[[j]] %in% colnames(temp_corr_rerun_info$data)) && sig[[j]] != "TOTAL")
      {
        print(paste("Significant missing:", sig[[j]]))
        print("Check varclus!:")
        
        variables_to_keep = names(temp_corr_rerun_info$data)
        
        variables_to_add=c("")
        variables_to_remove=c("")
        if(project=="" && model_type=="")
        {
          variables_to_add = sig
          variables_to_remove=c("")
        } else
        {
          print(paste("Significant ------ STILL ------ missing:", sig[[j]]))
          
          data_no_sig = data[,!(names(data) %in% sig )]
          temp_corr_rerun_info_no_sig = myCorrelationsReRun(data_no_sig)
          
          temp_corr_rerun_info_no_sig_names = names(temp_corr_rerun_info_no_sig$data)
          
          data_w_sig = data[(names(data) %in% c(temp_corr_rerun_info_no_sig_names, sig) )]
          temp_corr_rerun_info_new = myCorrelationsReRun(data_w_sig)
          variables_to_keep = names(temp_corr_rerun_info_new$data)
          
          if(project=="" && model_type=="")
          {
            variables_to_add = sig
            variables_to_remove=c("")
          }
          else
          {
            print(paste("Significant ------ STILL ------ missing:", sig[[j]]))
          }
        }
        
        data_manual_fixed = data[(names(data) %in% c(variables_to_keep, variables_to_add) )]
        data_manual_fixed = data_manual_fixed[,!(names(data_manual_fixed) %in% variables_to_remove )]
        print(paste(names(data_manual_fixed),collapse=" + "))
        
        temp_corr_rerun_info$data = data_manual_fixed
      }
    
    print(paste("NumberOfMetricsInitial:", temp_corr_rerun_info$nvar_initial, "Budget:", temp_corr_rerun_info$budget,"Over Budget:", temp_corr_rerun_info$over_budget, "NumberOfMetricsKept:", temp_corr_rerun_info$nvar_final, "CorrelationCutoff:", temp_corr_rerun_info$cutoff))
    data_list_model <- c(data_list_model, list(list(name=name, project=project, info=temp_corr_rerun_info)))
  }
  
  data_list_model
}

dataSetupFormulasBinary <- function(data_list) { 
  ### Setup formulas
  form_bin_dep = "Distinct.count.of.Issue.Key.POST>0~"
  
  form_list_bin = vector("list", 0)
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_redun = data_list[[i]]$data
    temp_names = names(temp_redun)
    temp_names_ind = temp_names[!(temp_names %in% dropToPredict)]
    
    temp_form_ind = paste(temp_names_ind,collapse=" + ")
    print(temp_form_ind)
    
    form_bin = paste(form_bin_dep,temp_form_ind,sep="")
    
    form_list_bin <- c(form_list_bin, list(form_bin))
  }
  form_list_bin
}

dataSetupFormulasBinaryByModel <- function(data_list) { 
  ### Setup formulas
  form_bin_dep = "Distinct.count.of.Issue.Key.POST>0~"
  
  form_list_bin = vector("list", 0)
  for (i in 1:length(data_list)) {
    name = data_list[[i]]$name
    project = data_list[[i]]$project
    data = data_list[[i]]$info$data
    print(paste("Project:", project, "Name:", name))
    
    temp_redun = data
    temp_names = names(temp_redun)
    temp_names_ind = temp_names[!(temp_names %in% dropToPredict)]
    
    temp_form_ind = paste(temp_names_ind,collapse=" + ")
    print(temp_form_ind)
    
    form_bin = paste(form_bin_dep,temp_form_ind,sep="")
    
    form_list_bin <- c(form_list_bin, list(list(name=name, project=project,form=form_bin)))
  }
  form_list_bin
}

dataSetupFormulas <- function(data_list) { 
  ### Setup formulas
  form_log_dep = "log(Distinct.count.of.Issue.Key.POST + 1, base = 10)~"
  
  form_list_log = vector("list", 0)
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_redun = data_list[[i]]$data
    temp_names = names(temp_redun)
    temp_names_ind = temp_names[!(temp_names %in% dropToPredict)]
    
    temp_form_ind = paste(temp_names_ind,collapse=" + ")
    print(temp_form_ind)
    
    form_log = paste(form_log_dep,temp_form_ind)
    sp = spearman2(as.formula(form_log), data=temp_redun, p=2)
    plot(sp)
    
    form_list_log <- c(form_list_log, list(form_log))
  }
  form_list_log
}

modelFitLogistic <- function(data_list, form_list, model_name) { 

  # Fit regression model (MC7)
  models = vector("list", 0)
  for (i in 1:length(projects)) {
    print(paste("Project:", projects[i]))
    temp_redun = data_list[[i]]$data
    dd <<- datadist(temp_redun); options(datadist='dd')
    
    form_bin = as.character(form_list[i])
    dropText = c("File.Path",  "Project", "Language", "Table.Name","Name","Kind", "File")
    temp_redun_nn = temp_redun[,!(names(temp_redun) %in% dropText )]
    temp_data_log = log10(temp_redun_nn+1)
    fit_lrm = try(lrm(formula=as.formula(form_bin), data=temp_data_log,x=T,y=T), silent=TRUE)
    print(fit_lrm, latex = TRUE)
    models <- c(models, list(list(name=model_name,project=projects[i],fit=fit_lrm,data=temp_data_log, nvar_initial=data_list[[i]]$nvar_initial, cutoff=data_list[[i]]$cutoff)))
  }
  models
}

modelFitLogisticByModel <- function(data_list, form_list, model_name) { 
  # data_list = all_list_model_m2
  # form_list = form_list_bin_m2
  # model_name = "CAT.BSAP"
  
  # Fit regression model (MC7)
  models = vector("list", 0)
  for (i in 1:length(data_list)) {
    name = data_list[[i]]$name
    project = data_list[[i]]$project
    data = data_list[[i]]$info$data
    print(paste("Project:", project, "Name:", name))
    
    temp_redun = data
    dd <<- datadist(temp_redun); options(datadist='dd')
    
    form_bin = form_list[[i]]$form
    dropText = c("File.Path",  "Project", "Language", "Table.Name","Name","Kind", "File")
    temp_redun_nn = temp_redun[,!(names(temp_redun) %in% dropText )]
    temp_data_log = log10(temp_redun_nn+1)
    fit_lrm = try(lrm(formula=as.formula(form_bin), data=temp_data_log,x=T,y=T), silent=TRUE)
    print(fit_lrm, latex = TRUE)
    models <- c(models, list(list(name=model_name,project=project,fit=fit_lrm,data=temp_data_log, nvar_initial=data_list[[i]]$info$nvar_initial, cutoff=data_list[[i]]$info$cutoff)))
  }
  models
}

