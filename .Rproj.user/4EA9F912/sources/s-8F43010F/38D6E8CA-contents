suppressMessages(require(rms, quietly = TRUE, warn.conflicts = FALSE))
require(splines, quietly = TRUE)
require(plotly, quietly = TRUE, warn.conflicts = FALSE)
require(Hmisc, quietly = TRUE)
require(e1071, quietly = TRUE)
require(caret, quietly = TRUE)
require(BiodiversityR, quietly = TRUE)
require(logistf, quietly = TRUE)

# Functions

findModel <- function(model_list, model_name, model_project) {
  
  for (i in 1:length(model_list)) {  
    m_name = model_list[[i]]$name
    m_project = model_list[[i]]$project
    
    if (m_name == model_name && m_project == model_project) {
      index = i
      break
    }
  }
  return (index)
}

findProjectData <- function(data_list, model_project) {
  
  for (i in 1:length(data_list)) {  
    m_project = data_list[[i]]$project
    
    if (m_project == model_project) {
      index = i
      break
    }
  }
  return (index)
}

# Basic stats of the initial full model (MC7)
modelStats <- function(model_list) {

  my_model_list = vector("list", 0)

  for (i in 1:length(model_list)) {  
    fit = model_list[[i]]$fit
    project = model_list[[i]]$project
    name = model_list[[i]]$name
    nvar_initial = model_list[[i]]$nvar_initial
    cutoff = model_list[[i]]$cutoff
    file_set = strsplit(name, ".", fixed = TRUE)[[1]][1]
    metric_set = strsplit(name, ".", fixed = TRUE)[[1]][2]
    model_type = strsplit(name, ".", fixed = TRUE)[[1]][3]
    class = as.character(class(fit)[1])
    
    my_model_list <- c(my_model_list,list(list(name=name,project=project,file_set=file_set,metric_set=metric_set,model_type=model_type,nvar_initial=nvar_initial,cutoff=cutoff)))
    index = findModel(model_list,name,project)
    
    if(class != "try-error") 
    {
      my_model_list[[index]][["error"]] <- FALSE
      my_model_list[[index]][["error_msg"]] <- NA
      my_model_list[[index]][["n"]]  <- as.numeric(switch(class
                                                  ,"ols" = fit$stats[[1]]
                                                  ,"lrm" = fit$stats[[1]]
                                                  ,"glm" = nrow(fit$data)
                                                  ,"logistf" = fit$n
                                                  ,NA))
      my_model_list[[index]][["budget"]] <- floor(my_model_list[[index]][["n"]] / 15)
      my_model_list[[index]][["df"]]  <- as.numeric(switch(class
                                                   ,"ols" = fit$stats[[3]]
                                                   ,"lrm" = fit$stats[[4]]
                                                   ,"glm" = sum(anova(fit)$Df, na.rm = TRUE)
                                                   ,"logistf" = fit$df
                                                   ,NA))
      my_model_list[[index]][["over_budget"]]=(my_model_list[[index]][["df"]]>my_model_list[[index]][["budget"]])
      my_model_list[[index]][["lr"]]  <- as.numeric(switch(class
                                                   ,"ols" = fit$stats[[2]]
                                                   ,"lrm" = fit$stats[[3]]
                                                   ,"glm" = NA
                                                   ,"logistf" = as.numeric(strsplit(strsplit(capture.output(fit), " ")[[length(capture.output(fit))-1]][3], "=")[[1]][2])
                                                   ,NA))
      my_model_list[[index]][["r2"]]  <- as.numeric(switch(class
                                                   ,"ols" = fit$stats[[4]]
                                                   ,"lrm" = fit$stats[[10]]
                                                   ,"glm" = as.numeric((strsplit((capture.output(deviancepercentage(fit,data=fit$data, test="Chisq"))[1])," "))[[1]][8])/100
                                                   ,"logistf" = NA
                                                   ,NA))
      my_model_list[[index]][["aic"]]  <- switch(class
                                         ,"ols" = AIC(fit)
                                         ,"lrm" = AIC(fit)
                                         ,"glm" = AIC(fit) # fit$aic
                                         ,"logistf" = extractAIC(fit)[[2]]
                                         ,NA)
      my_model_list[[index]][["form"]] = Reduce(paste, deparse(switch(class
                                                              ,"ols" = fit$sformula[[3]]
                                                              ,"lrm" = fit$sformula[[3]]
                                                              ,"glm" = fit$formula[[3]]
                                                              ,"logistf" = fit$formula[[3]]
                                                              ,NA)))
      
    } else
    {
      my_model_list[[index]][["error"]] <- TRUE
      my_model_list[[index]][["error_msg"]] <- NA
      my_model_list[[index]][["n"]] <- NA
      my_model_list[[index]][["budget"]] <- NA
      my_model_list[[index]][["df"]] <- NA
      my_model_list[[index]][["lr"]] <- NA
      my_model_list[[index]][["r2"]] <- NA
      my_model_list[[index]][["aic"]] <- NA
      my_model_list[[index]][["form"]] <- NA
    }
  }
  return(my_model_list)
}

# Bootstrap - Model stability assessment (MA1)
modelValidate <- function(model_list, my_model_list) {
  
  for (i in 1:length(model_list)) {  
    fit = model_list[[i]]$fit
    
    project = model_list[[i]]$project
    name = model_list[[i]]$name
    class = as.character(class(fit)[1])
    
    index = findModel(model_list,name,project)
    
    if ((class == "ols" || class == "lrm" ) && !my_model_list[[index]][["over_budget"]] )
    {
      v = try(validate(fit, B=1000), silent=TRUE)
      if (class(v) == "validate") {
        my_model_list[[index]][["r2_c"]] <- switch(class
                                               ,"ols" = v["R-square","index.corrected"]
                                               ,"lrm" = v["R2","index.corrected"]
                                               ,NA)
        my_model_list[[index]][["slope_o"]] <- switch(class
                                                  ,"ols" = v["Slope","index.orig"]
                                                  ,"lrm" = v["Slope","index.orig"]
                                                  ,NA)
        my_model_list[[index]][["slope_c"]] <- switch(class
                                                  ,"ols" = v["Slope","index.corrected"]
                                                  ,"lrm" = v["Slope","index.corrected"]
                                                  ,NA)
      }
    } else
    {
      my_model_list[[index]][["r2_c"]] <- NA
      my_model_list[[index]][["slope_o"]] <- NA
      my_model_list[[index]][["slope_c"]] <- NA
    }
  }
  return(my_model_list)
}

# Extract significant metrics from ANOVA using p-value below 0.05
modelSignificance <- function(model_list, my_model_list) {
  
  for (i in 1:length(model_list)) {  
    fit = model_list[[i]]$fit
    
    project = model_list[[i]]$project
    name = model_list[[i]]$name
    class = as.character(class(fit)[1])
    
    index = findModel(model_list,name,project)
    
    if ((class == "ols" || class == "lrm" || class == "glm") && !my_model_list[[index]][["over_budget"]] )
    {
      fitx = fit
      classx = class
      an = anova(fitx, test="Chisq")
      print(paste("project: ", project,"model: ", name))
      print(an)
      signifcant = switch(classx
                          ,"glm" = row.names(an[an[5] <= 0.05 ,])
                          ,"ols" = row.names(an[an[,5] <= 0.05 ,])
                          ,"lrm" = row.names(an[an[,3] <= 0.05 ,])
                          ,NA)
      
      signifcant = paste(signifcant,collapse = ", " )
      
    } else
    {
      signifcant = NA
    }
    
    my_model_list[[index]][["signifcant"]] <- signifcant
  }
  return(my_model_list)
}

# Model simplification (MA2): simplify model using backward selection, re-run fit, re-assess stability, get significant power, get effect on outcome adjusting the means
modelSimplification <- function(model_list, my_model_list) {
  
  # model_list = models_2_BASE
  # my_model_list = model_things_2_BASE
  # # 
  for (i in 1:length(model_list)) {  
    
    fit = model_list[[i]]$fit
    project = model_list[[i]]$project
    temp_data = model_list[[i]]$data
    name = model_list[[i]]$name
    
    file_set = strsplit(name, ".", fixed = TRUE)[[1]][1]
    metric_set = strsplit(name, ".", fixed = TRUE)[[1]][2]
    model_type = strsplit(name, ".", fixed = TRUE)[[1]][3]
    
    class = as.character(class(fit)[1])
    
    print(paste("project: ", project,"model: ", name, "Refit"))
    
    index = findModel(model_list,name,project)
    
    if ((class == "ols" || class == "lrm" ) && !my_model_list[[index]][["over_budget"]])
    {
      # v = try(validate(fit, B=50, bw=T, rule="p", sls=.05, type="individual" ), silent=TRUE)
      # if (class(v) == "validate") {
      # s_fastbw = fastbw(fit, type="individual")
      s_fastbw = fastbw(fit, rule="p", sls=.05, type="individual")
      
      if (length(s_fastbw$names.kept) > 0)
      {
        names = paste(s_fastbw$names.kept, collapse="+")
        if (class == "lrm" )
        {
          form_dep = "Distinct.count.of.Issue.Key.POST>0~"
          final_form = as.formula(paste(form_dep, names, sep=""))
          dd <<- datadist(temp_data); options(datadist='dd')
          refit = try(lrm(formula=as.formula(final_form), data=temp_data,x=T,y=T), silent=TRUE)
        }
        else if (class == "ols")
        {
          form_dep = "log(Distinct.count.of.Issue.Key.POST + 1, base = 10)~"
          final_form = as.formula(paste(form_dep, names, sep=""))
          dd <<- datadist(temp_data); options(datadist='dd')
          refit = try(ols(as.formula(final_form),data=temp_data,x=T,y=T), silent=TRUE)
        }
        if (class(refit)[1] != "try-error") {
          coef = try(refit$coefficients, silent=TRUE)
        } else {
          coef = NA
        }
      }
    } else if ((class == "glm") && !my_model_list[[index]][["over_budget"]])
    {
      temp_data_log = temp_data
      s_step = step(fit, trace=0)
      final_form = s_step$formula
      dd <<- datadist(temp_data_log); options(datadist='dd')
      refit_glm = glm(formula=as.formula(final_form), data=temp_data_log,family=binomial(link = "logit"),control = list(maxit = 50) )
      refit = refit_glm
      if (class(refit_glm)[1] != "try-error") {
        coef = try(refit_glm$coefficients, silent=TRUE)
      } else {
        coef = NA
      }
    } else {
      final_form = NA
      refit = NA
      coef = NA
    }
    
    print(refit)
    class_refit = as.character(class(refit)[1])
    # 
    # print(class_refit)
    
    if(class_refit != "try-error" && !is.na(refit)) {
      error = FALSE
      error_msg_na = NA
      print(paste("project: ", project,"model: ", name, "Refit - summary"))
      
      df_r = as.numeric(switch(class_refit
                               ,"ols" = refit$stats[[3]]
                               ,"lrm" = refit$stats[[4]]
                               ,"glm" = sum(anova(refit)$Df, na.rm = TRUE)
                               ,"logistf" = NA #refit$df
                               ,NA))
      lr_r = as.numeric(switch(class_refit
                               ,"ols" = refit$stats[[2]]
                               ,"lrm" = refit$stats[[3]]
                               ,"glm" = NA
                               ,"logistf" = NA #as.numeric(strsplit(strsplit(capture.output(refit), " ")[[length(capture.output(refit))-1]][3], "=")[[1]][2])
                               ,NA))
      r2_r = as.numeric(switch(class_refit
                               ,"ols" = refit$stats[[4]]
                               ,"lrm" = refit$stats[[10]]
                               ,"glm" = as.numeric((strsplit((capture.output(deviancepercentage(refit,data=refit$data, test="Chisq"))[1])," "))[[1]][8])/100
                               ,"logistf" = NA
                               ,NA))
      aic_r = switch(class_refit
                     ,"ols" = AIC(refit)
                     ,"lrm" = AIC(refit)
                     ,"glm" = AIC(refit) # refit$aic
                     ,"logistf" = NA #extractAIC(refit)[[2]]
                     ,NA)
      
      if (class_refit == "ols" || class_refit == "lrm" ) {
        print(paste("project: ", project,"model: ", name, "Refit - validate"))
        v = try(validate(refit, B=1000), silent=TRUE)
        if (class(v) == "validate") {
          r2_c_r = switch(class_refit
                          ,"ols" = v["R-square","index.corrected"]
                          ,"lrm" = v["R2","index.corrected"]
                          ,NA)
          slope_o_r = switch(class_refit
                             ,"ols" = v["Slope","index.orig"]
                             ,"lrm" = v["Slope","index.orig"]
                             ,NA)
          slope_c_r = switch(class_refit
                             ,"ols" = v["Slope","index.corrected"]
                             ,"lrm" = v["Slope","index.corrected"]
                             ,NA)
        }
      } else {
        r2_c_r = NA
        slope_o_r = NA
        slope_c_r = NA
      }
      
      if (class_refit == "ols" || class_refit == "lrm" || class_refit == "glm") {
        print(paste("project: ", project,"model: ", name, "Refit - anova"))
        fitx = refit
        classx = class_refit
        an = anova(fitx, test="Chisq")
        print(an)
        signifcant_r = switch(classx
                              ,"glm" = row.names(an[an[5] <= 0.05 ,])
                              ,"ols" = row.names(an[an[,5] <= 0.05 ,])
                              ,"lrm" = row.names(an[an[,3] <= 0.05 ,])
                              ,NA)
        
        newdata =  data.frame(project=project, name=name)
        coefs = names(fitx$coefficients)
        for( n in 1:length(coefs)){
          temp_coef = coefs[n]
          if(temp_coef != "Intercept"){
            print(temp_coef)
            newdata[,temp_coef] <- mean(10^temp_data[[temp_coef]]-1)
          }
        }
        newdata = newdata[!(names(newdata) %in% c("project", "name"))]
        print(newdata)
        fixed_at_mean_predicted = predict(fitx, log10(newdata+1), type="fitted.ind")    # Prob(Y=j)
        print(paste("Fixed at Mean:", fixed_at_mean_predicted))
        
        for( n in 1:length(coefs)){
          temp_coef = coefs[n]
          if(temp_coef != "Intercept"){
            original = newdata[,temp_coef]
            newdata[,temp_coef] <- 1.1*original
            temp_predicted = predict(fitx, log10(newdata+1), type="fitted.ind")    # Prob(Y=j)
            print(paste(temp_coef, " Coef at Mean + 10%:", temp_predicted))
            newdata[,temp_coef] = original
          }
        }
        signifcant_r = paste(signifcant_r,collapse = ", " )
      } else
      {
        signifcant_r = NA
      }
    } else
    {
      error = TRUE
      error_msg = capture.output(fit)[1]
      df_r = NA
      lr_r = NA
      r2_r = NA
      aic_r = NA
      
      r2_c_r = NA
      slope_o_r = NA
      slope_c_r = NA
      
      signifcant_r = NA
    }
    my_model_list[[index]][["final_form"]] <- as.character(list(final_form))
    my_model_list[[index]][["coef"]] <- as.character(list(coef))
    my_model_list[[index]][["over_budget_r"]] <- (df_r>my_model_list[[index]][["budget"]])
    my_model_list[[index]][["df_r"]] <- df_r
    my_model_list[[index]][["lr_r"]] <- lr_r
    my_model_list[[index]][["r2_r"]] <- r2_r
    my_model_list[[index]][["aic_r"]] <- aic_r
    my_model_list[[index]][["r2_c_r"]] <- r2_c_r
    my_model_list[[index]][["slope_o_r"]] <- slope_o_r
    my_model_list[[index]][["slope_c_r"]] <- slope_c_r
    my_model_list[[index]][["signifcant_r"]] <- signifcant_r
  }
  return(my_model_list)
}


source_rmd <- function(file, local = FALSE, ...){
  options(knitr.duplicate.label = 'allow')
  
  tempR <- tempfile(tmpdir = ".", fileext = ".R")
  on.exit(unlink(tempR))
  knitr::purl(file, output=tempR, quiet = TRUE)
  
  envir <- globalenv()
  source(tempR, local = envir, ...)
}