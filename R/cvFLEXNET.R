
cvFLEXNET <- function(formula, pro.time=NULL, data, ratetable, cv=10, 
                       m = 2, mpos = NULL, mquant = NULL, init = NULL, delta_th = 0,
                       weights = NULL, m_s = NULL, Kref = NULL, stand = NULL){
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a ratetable argument is required")
  if(!is.null(m)){
    if(length(m) != 1 && !is.null(mpos)) stop(
      "The cross-validation function can't handle multiple 'm' number of knots and multiple 'mpos' knots positions. Please input either multiple 'm' and 'mpos = NULL' or only one 'm' and multiple knots positions.")
  }
  if(!is.null(m)){
    if(length(m) != 1 && !is.null(mquant)) stop(
      "The cross-validation function can't handle multiple 'm' number of knots and multiple 'mquant' knots quantile positions. Please input either multiple 'm' and 'mquant = NULL' or only one 'm' and multiple knots positions.")
  }
  if(!is.null(mpos) & !is.null(mquant))warning("'mpos' and 'mquant' have both been specified. 'mpos' values have been chosen over 'mquant' quantiles.") #ligne : mquant = NULL
  if(!is.null(mpos) & !is.null(mquant)){mquant = NULL} #ligne : mquant = NULL
  
  if( !is.null(mpos) ){
    if(!is.list(mpos))stop("'mpos' must be a list.")
  }
  
  if( !is.null(mquant) ){
    if(!is.list(mquant))stop("'mquant' must be a list.")
  }
  if(!is.null(mpos)){
    if(all(sapply(mpos, length) != length(mpos[[1]])))stop("The different knots values must have the same length across the list.")
  }
  if(!is.null(mquant)){
    if(all(sapply(mquant, length) != length(mquant[[1]])))stop("The different knots quantile values must have the same length across the list.")
  }
  if(length(m) == 1 && !is.null(mpos) && unique(sapply(mpos, length)) != m+2 )stop(
    "The different 'mpos' values must be of the same length as the number of internal knots+2 (m+2=",m+2,")"
  )
  if(length(m) == 1 && !is.null(mquant) && length(unique(sapply(mquant, length))) != m+2 )stop(
    "The different 'mquant' quantile positions must be of the same length as the number of internal knots +2 (m+2=",m+2,")"
  )
  if (is.null(m) && is.null(mpos) && is.null(mquant)) {
    stop("No list of hyperparameters ('m', 'mpos', or 'mquant') given to do cross-validation.")
  }  #######
  
  if(is.null(m) && length(unique(sapply(mpos, length))) == 1){ 
    m <- unique(sapply(mpos, length))
  }
  if(is.null(m) && length(unique(sapply(mquant, length))) == 1){ 
    m <- unique(sapply(mquant, length))
  }
  
  times <- as.character(formula[[2]][2])
  failures <- as.character(formula[[2]][3])
  
  all_terms <- attr(terms(formula), "term.labels")
  group_term <- grep("group\\(", all_terms, value = TRUE)
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms)>1) stop("More than one 'strata' term found in  the formula. Only one variable at a time can be stratified")
  ratetable_terms <- grep("^ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  if(length(ratetable_terms)>1) stop("More than one 'ratetable' term found in  the formula.")
  CV <- setdiff(all_terms, c(group_term, strata_terms, ratetable_terms))
  if(length(CV) == 0){covnames = "1"} else{covnames <- CV}
  
  cova <- as.matrix(data[,CV])
  for(i in colnames(cova)){ if(!is.numeric(cova[,i])) stop("All covariates must be numeric")  }
  
  extract_vars <- function(term) {
    var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
    vars <- trimws(unlist(strsplit(var_string, ",")))
    return(vars)
  }
  assign_ratetable_vars <- function(vars) {
    age <- year <- sex <- NULL
    for (var in vars) {
      if (grepl("age = ", var)) {
        age <- sub("age = ", "", var)
      } else if (grepl("year = ", var)) {
        year <- sub("year = ", "", var)
      } else if (grepl("sex = ", var)) {
        sex <- sub("sex = ", "", var)
      }
    }
    unnamed_vars <- setdiff(vars, c(age, sex, year))
    if (length(unnamed_vars) > 0) {
      if (is.null(age) && length(unnamed_vars) >= 1) age <- unnamed_vars[1]
      if (is.null(year) && length(unnamed_vars) >= 2) year <- unnamed_vars[2]
      if (is.null(sex) && length(unnamed_vars) >= 3) sex <- unnamed_vars[3]
    }
    return(list(age = age, year = year, sex = sex))
  }
  
  if(!is.null(init)){
    if(!is.list(init))stop("'init' must be a list of length  ", length(m))
    if(length(init) != length(m))stop("'init' must be a list of length  ", length(m))
  }
  
  if(length(group_term) == 0){
    group = NULL
  }else{
    group <- unlist(lapply(group_term, extract_vars))
    if(!all(unique(data[,group]) %in% c(0, 1))) stop("The ", group," covariate can only take the values 0 or 1.")
  }
  
  if(length(CV) == 0){
    CV = NULL
  }
  ####
  strata_var = unlist(lapply(strata_terms, extract_vars))
  if(!is.null(strata_var) && strata_var %in% covnames) stop("The stratified covariate also appears as a covariate in the formula.")
  if(is.null(strata_var)){
    timevar = strata_var
    xlevels = NULL
  }
  if(!is.null(strata_var)){
    if(is.null(m_s))stop("The number of internal knots 'm_s' for the stratifed splines needs to be specified")
    timevar <- data[,strata_var]
    xlevels <- list(levels(as.factor(timevar)))
    names(xlevels) <- c(strata_var)
  }
  if(!is.null(timevar)){
    if(length(unique(timevar))>15) stop("The variable with a time-dependant effect has too many categories (>15)")
  }
  
  #######
  
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
  age <- ratetable_vars$age
  year <- ratetable_vars$year
  sexchara <- ratetable_vars$sex
  data.net <- data[,c(times, failures, age, year, sexchara, group, strata_var, CV)]
  
  if(is.null(pro.time)) {pro.time <- median(data[,times])}
  
  sample_id <- sample(nrow(data.net))
  folds <- cut(seq(1,nrow(data.net)), breaks=cv, labels=FALSE)
  folds_id <- folds[sample_id]
  data.net$folds <- folds_id
  
  if(!is.null(group)){
    .outcome <- paste("Surv(", times, ",", failures, ")")
    .f <- as.formula( paste(.outcome, "~", paste( CV,  collapse = " + "), "+", group) )
  }else{
    .f <- formula
  }
  
  .time <- sort(unique(data.net[,times]))
  
  if(is.null(m_s)){
    if(is.null(mpos) && is.null(mquant)){
      .grid <-  expand.grid(m = m)}else if(!is.null(mpos)){
        .grid <- expand.grid(m = m, mpos  = mpos)}else{
          .grid <- expand.grid(m = m, mquant  = mquant)
        }
  }else{
    if(is.null(mpos) && is.null(mquant)){
      .grid <-  expand.grid(m = m, m_s = m_s)}else if(!is.null(mpos)){
        .grid <- expand.grid(m = m, mpos  = mpos, m_s = m_s)}else{
          .grid <- expand.grid(m = m, mquant  = mquant, m_s = m_s)
        }
  } 
  
  if(!is.null(init)){
    correstab <- data.frame(
      grid = .grid
    )
    correstab$init <- init
  }else{
    correstab <- data.frame(init <- rep(list(NULL), length(m)))
  }
  
  .CVtune<-vector("list",cv*dim(.grid)[1])
  
  if (!is.null(stand)) {
    if (!is.character(stand))
      stop("'standardize' must be a character vector of variable names.")
    missing_vars <- setdiff(stand, names(data.net))
    if (length(missing_vars) > 0)
      stop("The following variables are not in 'data': ",
           paste(missing_vars, collapse = ", "))
  }

  
  l<-1
  for (k in 1:cv){
    trainSt <- data.net[data.net$folds != k, ]
    validSt <- data.net[data.net$folds == k, ]
    if (!is.null(stand)) {
      muSt <- sapply(trainSt[, stand, drop = FALSE], mean)
      sigmaSt <- sapply(trainSt[, stand, drop = FALSE], sd)
      trainSt[, stand] <- scale(trainSt[, stand, drop = FALSE], center = muSt, scale = sigmaSt)
      validSt[, stand] <- scale(validSt[, stand, drop = FALSE],center = muSt,scale = sigmaSt)
    }
    for (j in 1:dim(.grid)[1]){
      .CVtune[[l]]<-list(train=trainSt, valid=validSt, grid=.grid[j,])
      l=l+1
    }
  }
  
  net.time.par<-function(xx, times, failures, ratetable, age, year, 
                         sexchara, group, strata_var, CV, newtimes){
    
    if(length(xx$grid) == 1){
      m = xx$grid  
      knots = NULL
      quant = NULL
    }else if ("mpos" %in% names(xx$grid) ){
      m = xx$grid$m
      knots = unlist(xx$grid$mpos)
      quant = NULL
    }else{
      m = xx$grid$m
      knots = NULL
      quant = unlist(xx$grid$mquant) }
    if(!is.null(m_s)){
      m_s <- xx$grid$m_s
    }
    init = unlist(xx$init)
    
    data=xx$train
    newdata=xx$valid
    
    if(!(is.null(group))){
      .data <- data[,c(times, failures, age, year, sexchara, group, CV)]}   else{
        .data <- data[,c(times, failures, age, year, sexchara, CV)] }
    
    if(!(is.null(strata_var))){
      .data <- cbind(.data, data[, strata_var])
      colnames(.data)[length(.data)] <- strata_var
    }else{ 
      .data <- .data}
    
    .net <- survivalFLEXNET(.f, data = .data,
                            m = m, mpos = knots, mquant = quant, ratetable = ratetable, init = init,
                            delta_th = delta_th, weights = weights, m_s = m_s, Kref = Kref)
    
    .time<-sort(unique(.data[,times]))
    
    .newdata <- data.frame(newdata[,c(group, CV)])
    .pred.temp <- predict(.net, newdata=newdata)
    .pred <- .pred.temp$predictions
    .time.net <- .pred.temp$times
    
    if(!is.null(newtimes)) {
      .pred.net <- cbind(rep(1, dim(.pred)[1]), .pred)
      .time.net <- c(0, .time.net)
      
      idx=findInterval(newtimes, .time.net)
      .pred=.pred.net[,pmax(1,idx)]
      
      .time <- newtimes
    }
    
    keep_mod <- .net
    
    return(list(
      pred = as.matrix(.pred),
      model = .net
    ))
  }
  
  .preFIT<-list()
  .preFIT<-lapply(.CVtune, net.time.par, times=times, failures=failures, 
                  ratetable = ratetable, age= age, year = year, sexchara = sexchara,
                  group=group, strata_var=strata_var, CV = CV, newtimes=.time)
  
  .FitCV <- replicate(dim(.grid)[1], matrix(NA, nrow = length(data[,times]),
                                            ncol = length(.time)), simplify=F)
  .modelsCV <- vector("list", dim(.grid)[1])
  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .FitCV[[j]][data.net$folds==k,] <- .preFIT[[l]]$pred
      .modelsCV[[j]][[k]] <- .preFIT[[l]]$model
      l<-l+1
    }
  }
  
  
  net.loglik <- function(model, data, formula, ratetable) {
    
    times <- as.character(formula[[2]][2])
    failures <- as.character(formula[[2]][3])
    
    time <- data[, times]
    event <- data[, failures]
    
    
    hP <- sapply(seq_along(time), function(i) {
      survivalNET::expectedhaz(ratetable = ratetable, age = data[i, age], sex = data[i, sexchara], year = data[i, year], time = time[i])
    })
    
    coeff <- model$coefficients
    
    m <- model$m  
    m_s <- model$m_s
    
    beta_est <- coeff[1:(length(coeff)-(m+2))]
    gamma_est <- tail(coeff, m+2)
    
    if(!is.null(m_s)){
      beta_est <-  coeff[(1:(dim(model$x)[2]) )] 
      gamma_est <- coeff[((dim(model$x)[2]+1): (length(coeff)))]  
    }
    
    cova <- as.matrix(data[, names(beta_est), drop = FALSE])
    
    
    if(is.null(m_s)){
      spln <- survivalNET::splinecube(time, gamma_est, m, model$mpos)$spln
      spln_P <- survivalNET::splinecubeP(time, gamma_est, m, model$mpos)$spln
      
      hE <- (1/time) * spln_P * exp(spln + cova %*% beta_est)
      
      loglik <- sum(event * log(hP + hE) - exp(spln + cova %*% beta_est))
      
    }else{
      
      value <- c()
      K <- sort(unique(data[,strata_var]))
      
      gamma_base <- gamma_est[1:(m+2)]
      
      for(k in K){
        betak <- beta_est
        gammak <- gamma_est[((m+2)+1+(k-1)*(m_s+2)):((m+2)+(k)*(m_s+2))]
        idx <- data[,strata_var] == k
        
        timek <- time[idx]
        eventk <- event[idx]
        hPk <- hP[idx]
        covak <- cova[idx, , drop = FALSE]
        
        splk_base <- survivalNET::splinecube(timek, gamma_base, m, model$mpos)$spln
        splkP_base <- survivalNET::splinecubeP(timek, gamma_base, m, model$mpos)$spln
        
        Kref <- model$Kref
        
        if (k != Kref) {
          splk  <- survivalNET::splinecube(timek, gammak, m_s, model$mpos_s[[k]])$spln
          splkP <- survivalNET::splinecubeP(timek, gammak, m_s, model$mpos_s[[k]])$spln
        } else {
          splk <- 0
          splkP <- 0
        }
        linpred <- splk_base + splk + covak %*% betak
        
        calc <- hPk + (1/timek) * (splkP_base + splkP) * exp(linpred)
        
        # value_k <- sum( (
        #   eventk[calc >= 0] * log( calc[calc >= 0] )) -
        #     exp(linpred[calc >= 0]))
        value_k <- sum( (
          eventk * log( calc )) -
            exp(linpred))
        value <- c(value, value_k)
        loglik <- sum(value)
      }
    }
    
    return(loglik)
  }  
  
  
  .measure <- matrix(NA, nrow = dim(.grid)[1], ncol = cv)
  
  for (j in 1:dim(.grid)[1]) {
    for (k in 1:cv) {
      
      model <- .modelsCV[[j]][[k]]
      
      valid_data <- .CVtune[[(k - 1) * nrow(.grid) + 1]]$valid
      
      .measure[j, k] <- net.loglik(
        model = model,
        data = valid_data,
        formula = formula,
        ratetable = ratetable
      )
    }
  }
  
  .measure = rowMeans(.measure)
  
  if(is.null(m_s)){
    if(is.null(mpos) && is.null(mquant)){
      .res <- data.frame(m = .grid[,1], measure = .measure)
    }else if(!is.null(mpos)){
      .res <- data.frame(m = .grid[,1], mpos = as.character(.grid[,2]), measure = .measure)
    }else{
      .res <- data.frame(m = .grid[,1], mquant = as.character(.grid[,2]), measure = .measure)
    }
  }else{
    if(is.null(mpos) && is.null(mquant)){
      .res <- data.frame(m = .grid[,1], m_s = .grid[,2], measure = .measure)
    }else if(!is.null(mpos)){
      .res <- data.frame(m = .grid[,1], mpos = as.character(.grid[,2]), m_s = .grid[,3], measure = .measure)
    }else{
      .res <- data.frame(m = .grid[,1], mquant = as.character(.grid[,2]), m_s = .grid[,3], measure = .measure)
    }
  }
  .maxi<-.res[which(.res$measure==max(.res$measure, na.rm=TRUE) & is.na(.res$measure)==FALSE),]
  .maxi<-.maxi[1,]
  maxi_mpos <- as.numeric(unlist(strsplit(sub("^[^\\(]+\\((.*)\\)$", "\\1",
                                              as.character(.maxi$mpos)), "," ))) 
  maxi_mquant <- as.numeric(unlist(strsplit(sub("^[^\\(]+\\((.*)\\)$", "\\1",
                                                as.character(.maxi$mquant)), "," ))) 
  
  res_list <- list(optimal=list(m=.maxi$m,
                                knots = maxi_mpos, quants = maxi_mquant
  ), results=.res )
  if(!(is.null(m_s))){
    res_list$optimal$m_s <- .maxi$m_s
  }
  return(res_list)
}