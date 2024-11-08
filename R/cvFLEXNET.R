
cvFLEXNET <- function(formula, pro.time=NULL, data, ratetable, cv=10, 
                             m = 2, mpos = NULL){
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a ratetable argument is required")
  if(!is.null(m)){
    if(length(m) != 1 && !is.null(mpos)) stop(
  "The cross-validation function can't handle multiple 'm' number of knots and multiple 'mpos' knots positions. Please input either multiple 'm' and 'mpos = NULL' or only one 'm' and multiple knots positions.")
  }
  if( !is.null(mpos) ){
  if(!is.list(mpos))stop("'mpos' must be a list.")
  }
  if(length(m) == 1 && !is.null(mpos) && length(unique(sapply(mpos, length))) != 1 )stop(
    "The different 'mpos' positions must be of the same length (m=",m,")"
  )
  if(is.null(m) && is.null(mpos))stop("No list of hyper parameters ('m' or 'mpos') given to do cross-validation.")
  #######

  if(is.null(m) && length(unique(sapply(mpos, length))) == 1){ 
    m <- unique(sapply(mpos, length))+2
    }

  times <- as.character(formula[[2]][2])
  failures <- as.character(formula[[2]][3])

  all_terms <- attr(terms(formula), "term.labels")
  group_term <- grep("group\\(", all_terms, value = TRUE)
  ratetable_terms <- grep("ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  if(length(ratetable_terms)>1) stop("More than one 'ratetable' term found in  the formula.")
  CV <- setdiff(all_terms, c(group_term, ratetable_terms))
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
  
  ##diff quanti/quali
  quali_col <- c()
  quanti_col <- c()
  warn <- 0
  col_warn <- c()

  for (col in CV) {

    unique_values <- unique(as.data.frame(cova)[[col]])

    if (length(unique_values) == 2 && !all(unique_values %in% c(0, 1))) {
      warn <- warn + 1
      col_warn <- c(col_warn, col)
    }

    if (all(unique_values %in% c(0, 1))) {
      quali_col <- c(quali_col, col)
    } else if (length(unique_values) > 2) {
      quanti_col <- c(quanti_col, col)
    }
  }
  if (warn > 0) {
    warning(paste(warn, "columns have exactly 2 modalities but are not 0 and 1. (",col_warn,"). Those
                  columns have been considered as quantitative variables."))
  }

  cov.quali <- quali_col
  cov.quanti <- quanti_col


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
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
  age <- ratetable_vars$age
  year <- ratetable_vars$year
  sexchara <- ratetable_vars$sex
    data.net <- data[,c(times, failures, age, year, sexchara, group, CV)]

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

  if(is.null(mpos)){
    .grid <-  expand.grid(m = m)}else{
    .grid <- expand.grid(m = m, mpos  = mpos)
  }

  .CVtune<-vector("list",cv*dim(.grid)[1])

  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .CVtune[[l]]<-list(train=data.net[data.net$folds!=k, ], valid=data.net[data.net$folds==k, ], grid=.grid[j,])
      l=l+1
    }
  }

  net.time.par<-function(xx, times, failures, ratetable, age, year, 
                         sexchara,  group, CV, newtimes){

    if(length(xx$grid) == 1){
      m = xx$grid  
      mpos = NULL
    }else{
      m = xx$grid$m
      mpos = unlist(xx$grid$mpos)
    }
    data=xx$train
    newdata=xx$valid

    if(!(is.null(group))){
      .data <- data[,c(times, failures, age, year, sexchara, group, CV)]}   else{
        .data <- data[,c(times, failures, age, year, sexchara, CV)] }

    if(!is.null(group)){
      .outcome <- paste("Surv(", times, ",", failures, ")")
      .f <- as.formula( paste(.outcome, "~", paste( CV,  collapse = " + "), "+", group) )
    }else{
      .f <- formula
    }

    .net <- survivalFLEXNET(.f, data = .data,
                    m = m, mpos = mpos, ratetable = ratetable)

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

    return(as.matrix(.pred))
  }

  .preFIT<-list()
  .preFIT<-lapply(.CVtune, net.time.par, times=times, failures=failures, 
                  ratetable = ratetable, age= age, year = year, sexchara = sexchara,
                  group=group, CV = CV, newtimes=.time)

  .FitCV <- replicate(dim(.grid)[1], matrix(NA, nrow = length(data[,times]),
                                            ncol = length(.time)), simplify=F)
  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .FitCV[[j]][data.net$folds==k,] <- .preFIT[[l]]
      l<-l+1
    }
  }

  net.best.measure <- function(prediction.matrix, formula, data, prediction.times){
    .times <- as.character(formula[[2]][2])
    .failures <- as.character(formula[[2]][3])
    .outcome <- paste("Surv(", .times, ",", .failures, ")")
    .predformula <- as.formula(paste(.outcome, "~ 1"))
    return(metrics(formula = .predformula, prediction.matrix =
                    as.matrix(as.data.frame(prediction.matrix)),
                   data=data, prediction.times=prediction.times,
                   pro.time=pro.time, metric="ci"))
  }

  .measure<-sapply(.FitCV, net.best.measure, formula = formula , data=data.net, prediction.times=.time)

  if(is.null(mpos)){
    .res <- data.frame(m = .grid[,1], measure = .measure)
    }else{
      
    .res <- data.frame(m = .grid[,1], mpos = as.character(.grid[,2]), measure = .measure)
  }
  .maxi<-.res[which(.res$measure==max(.res$measure, na.rm=TRUE) & is.na(.res$measure)==FALSE),]
  .maxi<-.maxi[1,]
  maxi_mpos <- as.numeric(unlist(strsplit(sub("^[^\\(]+\\((.*)\\)$", "\\1",
                                              as.character(.maxi$mpos)), "," ))) 
  return( list(optimal=list(m=.maxi$m,
                            mpos = maxi_mpos
                            ),
               results=.res ))
}

