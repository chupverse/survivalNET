
# Expected instantaneous hazard, age in days
# date in number of days since 1960 / age in days / time in days

expectedhaz <- function(ratetable, age, year, sex, time, max_age = NULL, max_year = NULL) 
{
  time <- pmin(time, 1000000)
  .year <- as.numeric(format( as.Date(time + year, origin = "1960-01-01"), "%Y" ) ) 
  .age <- floor((age+time)/365.24)
  if(is.null(max_age)){
    max_age <- max(as.numeric(dimnames(ratetable)$age))
  }
  
  if(is.null(max_year)){
    max_year <- max(as.numeric(dimnames(ratetable)$year))
  }
  
  values <- list(
    age  = as.character(pmin(.age, max_age)),
    year = as.character(pmin(.year, max_year)),
    sex  = sex
  )
  
  idx <- do.call(cbind, values[names(dimnames(ratetable))])
  
  ratetable[idx]
}
