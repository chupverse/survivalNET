
# Expected instantaneous hazard, age in days
# date in number of days since 1960 / age in years / time in days

expectedhaz <- function(ratetable, age, year, sex, time, max_age = NULL, max_year = NULL) 
{
  time <- min(time, 1000000)
  .year <- as.numeric(format( as.Date(time + year, origin = "1960-01-01"), "%Y" ) ) 
  
  if(is.null(max_age)){
  max_age <- max(as.numeric(dimnames(ratetable)[[1]]))
  }
  
  if(is.null(max_year)){
  max_year <- max(as.numeric(dimnames(ratetable)[[2]]))
  }
  
  ratetable[as.character( min( floor((age+time)/365.24), max_age ) ),
            as.character( min( .year, max_year ) ),
            sex]
}