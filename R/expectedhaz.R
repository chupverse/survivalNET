
# Expected instantaneous hazard, age in days
# date in number of days since 1960 / age in years / time in days

expectedhaz <- function(ratetable, age, year, sex, time) 
{
  time <- min(time, 1000000)
  .year <- as.numeric(format( as.Date(time + year, origin = "1960-01-01"), "%Y" ) ) 
  
  ratetable[as.character( min( floor((age+time)/365.24), max(as.numeric(names(ratetable[, "2000", "male"]))) ) ),
            as.character( min( .year, max(as.numeric(names(ratetable["51", , "male"]))) ) ),
            sex]
}