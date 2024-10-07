
# Expected cummulative hazard, age in days

expectedcumhaz <- function(ratetable, age, year, sex, time, method="exact", subdivisions = 100)
{
  if(method=="exact") {
    .year <- date.mdy(year+(1:time))$year
    .temp <- diag(
      ratetable[as.character( pmin( floor( (age + (1:time) )/365.24), max(as.numeric(names(ratetable[, "2000", "male"]))) ) ),
            as.character( pmin( .year, max(as.numeric(names(ratetable["51", , "male"]))) ) ),
            sex] )
    
    return(sum(.temp))
  }
  
  if(method=="trapezoidal") {
    .f <- function(x)   { expectedhaz(ratetable=ratetable, age=age, year=year, sex=sex, time=x)}
    
    integrateA <- function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25,
                           abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL)
    {
      r <- integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol,
                            abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
      if ( !(r[['message']] %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) { stop(r$message) }  }
      return(r)
    }
    
    return(integrateA(Vectorize(.f), lower=0, upper=time, subdivisions = subdivisions)$value)
  }
  
  ### under development
  if(method == "table"){

    birth_md <- format(as.Date(year-age,
                  origin = "1960-1-1"), "%m-%d")
    
    if(birth_md == "02-29"){birth_md <- "03-01"}
    
    bday <- as.Date(paste0(as.numeric(format(as.Date(year,
                    origin = "1960-01-01"), "%Y")),
                     paste0("-", birth_md)) ) 
    
    next_y <- as.Date(paste0(as.numeric(format(as.Date(year,
                              origin = "1960-01-01"),"%Y")) + 1, "-01-01") ) 
                      
    process_dates <- function(bday, next_y, end_date){
      bdays <- c()
      new_years <- c()
      y10 = as.numeric(difftime("1970-01-01","1960-01-01"))
      
      while (bday <= end_date | next_y <= end_date) {
        if (bday <= end_date) {
          bdays <- c(bdays, bday)
          bday <- bday + 31557600
        }
        if (next_y <= end_date) {
          new_years <- c(new_years, next_y)
          next_y <- next_y + 31557600
        }
      }
      
      return(list(birthdays = bdays+y10, new_years = new_years+y10))
    }
    
    results <- c()
    delta <- c()

    end_date <- as.Date(year + time)
      
    result <- process_dates(bday, next_y, end_date)

    results <- result
    results <- c(as.Date(year),as.Date(sort(c(results$birthdays,
                 results$new_years, end_date+ as.numeric(difftime("1970-01-01",
                 "1960-01-01")) )), origin = "1960-01-01"))
   
    delta <- as.numeric(difftime(results[-1], 
                 results[-length(results)], units = "days"))
    
    if( length(delta) %% 2 == 0){
      cond <- length(delta)/2 - 1 
      pair = TRUE
      }else{
      cond <- floor(length(delta)/2) 
      pair = FALSE}
    ##first part, if the anniversary is before the diag date
    if(bday < next_y ){
      haz_values <- c() 
   
      for (i in 0:cond ) {
        for (j in 0:1) {
         haz_values <- c(haz_values,
                       ratetable[as.character( pmin(floor(age/365.24 + i +j ),
                                max(as.numeric(names(ratetable[, "2000", "male"])))) ),
                                as.character( pmin(date.mdy(year)$year + i,
                                max(as.numeric(names(ratetable["51", , "male"])))) ),
                                sex])
       }
      }
      if (pair == FALSE){haz_values <- haz_values[-length(haz_values)]}  
   } 
    ## deuxieme partie de la fonction
    if(bday > next_y){
      haz_values <- c() 
      
      for (i in 0:cond) {
        for (j in 0:1) {
          haz_values <- c(haz_values,
                          ratetable[as.character( floor(age/365.24 + i ) ),
                                    as.character( date.mdy(year)$year + i+j) ,sex])
        }
      }
      if (pair == FALSE){haz_values <- haz_values[-length(haz_values)]}  
      
    }
    
    cumhaz <- haz_values%*%delta
    return(cumhaz)
  }
}


