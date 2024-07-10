
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
      r <- stats::integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol,
                            abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
      if ( !(r[['message']] %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) { stop(r$message) }  }
      return(r)
    }
    
    return(integrateA(Vectorize(.f), lower=0, upper=time, subdivisions = subdivisions)$value)
  }
  ### en développement
  # if(method == "table"){
  #   # year_seq = seq(0,t, by = 365.24)
  #   # birthday_seq = sapply(1:nrow(dataK), function(i)
  #   #                     {tail(seq.Date(as.Date(dataK$year[i]-dataK$age[i],
  #   #                     origin = "1960-1-1"),as.Date(t, origin =
  #   #                     as.Date(dataK$year[i])),"years"),floor(t/365.24))
  #   # }
  #   # )
  # 
  #   birth_date <- format(as.Date(dataK$year[2]-dataK$age[2],
  #                 origin = "1960-1-1"), "%m-%d")
  #   # year_date <- format(as.Date(dataK$year[2], origin = "1960-1-1"), "%m-%d")
  # 
  #   this_birthday <- as.Date(paste0(as.numeric(format(as.Date(dataK$year[2],
  #                                 origin = "1960-1-1"), "%Y")),
  #                                 paste0("-",birth_date)))
  # 
  #   next_year_date <- as.Date(paste0(as.numeric(format(as.Date(dataK$year[2],
  #                                 origin = "1960-1-1"),"%Y")) + 1, "-01-01"))
  # 
  #   if(this_birthday>as.Date(dataK$year[2], origin = "1960-1-1")){
  #     as.numeric(difftime(next_year_date, this_birthday, units = "days"))}
  #   if(this_birthday<as.Date(dataK$year[2], origin = "1960-1-1")){}
  # 
  # 
  # }
}