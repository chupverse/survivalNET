
# Expected cummulative hazard, age in days

expectedcumhaz <- function(ratetable, age, year, sex, time, method="exact", subdivisions = 100)
{
  if(method=="exact") {
    .year <- as.numeric(format( as.Date((1:time) + year, origin = "1960-01-01"), "%Y" ) ) 
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
  
  if (method == "table") {
    
    .yrs   <- 1800:2600
    .jan1  <- as.numeric(as.Date(paste0(.yrs, "-01-01")) - as.Date("1960-01-01"))
    .year_of <- function(d) .yrs[findInterval(d, .jan1)]
    
    table_scalar <- function(ratetable, age, year, sex, t) {
      if (t <= 0) return(0)
      
      max_age  <- max(as.numeric(dimnames(ratetable)$age))
      max_year <- max(as.numeric(dimnames(ratetable)$year))
      
      lo <- floor(age / 365.24) + 1
      hi <- floor((age + t) / 365.24)
      b_age <- if (hi >= lo) (lo:hi) * 365.24 - age else numeric(0)
      
      ## ruptures d'année : 1er janvier du calendrier
      b_year <- .jan1[.jan1 > year & .jan1 < year + t] - year
      
      cuts  <- sort(unique(c(0, b_age, b_year, t)))
      delta <- diff(cuts)
      mid   <- cuts[-length(cuts)] + delta / 2    
      
      a_i <- pmin(floor((age + mid) / 365.24), max_age)
      y_i <- pmin(.year_of(year + mid), max_year)
      
      values <- list(age  = as.character(a_i),
                     year = as.character(y_i),
                     sex  = rep(sex, length(a_i)))
      idx <- do.call(cbind, values[names(dimnames(ratetable))])   
      
      sum(ratetable[idx] * delta)
    }
    
    results <- vapply(time, function(t) table_scalar(ratetable, age, year, sex, t), numeric(1))
    
    return(matrix(results, ncol =1))
  }
  
}


