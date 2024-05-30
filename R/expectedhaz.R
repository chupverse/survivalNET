
# Expected instantaneous hazard, age in days

expectedhaz <- function(ratetable, age, sex, year, time)
{
  age.year<-age/365.24 
  t.year<-time/365.24

  maxyear.ratetable<-max(as.numeric(attributes(ratetable)$dimnames[[3]]))
  minage.year.ratetable<-round(min(as.numeric(attributes(ratetable)$dimnames[[1]])/365.24))

  return(
      mapply(FUN = function(age, sex, year) { ratetable[age, sex, year] }, 
             trunc(trunc(age.year)+t.year)-minage.year.ratetable+1,
             sex, as.character(pmin(maxyear.ratetable,trunc(year+t.year)))
             )
  )
}

