
# Expected cummulative hazard, age in days

expectedcumhaz <- function(ratetable, age, sex, year, time)
{
  age.year<-age/365.24
  t.year<-time/365.24

  sumcum<-rep(0,length(age))
  j<-0   #time loop
  i<-1   #subject loop
  while (i<=length(age)){
    while (j<trunc(t.year[i])) {
      sumcum[i]<-sumcum[i]+365.24*expectedhaz(ratetable, age.year[i]*365.24, sex[i], year[i], j*365.24)
      j<-j+1
    } #end of time loop

  sumcum[i]<-sumcum[i]+(t.year[i]-trunc(t.year[i]))*expectedhaz(ratetable, age[i], sex[i], year[i],j)
  i<-i+1
  j<-0
  } #end of subject loop
return(sumcum)
}
