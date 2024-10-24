
spline <- function(time, gamma, m, mpos = NULL)
{
  #number of internal knots
  if(length(gamma)!=(m+2))(stop("The length of the gamma's coefficients vector
            should be ", m+2))
  #quantiles positions
  if(!is.null(mpos) & length(mpos) != m+2)
    (stop("Number of knots positions must be equal to ", m+2, " with ", m, 
          " internal knots & 2 boundary knots"))
  if(is.null(mpos)){
    a <- c()
    for(i in (0:(m+1))){
      a <- c(a,i/(m+1))}
    mpos <- quantile(time, probs = a)
    mpos <- as.numeric(mpos)}
  else{
    a <- sort(mpos)
    mpos <- quantile(time, probs = a)
    mpos <- as.numeric(mpos)
  }
  
  if(m==0){
    spln <- gamma[1]+gamma[2]*time
    list <- list(
      time = time,
      spln = spln,
      mpos = mpos)
  }
  else{
    phi <- c()
    nu <- c()
    spln <- 0
    for(i in 2:(length(gamma)-1)){
      phi <- c(phi, (mpos[m+2]-mpos[i])/(mpos[m+2]-mpos[1]))
      nu <- cbind(nu,pmax(0,(time-mpos[i]))^3-phi[i-1]*pmax(0,(time-mpos[1]))^3
                 -(1-phi[i-1])*pmax(0,(time-mpos[m+2]))^3)
      spln <- spln + gamma[i+1]*nu[,i-1] }
    spln <- spln + gamma[1]+gamma[2]*time
    list <- list(time = time,
                 spln = spln,
                 mpos = mpos,
                 phi = phi,
                 nu = nu)
  }
  return(list)
  }



