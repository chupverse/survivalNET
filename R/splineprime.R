
splineprime <- function(time, gamma, m, mpos = NULL)
{ #nombre de noeuds internes
  if(length(gamma)!=(m+2))(stop("The length of the gamma's coefficients vector
                                should be ", m+2))
  #position des quantiles
  if(!is.null(mpos) & length(mpos) != m+2)
    (stop("Number of knots positions must be equal to ", m+2, " with ", m,
          " internal knots & 2 boundary knots"))
  if(is.null(mpos) == TRUE){a = c()
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
    spln <- gamma[2]
    list <- list(
      time = time, 
      spln = spln, 
      mpos = mpos)
  }
  else{
    phi <- c()
    nu_prime <- c()
    spln <- 0
    for(i in 2:(length(gamma)-1)){
      phi <- c(phi, (mpos[m+2]-mpos[i])/(mpos[m+2]-mpos[1]))
      nu_prime <- cbind(nu_prime,
                      3*pmax(0,(time-mpos[i]))^2-3*phi[i-1]*pmax(0,(time-mpos[1]))^2
                      -3*(1-phi[i-1])*pmax(0,(time-mpos[m+2]))^2)
      spln <- spln + gamma[i+1]*nu_prime[,i-1] }
    spln <- spln + gamma[2]
    list <- list(
      time = time,
      spln = spln,
      mpos = mpos,
      phi = phi,
      nu_prime = nu_prime)
  }
  return(list)
  
}

