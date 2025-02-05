
splinecube <- function(time, gamma, m, mpos = NULL)
  
{
  
 ##error check
  if(length(time)<2)stop("length(time)<2. At least 2 times are necessary 
                                 to compute the spline")
  
  if(length(gamma)!=(m+2))stop("The length of the gamma's coefficients vector
                                should be ", m+2)
  if(!is.null(mpos) & m != length(mpos))stop("number of internal knots
                                positions must be equal to m=", m)

  ##
  
  x <- log(time)

  
  if(is.null(mpos) == TRUE){
      a <- c()
      for(i in (0:(m+1))){
          a <- c(a,i/(m+1))}
          mpos <- quantile(x, probs = a)
          mpos <- as.numeric(mpos)}
      else{
          a <- c(0,mpos,1)
          mpos <- quantile(x, probs = a)
  }
  
  if(m==0){
      spln <- gamma[1]+gamma[2]*x
      
      res <- list(
        spln = spln,
        mpos = mpos)
  }
  else{
      phi <- c()
      nu <- c()
      spln <- 0
      for(i in 2:(length(gamma)-1)){
        phi <- c(phi, (mpos[m+2]-mpos[i])/(mpos[m+2]-mpos[1]))
        nu <- cbind(nu,pmax(0,(x-mpos[i]))^3-phi[i-1]*pmax(0,(x-mpos[1]))^3-(1-phi[i-1])*pmax(0,(x-mpos[m+2]))^3)
        spln <- spln + gamma[i+1]*nu[,i-1] }
      spln <- spln + gamma[1]+gamma[2]*x
      
      res = list(
        spln = spln,
        knots = mpos,
        phi = phi,
        nu = nu)
  }
  return(res)
}