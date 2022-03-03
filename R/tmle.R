#Function for computing ATE and VTE using TMLE
#
#df is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)
#
#ab - (optional) Min and Max of range of Y. Defaults to min/max of (Y,mu0_hat,mu1_hat)

TMLE_VTE <- function(df,ab=NULL){
  ## Some params
  max.it <- 20 #maximum number of iterations in targetting step
  eps <- 0.0001 #TMLE target step size
  
  
  if(!is.null(ab)){
    a = ab[1]
    b = ab[2]
    if(a>b) {
      warning("ab values incorrect - using default method")
      ab <- NULL
    }
  }
  if(is.null(ab)){
    a = with(df,min(Y,mu1_hat,mu0_hat))
    b = with(df,max(Y,mu1_hat,mu0_hat))
  }

  N = NROW(df)
  A <- df$A
  beta1 <- with(df,  1/pi_hat )
  beta0 <- with(df,  -1/(1-pi_hat) )
  beta <- ifelse(A,beta1,beta0)
  y <- (df$Y -a)/(b-a)
  q_1 <- (df$mu1_hat-a)/(b-a) #initial fit
  q_0 <- (df$mu0_hat-a)/(b-a)
  

  
  retarget <- TRUE
  it <- 1
  while(retarget & it <= max.it){
    cate <- q_1 - q_0
    #q_m   <- q_0 + A*cate
    #loss  <- -sum(y*log(q_m) + (1-y)*log(1-q_m))/N  
    po <- beta*(y-q_0 - A*cate) + cate
    
    a0 <- sum(cate)/N
    a1 <- sum(po)/N
    b0 <- sum(cate^2)/N
    b1 <- sum(po^2)/N
    b2 <- sum(po*cate)/N

    ATE <- a0
    VTE <- b0 - a0^2
    pD1 <- a1 - a0
    pD2 <- 2*(b2-b0+a0^2-a0*a1)
    pD_norm <- sqrt(pD1^2 + pD2^2)
    
    Sig1 <- b1 - 2*a1*a0 + a0^2
    Sig2 <- sum(((po-a0)^2 -(po-cate)^2 - VTE)^2  )/N 

    if( (N^2*pD1^2 >= Sig1) | (N^2*pD2^2 >= Sig1)){
      q_1 <- plogis(qlogis(q_1)  + eps*beta1*(pD1+2*(cate-ATE)*pD2)/pD_norm )
      q_0 <- plogis(qlogis(q_0)  + eps*beta0*(pD1+2*(cate-ATE)*pD2)/pD_norm ) 
      
    }else{
      retarget <- FALSE
    }
    it <- it+1
  }
  if(it>=max.it) warning("Max iterations reached in TMLE")
  
  rootV <- ifelse(VTE>=0,sqrt(VTE),NA)
  ss <- sqrt(Sig2/N)
  #return output
  out <- matrix(
    c(ATE*(b-a),sqrt(Sig1/N)*(b-a),
      rootV*(b-a),ss*(b-a)/rootV,
      VTE*(b-a)^2,ss*(b-a)^2),
      ncol=2,byrow=TRUE)
  rownames(out) <- c("ATE","rootVTE","VTE")
  colnames(out) <- c("Estimate","Std_err")
  return(out)
}

