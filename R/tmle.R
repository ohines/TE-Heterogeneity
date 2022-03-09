#Function for computing ATE and VTE using TMLE
#
#data is a dataframe or tibble containing:
#Y - outcome of interest
#A - treatment indicator
#pi_hat - estimate of the propensity score E(A|X)
#mu1_hat - estimate of E(Y|A=1,X)
#mu0_hat - estimate of E(Y|A=0,X)
#
#ab - (optional) Min and Max of range of Y. Defaults to min/max of (Y,mu0_hat,mu1_hat)

TMLE_VTE <- function(data,ab=NULL){
  ## Some params
  max.it <- 600 #maximum number of iterations in targetting step
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
    a = with(data,min(Y,mu1_hat,mu0_hat))
    b = with(data,max(Y,mu1_hat,mu0_hat))
  }
  
  N = NROW(data)
  A <- data$A
  beta1 <- with(data,  1/pi_hat )
  beta0 <- with(data,  1/(pi_hat-1) )
  beta <- ifelse(A,beta1,beta0)
  y <- (data$Y -a)/(b-a)

  expitq1 <- qlogis((data$mu1_hat-a)/(b-a))
  expitq0 <- qlogis((data$mu0_hat-a)/(b-a))
  
  retarget <- TRUE
  it <- 1
  while(retarget & it <= max.it){
    q_1 <- plogis(expitq1)
    q_0 <- plogis(expitq0)
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
    
    if( (N^2*pD1^2 >= Sig1) | (N^2*pD2^2 >= Sig2)){
      expitq1 <- expitq1  + eps*beta1*(pD1+2*(cate-ATE)*pD2)/pD_norm
      expitq0 <- expitq0  + eps*beta0*(pD1+2*(cate-ATE)*pD2)/pD_norm 
      
    }else{
      retarget <- FALSE
    }
    it <- it+1
  }
  if(it>=max.it) warning("Max iterations reached in TMLE")
  
  rootV <- ifelse(VTE>=0,sqrt(VTE),NA)
  ss <- sqrt(Sig2/N)

  #out <- matrix(
  #  c(ATE*(b-a),sqrt(Sig1/N)*(b-a),
  #    rootV*(b-a),ss*(b-a)/rootV,
  #    VTE*(b-a)^2,ss*(b-a)^2),
  #  ncol=2,byrow=TRUE)
  #rownames(out) <- c("ATE","rootVTE","VTE")
  #colnames(out) <- c("Estimate","Std_err")
  
  coef <- c(ATE*(b-a),rootV*(b-a),VTE*(b-a)^2)
  std.err <- c(sqrt(Sig1/N)*(b-a),ss*(b-a)/rootV,ss*(b-a)^2)
  names(coef) <- names(std.err) <- c("ATE","rootVTE","VTE")
  names(cate) <- NULL
  
  out<- list(
    coef = coef,
    std.err = std.err,
    CATE = cate*(b-a)
  )
  
  return(out)
}

