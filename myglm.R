binomial_abc <- function(linkfun){
  a   <- function(y){return(y)}
  db  <- function(theta){return(1/(theta*(1-theta)))}
  ddb <- function(theta){return((-1+2*theta)/(theta^2*(1-theta)^2))}
  dc  <- function(theta){return(-1/(1-theta))} #we assumed n_i=1
  ddc <- function(theta){return(-1/(1-theta)^2)}#we assumed n_i=1
  invlink <- function(linkfun){
    if (linkfun=="logit"){
      fun <- function(nu){return(1/(1+exp(-nu)))}
      return(fun)
    }
    else if (linkfun=="probit"){
      fun <- function(nu){return(pnorm(nu))}
      return(fun)
    }
  }
  dinvlink <- function(linkfun){
    if (linkfun=="logit"){
      fun <- function(nu){return(exp(-nu)/(1+exp(-nu))^2)}
      return(fun)
    }
    else if (linkfun=="probit"){
      fun <- function(nu){return(dnorm(nu))}
      return(fun)
    }
  }
  loglike <- function(Y,theta){
    return(Y*log(theta/(1-theta))+log(1-theta))
  }#we assumed n_i=1
  return(list("a"=a,"db"=db,"ddb"=ddb,"dc"=dc,"ddc"=ddc,"loglike"=loglike,"invlink"=invlink(linkfun),"dinvlink"=dinvlink(linkfun)))
}
gaussian_abc <- function(linkfun,sigma=1.0){
  a   <- function(y){return(y)}
  db  <- function(theta){return(matrix(1.0,length(theta))/sigma^2)}
  ddb <- function(theta){return(matrix(0.0,length(theta)))}
  dc  <- function(theta){return(-theta/sigma^2)} #we assumed n_i=1
  ddc <- function(theta){return(-matrix(1.0,length(theta))/sigma^2)}#we assumed n_i=1
  invlink <- function(linkfun){
    if (linkfun=="identity"){
      fun <- function(nu){return(nu)}
      return(fun)
    }
  }
  dinvlink <- function(linkfun){
    if (linkfun=="identity"){
      fun <- function(nu){return(matrix(1.0,length(nu)))}
      return(fun)
    }
  }
  loglike <- function(Y,theta,sigmaint=sigma){
    return(log(1/(sqrt(2*pi*(sigmaint^2)))) - 0.5*((Y- theta)^2)/(sigmaint^2))
  }#
  return(list("a"=a,"db"=db,"ddb"=ddb,"dc"=dc,"ddc"=ddc,"loglike"=loglike,"invlink"=invlink(linkfun),"dinvlink"=dinvlink(linkfun)))
}


poisson_abc <- function(linkfun){
  a   <- function(y){return(y)}
  db  <- function(theta){return(1/theta)}
  ddb <- function(theta){return(-1/((theta)^2))}
  dc  <- function(theta){return(-1)} 
  ddc <- function(theta){return(0)}
  invlink <- function(linkfun){
    if (linkfun=="log"){
      fun <- function(nu){return(exp(nu))}
      
      return(fun)
    }
  }
  dinvlink <- function(linkfun){
    if (linkfun=="log"){
      fun <- function(nu){return(exp(nu))}
      return(fun)
    }
  }
  
  loglike <- function(Y,theta){
    return(Y * log(theta) - theta - lgamma(Y + 1))
  }
  return(list("a"=a,"db"=db,"ddb"=ddb,"dc"=dc,"ddc"=ddc,"loglike"=loglike,"invlink"=invlink(linkfun),"dinvlink"=dinvlink(linkfun)))
}

negbinomial_abc <- function(linkfun, r=1.0){
  a   <- function(y){return(y)}
  db  <- function(theta){return(1/(theta*(theta/r+1)))}
  ddb <- function(theta){
    V = theta + theta^2/r
    dV = 1+2*theta/r
    return(-dV/(V^2))}
  dc  <- function(theta){return(-theta * db(theta))} 
  ddc <- function(theta){return(-(db(theta) + theta * ddb(theta)))}
  
  invlink <- function(linkfun){
    if (linkfun=="log"){
      fun <- function(nu){return(exp(nu))}
      
      return(fun)
    }
  }
  dinvlink <- function(linkfun){
    if (linkfun=="log"){
      fun <- function(nu){return(exp(nu))}
      return(fun)
    }
  }
  
  loglike <- function(Y,theta, r0=r){
    return(Y*log(theta/(theta+r0)) + (r0*log(r0/(theta +r0))) + lgamma(Y+r0)-lgamma(r0)- lgamma(Y+1))
  }
  return(list("a"=a,"db"=db,"ddb"=ddb,"dc"=dc,"ddc"=ddc,"loglike"=loglike,"invlink"=invlink(linkfun),"dinvlink"=dinvlink(linkfun)))
}

compute_abc_func <- function(family,linkfun,other=1.0){
  if (family=="binomial")
  {
    return(binomial_abc(linkfun))
  }
  else if (family=="gaussian")
  {
    return(gaussian_abc(linkfun,other))
  }
  else if (family=="poisson")
  {
    return(poisson_abc(linkfun))
  }
  else if (family=="negbinomial")
  {
    return(negbinomial_abc(linkfun,other))
  }
}

# mu and variance
mu <- function(theta,funs){return(-funs$dc(theta)
                                  /funs$db(theta))}
vary <- function(theta,funs){
  return(( funs$ddb(theta)*funs$dc(theta)
           -funs$ddc(theta)*funs$db(theta))
         /(funs$db(theta))^3)}

Define_Gradient_Hessian <- function(X,Y,funs){
  Gradient <- function(beta){
    theta = funs$invlink(X%*%beta)
    w = (funs$a(Y)-mu(theta,funs))/vary(theta,funs)*funs$dinvlink(X%*%beta)
    return(-t(X)%*%w)
  }
  Hessian <- function(beta){
    theta = funs$invlink(X%*%beta)
    d = (-(1/vary(theta,funs))*(funs$dinvlink(X%*%beta))^2)[,1]
    D=diag(d)#d is a list
    return(-t(X)%*%D%*%X)
  }
  return (list("Gradient"=Gradient, "Hessian"=Hessian))
}

# Newton-Raphson
newton.raphson <- function(Gradient, Hessian, x0 , tol = 1e-5, n = 5000) {
  # x0 start value, we called it theta_0
  # Gradient: function that returns gradient (that is a vector)
  # Hessian: function that returns Hessian matrix
  for (i in 1:n) {
    
    x1 <- x0 - solve( Hessian(x0) , Gradient(x0)) # Calculate next value x1
    #stopping criterion
    if (norm(x1 - x0, type="2") < tol) {
      return(x1)
    }
    # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
    #print(i)
  }
  print('Too many iterations in method')
}

deviance_null <- function(Y,family,linkfun){
  X=matrix(1,length(Y))
  out=myglm(X,Y,family,linkfun)
  return(out$deviance)
}

#This function can be improved by adding a stopping criterion
alteroptim <- function(Y, X, family, linkfun, hat_beta, tol = 1e-5, n = 200) {
  beta = hat_beta
  other_mle_old = 1  
  
  for (iter in 1:n) {
    beta_old = beta
    if (iter == 1) {
      other_mle_old = 1  # Initial value for first iteration
    } else {
      other_mle_old = other_mle  # Previous iteration's value
    }  

    funs = compute_abc_func(family, linkfun, 1) 
    
    # Negative loglikelihood to be minimised
    fr <- function(r) { 
      mu = funs$invlink(X %*% beta)
      return(-sum(funs$loglike(Y, mu, r)))
    }
    other_mle = optim(c(1), fr, method = "Brent", lower = 1e-6, upper = 1e6)$par
    
    # Step 2: Update beta coefficients
    funs = compute_abc_func(family, linkfun, other_mle) 
    out = Define_Gradient_Hessian(X, Y, funs)
    Gradient = out$Gradient
    Hessian = out$Hessian
    
    beta = newton.raphson(Gradient, Hessian, beta, tol = 1e-5, n = 1000)
    
    beta_change = max(abs(beta - beta_old))
    disp_change = abs(other_mle - other_mle_old)
    
    
    if (beta_change < tol && disp_change < tol) {
      #cat("Converged at iteration", iter, "\n")
      return(list("beta" = beta, "other_mle" = other_mle))
    }
  }
  print("Failed to converge within ", n, " iterations")
  return(list("beta" = beta, "other_mle" = other_mle))
}

# Output
compute_outputs = function(family,linkfun,hat_beta,Y,X){
  if (family=="gaussian" |family=="negbinomial"){
    out=alteroptim(Y,X,family,linkfun, hat_beta, tol = 1e-5, n = 200)
    hat_beta=out$beta
    other_mle=out$other_mle
    
    if (family=="gaussian"){
      other_par = 1
    }
    else
    {other_par=other_mle}
    funs = compute_abc_func(family,linkfun,other_par) 
    #dispersion residuals
    theta = funs$invlink(X%*%hat_beta)
    dispersion = sum((funs$a(Y)-mu(theta,funs))^2/vary(theta,funs))/(length(Y)-length(hat_beta))
    # recompute functions
    funs = compute_abc_func(family,linkfun,other_par)
    out = Define_Gradient_Hessian(X,Y,funs)
    Hessian  = out$Hessian
    # standard deviation of betas
    IH = solve(Hessian(hat_beta)) # inverse Hessian
    stdv = sqrt(dispersion)*sqrt(diag(IH))
    #t value 
    z=hat_beta/stdv
    #p-value using t-student
    pvalz = (1-pt(abs(z),length(Y)-length(hat_beta)))*2
    #Deviance 
    deviance = sum(2*funs$loglike(Y,Y,other_par)-2*funs$loglike(Y,funs$invlink(X%*%hat_beta),other_par))
    # Deviance Null-model
    dev_null = NaN
    if (!all(X == 1)){
      dev_null = sum(deviance_null(Y,family, linkfun))
    }
    #AIC: it uses the log-likelihood of  MLE
    aic = 2*(length(hat_beta)+1)-2*sum(funs$loglike(Y,funs$invlink(X%*%hat_beta),
                                                    other_mle))
    Ypred = funs$invlink(X%*%hat_beta)
  d = 2*funs$loglike(Y,Y,other_par) - 2*funs$loglike(Y,Ypred,other_par)
  deviance_residuals = sqrt(d)*sign(Y-Ypred)    
      }
  
  if (family=="binomial" |family=="poisson"){
    other_mle = NaN
    funs = compute_abc_func(family,linkfun)
    out = Define_Gradient_Hessian(X,Y,funs)
    Hessian  = out$Hessian
    
    # standard deviation
    IH = solve(Hessian(hat_beta)) # inverse Hessian
    stdv = sqrt(diag(IH))
    #AIC
    aic = 2*length(hat_beta)-sum(2*funs$loglike(Y,funs$invlink(X%*%hat_beta)))
    #z 
    z=hat_beta/stdv
    #p-value using the Normal dsitribution
    pvalz = (1-pnorm(abs(z)))*2
    #Deviance 
    deviance = sum(2*funs$loglike(Y,Y)-2*funs$loglike(Y,funs$invlink(X%*%hat_beta)))
    # Deviance Null-model
    dev_null = NaN
    if (!all(X == 1)){
      dev_null = sum(deviance_null(Y,family, linkfun))
    }
    dispersion=1
    
    Ypred = funs$invlink(X%*%hat_beta)
    d = 2*funs$loglike(Y,Y) - 2*funs$loglike(Y,Ypred)
    deviance_residuals = sqrt(d)*sign(Y-Ypred)
  }

  
  predict <- function(Xpred,type= "response"){
    mu = funs$invlink(Xpred%*%hat_beta)
    if(type == "response")
    {return(mu)}
  }
  return(list("aic"=aic,"z"=z,"stdv"=stdv,"pvalz"=pvalz,"deviance"=deviance,"dev_null"=dev_null,"hat_beta"=hat_beta,"other_parameter_mle"=other_mle,"dispersion"=dispersion,"predict" = predict,"deviance_residuals"=deviance_residuals))
  
}


# GLM function
myglm <- function (X,Y,family,linkfun){
  funs = compute_abc_func(family,linkfun)
  
  # Gradient and Hessian
  out = Define_Gradient_Hessian(X,Y,funs)
  Gradient = out$Gradient
  Hessian  = out$Hessian
  
  #initial beta
  beta0 = matrix(0.1,ncol(X))
  #MLE
  hat_beta = newton.raphson(Gradient, Hessian, beta0 , tol = 1e-5, n = 1000)
  
  
  
  outputs = compute_outputs(family,linkfun,hat_beta,Y,X)
  
  #print("Warning: the deviance of the null model is computed using a model including  the intercept only.\n
        #This is different from what the standard GLM does.")
  return(list("estimate"=outputs$hat_beta, 
              "stdError"=outputs$stdv,
              "z"=outputs$z,
              "pvalz"=outputs$pvalz,
              "other_parameter_mle"=outputs$other_parameter_mle,
              "dispersion"=outputs$dispersion,
              "deviance"=outputs$deviance,
              "deviance_null"=outputs$dev_null,
              "AIC"=outputs$aic,
              "predict" = outputs$predict,
              "deviance_residuals" = outputs$deviance_residuals))
}

