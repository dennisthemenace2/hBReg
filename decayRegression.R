
### decay regression


#y is predicted by the features which have decreasing influence on the prediction.
# instead of estimating one beta for each feature be assume one beta that decreases over the features

#y =x1 * beta1*exp(-lambda*0) + x2 * beta1*exp(-lambda*1)+..+ beta0 

#M = number of features.
#t = 0:(M-1)
# the first number in t has to be 0 so that exp(0)=1
#y = x * beta1*exp(-lambda *t) + beta0

#Lets write down the joint probability 
#p(y,x,beta, lambda, s, sigma) = p(y|x*beta*exp(-lambda*t), sigma )
#                                p(beta|0,s) p(lambda|e0,f0) p(sigma|a0,b0) p(s|c0,c0)   
#
#Distributions used
# p(Y|..) = Normal 
# p(beta|..) = Normal
# p(sigma|..) = Gamma
# p(s|..) = Gamma
# p(lambda|..) = Gamma

####Generate data
set.seed(1234)
n <- 100
lambda = 1.125
beta   <- c(1.4,-0.5)
sigma2 <- 0.1
n_time = 10

beta_new = beta[1]%*% exp(-lambda* (0:(n_time-1))  )
beta_new = cbind(beta_new, beta[2])

X<-cbind(matrix(rnorm(n*n_time,2,1),ncol = n_time) ,1)

y<- matrix(rnorm(n,X%*% t(beta_new) ,sigma2) ,ncol = 1)


colnames(y) = 'y'
x<- matrix(X[,1:(ncol(X)-1) ],ncol = n_time)
colnames(x)= paste('x_',1:n_time,sep = '')



### note that the feature contribution is decreasing from left to right, that is important.
# if your data is arranged the other way around take care of this, or edit this code


###Variational approximation
bfreg <- setRefClass("bfreg",
                    fields = list(lambda ='numeric' ,beta1='numeric',beta0='numeric',sigma='numeric',s='numeric',bounds='numeric' ),
                    methods = list(

                  train = function(x,y,maxiterations = 1000,epsilon =1e-10,initialGuess=F ){

                    N = nrow(x)
                    M = ncol(x)
                    
                    ##need some inital guess
                    if(initialGuess == T){
                      model = lm('y~.',as.data.frame(cbind(x,y)))
                      betas = model$coefficients[2]
                      beta0 = model$coefficients[1]
                      lambda =exp( abs(model$coefficients[3]) )
                    }else{
                      lambda =rgamma(1,1,1)
                      betas = rnorm(1) 
                      beta0 = rnorm(1)
                    }
                    
                    ##set up the priors
                    a0  =0.001
                    an = a0+0.5*N
                    b0 = 0.001
                    bn = b0
                    c0 = 0.001
                    cn = c0
                    d0 = 0.001
                    dn = d0
                    e0 = 0.001
                    en = e0
                    f0 = 0.001
                    fn = f0
                    
                    cn = c0+0.5*2
                    
                    sn = c0/d0
                    sigma = a0/b0
                    
                    ###function to optimize p(lambda|..) \prop p(Y|x (beta,exp(-lambda*k) ) ) p(lambda|e0,f0)
                    fnlambda = function(lambda){
                      w = matrix( betas* exp(- lambda * 0:(M-1) ),ncol=1)
                      err = sum( ( (x%*%w +beta0) -y)^2)
                      res = -(sigma/2)*err  + (e0 -1)*log(lambda) - f0 *lambda
                      res
                    }
                    ###gradient function
                    fnlamdagr = function(lambda){
                      w = matrix( betas* exp(- lambda * 0:(M-1) ),ncol=1)
                      k = c(0:(M-1))
                      err = (x%*%w+beta0 -y)
                      res = 0
                      for(i in 1:N){
                        res = res + sum( err[i] *((x[i,]%*%(w*k ) )) )
                      }
                      res = sigma*res  + (e0 -1)/lambda - f0 
                      res
                    }
                    
                    lastBound = -Inf
                    for(i in 1:maxiterations){
                    
                        k  = c( exp(- lambda * 0:(M-1) ) )
                    
                        xnew  = x%*%k 
                        xnew = cbind(xnew,1)
                        
                        xky = t(xnew)%*%(y)
                        
                        xkxk = t(xnew)%*%(xnew)
                        dcov = (sigma*  xkxk )  + diag(  sn,2)
                        
                        Icov = solve(dcov)
                        
                        w = Icov %*% ( sigma *xky) 
                        
                        betas = w[1]
                        beta0 = w[2]
                        
                        err = sum( (xnew%*%w -y)^2 )
                              
                        #minimize lambda
                        lambdaoptimres = optim(par = c(lambda), fn= fnlambda,method = "L-BFGS-B",gr  = fnlamdagr,
                                             lower = 10e-10, upper = Inf,
                                             control = list(fnscale = -1), hessian = T)  
                        lambda = as.numeric(lambdaoptimres$par)
                        lambda_var =  1/(-lambdaoptimres$hessian)
                        en = as.numeric( (lambdaoptimres$par)^2/ lambda_var )
                        fn = as.numeric( lambdaoptimres$par/ lambda_var     )   
                    
                        
                        bn = b0 + 0.5*sum( err + 
                                            sum( diag( Icov %*%xkxk)  )   + 
                                            sum( as.numeric(lambda_var)*diag(xkxk) )
                                           )  
                        sigma = an/bn
                        
                        ##update s
                        dn = d0+ 0.5*(sum(w^2)  +  sum(diag(Icov) ) )
                        sn= cn/dn
                        ####
                        
                        #calc lower bound
                        lpy = (N/2)*( digamma(an) - log(bn) )- (sigma/2)*(err +  sum(diag( Icov%*%xkxk) )   + sum( as.numeric(lambda_var) *diag(xkxk) ) ) 
                        plbeta = (2/2)*(digamma(cn) - log(dn)) - ((cn/dn)/2)*(sum(w^2)  +  sum(diag(Icov) ) )
                            
                        lpsigma = (a0-1)*( digamma(an) - log(bn) ) - b0*(sigma)
                        lps = (c0-1)*( digamma(cn) - log(dn) ) - d0*(cn/dn)
                        lplambda =(e0-1)*( digamma(en) - log(fn) ) - e0*(lambda)
                        #----
                        esigma = an -log(bn)+lgamma(an)+(1-an)*digamma(an)
                        es = cn -log(dn)+lgamma(cn)+(1-cn)*digamma(cn)
                        elambda = en -log(fn)+lgamma(en)+(1-en)*digamma(en)
                        ebeta = (1/2)*log(det(Icov))
                         
                        bound = lpy+plbeta+lpsigma+lps+lplambda-esigma-es-elambda-ebeta
                      
                        bounds<<- c(bounds,bound)
                      #  print(bound )
                        if(lastBound>bound){
                          print('Bound should increase')
                        }
                        if(abs(lastBound-bound)<epsilon ){
                          break
                        }
                        lastBound = bound
                    }
                    
                    beta1 <<- betas
                    beta0 <<- beta0
                    sigma <<- sigma
                    s <<- sn
                    lambda <<- lambda
                    
                  },
                  predict=function(x,y){
                    N = nrow(x)
                    M = ncol(x)
                    w = matrix( .self$beta1* exp(- .self$lambda * 0:(M-1) ),ncol=1)
                    
                    pred = x%*% w +  .self$beta0
                    pred
                    
                  }
                )
  )
  

reg = bfreg()

reg$train(x,y,initialGuess=F,epsilon = 10e-5,maxiterations = 100)
plot(reg$bounds)

pred = reg$predict(x,y)
sum((pred-y)^2)


######mcmcm


printf <- function(...){
  invisible(print(sprintf(...)))
}

library("mvtnorm")

##Gibbs sampler
sampleGibbs = function(x,y,iterations = 1000 ){
  
  N = nrow(x)
  M = ncol(x)
  
  a0 = 0.001
  b0 = 0.001
  an = a0 + 0.5*N

  e0 = 0.001
  f0 = 0.001
  en = e0 + 0.5*2
  
  ##set up 
  lambda = rgamma(1,1,1)
  sn = rgamma(1,1,1)
  sigma =  rgamma(1,1,1)
  w = rmvnorm(1,rep(1,2),diag(2))
  
  fnlambda = function(lambda){
    dx = matrix( w[1]* exp(- lambda * 0:(M-1) ),ncol=1)
    err = sum( ( (x%*%dx +w[2]) -y)^2)
    res = -(sigma/2)*err  + (e0 -1)*log(lambda) - f0 *lambda
    res
  }

  updateLambda = function(lambda){
    
    oldprob =  fnlambda(lambda)
    repeat{
      newvalue =lambda  + rnorm(1,0,0.1)
      if(newvalue>0){
        break
      }
    }
    
    newprop =fnlambda(newvalue)
    ratio<-(newprop)-(oldprob)
    if(log(runif(1))<ratio) {
 #     print('accept')
      return(newvalue)
    }else{
   #   print('reject')
      return(lambda) 
    }
    
  }
  
  keep = matrix(,ncol=5,nrow=iterations)
  for(i in 1:iterations){
    if(i%%100==0){
      printf("%d / %d",i,iterations)
    }
    
    k  = c( exp(- lambda * 0:(M-1) ) )
    xnew  = x%*%k 
    xnew = cbind(xnew,1)
    
    xky = t(xnew)%*%(y)
    
    xkxk = t(xnew)%*%(xnew)
    dcov = (sigma)*(  xkxk )  + diag(  sn,2)
    
    Icov = solve(dcov)

    w = rmvnorm(1,Icov %*%  (sigma*xky) ,Icov)
    
    err = sum( (xnew%*%t(w) -y)^2 )

    sigma = rgamma(1,an,b0+0.5*err )

    ##update s
    sn =rgamma(1,en,f0+0.5*sum(w^2) )

    lambda = updateLambda(lambda)
    
    ###store
    keep[i,] = c(as.numeric(w),lambda, sn, sigma)
  }
  colnames(keep) = c('beta1','beta0','lambda','s','sigma')
  return(keep)
  
}

samples = sampleGibbs(x,y,12000)
samples = samples[2000:12000,]

plot(samples[,'beta1'])
plot(samples[,'beta0'])

colMeans(samples)

printf("%f %f %f %f %f",reg$beta1,reg$beta0,reg$lambda,reg$s,reg$sigma)


###jags model
library(mcmcplots)
library(runjags)
library(rjags)


reg_jags<- function(x,y){
  modelstring = "
  model{
  for( j in 1:N ){
    for(k in 1:10){
      seq[j,k] = exp( -lambda * weights[k] )
    }
    mu[j] =  x[j, ] %*% ( beta1 * seq[j,] )  
    y[j] ~ dnorm( mu[j] +beta0 , sigma) 
  }
  

  lambda~dgamma(0.001,0.001)
  sigma~dgamma(0.001,0.001)
  beta0 ~ dnorm(0,s)
  beta1 ~ dnorm(0,s)
  s~dgamma(0.001,0.001)
  }"
  writeLines(modelstring,con="decreg.txt")
 
 set.seed(123)
  
  p = ncol(X) ## 
  
  jags_data <- list(y = as.numeric(y),
                    x = x,
                    N=nrow(x),weights = seq(0,ncol(x)-1)
                    )
  params <- c("beta1",'sigma','s','lambda','beta0')                   
  adapt <- 5000
  burn <- 5000
  iterations <- 10000
  inits <- list( )
  
  sample <- run.jags(model="decreg.txt", thin =4, monitor=params, data=jags_data, n.chains=1, inits=inits, adapt=adapt, burnin=burn, sample=iterations, summarise=T, method="parallel") 
  
  sample
}


sample = reg_jags(x,y)
sample$summaries
