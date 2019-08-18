##Implementation of the hierachical Bayesian Regression
## using Variational Approximation

hreg <- setRefClass("hreg",
                    fields = list(Dcov= 'matrix', cov='list',beta='list', delta = 'matrix',sigma='numeric',si='numeric', a0 = 'numeric',b0='numeric' ,W='matrix',c0 = 'numeric',d0='numeric',e0 = 'numeric',f0='numeric',mse='numeric',iterations='numeric',lowerbound ='numeric',bounds='numeric',an='numeric',bn='numeric' ),
                    methods = list(
                      
                      # X Values
                      # Y Values
                      # id user id
                      train = function(x,y,id) {
                        C = max(id)
                        N = nrow(x)##samples ## == length(id)
                        p = ncol(x)##features
                        
                        ###set priors ##set should be setable as parameters
                        .self$a0 = 0.001
                        .self$b0 = 0.001
                        
                        .self$c0=0.001
                        .self$d0=0.001
                        
                        .self$e0=0.001
                        .self$f0=0.001
                        ##########
                        
                        ###create inital setup
                        .self$sigma = .self$a0 / .self$b0
                        .self$W = (.self$e0 / .self$f0) * diag(p) 
                        
                        .self$delta = matrix(ncol=1 ,rep(0,p ))
                        
                        ####prepare pre computed lists
                        xclist = list()
                        yclist = list()
                        xxlist = list()
                        xylist = list()
                          
                         .self$si = .self$c0 /.self$d0 
                              
                        for(c in 1:C){ ##create each client matrix
                          .self$beta[[c]]= as.matrix(.self$delta) ## plus intercept 
                                                    
                          xc =  as.matrix ( x[id==c,])
                          if(ncol(xc)==1){
                            xc = t(xc)
                          }
                          yc = y[id==c]
                          xclist[[c]] = xc
                          yclist[[c]] = yc
                          xxlist[[c]] = t(xc)%*%xc
                          xylist[[c]] = t(xc)%*%yc
                        }
                       
                        
                        #### ESTIMATE GAMMA a parameters
                        cn = (.self$c0 + (C * p)/2)
                        an = .self$a0 + (N)/2
                        en = (.self$e0 + 1/2)
                     
                        ##compute variational bound   
                        lastbound = -Inf
                        .self$bounds= c(lastbound)    
                        for(i in 1:1000){
                          sse = 0
                          deltadiv = 0
                          sumVarBi = 0
                          sumVar = 0
                          ##
                          sumlogdetbi = 0
                          ##
                          for (c in 1:C) {
                            xc = xclist[[c]]
                            yc = yclist[[c]]
                            xx = xxlist[[c]]
                            xy = xylist[[c]]
                            AX = .self$si * diag(p)
  
                            cov = AX + .self$sigma * xx
                            Icov = solve(cov)
                            .self$cov[[c]] = Icov
                            .self$beta[[c]] = matrix( Icov %*% (AX %*% .self$delta + .self$sigma * xy) )
                            sse = sse + sum((xc %*% .self$beta[[c]] - yc)^2)
                            deltadiv = deltadiv + sum((.self$beta[[c]] - .self$delta)^2)
                          
                            sumVarBi = sumVarBi + diag(Icov %*%xx)
                            sumVar = sumVar +sum(diag(Icov))
                            ##
                            sumlogdetbi = sumlogdetbi +log(det(Icov))
                          #  stop("weit")
                            ##
                          }
                          mvar = .self$si * diag(p) 
                          allb = matrix(unlist(.self$beta), ncol = C)
                         
                          meanb = matrix(rowSums(allb))
                          
                          mmu = mvar %*% meanb
                          mvar = mvar *C
                      
                          Dcov = .self$W + mvar
                          DIcov = solve(Dcov)
                        
                         .self$delta = DIcov %*% (mmu)
                        
                          dn = (.self$d0 + 0.5 * (deltadiv+ sum(diag(DIcov)) +  sumVar ) )
                      
                         .self$si = as.numeric(cn/dn )
                   
                          bn = as.numeric(.self$b0 + 0.5 * (sse + sum(sumVarBi)  ))
                        
                          .self$sigma = an/bn

                          fn  =as.numeric(.self$f0 + 0.5 * (.self$delta^2 +  diag(DIcov)))
                          tw = en/ fn
                          .self$W = diag(p) * tw
                                              
                          
                         ##calculate lower bound.
                         lpsigma = (.self$a0-1)*( digamma(an) - log(bn) ) - .self$b0*(.self$sigma)
                         lpsi = (.self$c0-1)*( digamma(cn) - log(dn) ) - .self$d0*(.self$si)
                         lpw =sum( (.self$e0-1)*( digamma(en) - log(fn) ) - .self$e0*(diag(.self$W) ) )
                         lpY= (N/2)*( digamma(an) - log(bn) )- (.self$sigma/2)*(sse + sum(sumVarBi)  )
                         plDelta = sum( (p/2)*(digamma(en) - log(fn))-sum( (diag(.self$W)/2)*( .self$delta^2 +  diag(DIcov)) ) )
                         plbeta = (C/2)*(digamma(cn) - log(dn)) - (.self$si/2)*(deltadiv+ sum(diag(DIcov)) +  sum(sumVar) )
                         #----
                         esigma = an -log(bn)+lgamma(an)+(1-an)*digamma(an)
                         esi = cn -log(dn)+lgamma(cn)+(1-cn)*digamma(cn)
                         eW = sum(en - log(fn)+lgamma(en)+(1-en)*digamma(en) )
                         eDelta = (1/2)*log(det(DIcov))
                         ebeta = (1/2) *sumlogdetbi
                      
                       
                         bound = lpsigma+lpsi+lpw+lpY+plDelta+plbeta -esigma - esi - eW -eDelta - ebeta
                       
                         .self$bounds= c(.self$bounds,bound)
                         if(lastbound > bound){
                           print('warning: lower bound should increase') 
                         }
                       #  print(bound)
                         if( abs(lastbound - bound) < 1e-10 ){
                        #   cat(c('converged:',i ,'\n' ))
                           break;
                         }

                         lastbound = bound
                        }
                        ##some statistics and bound
                        .self$mse = sse
                        .self$iterations = i
                        .self$lowerbound = bound
                        .self$an = an
                        .self$bn = bn
                        .self$Dcov =DIcov
                      },
                      ##predict unseen data
                      #x value
                      #id id of the user
                      # this function does not implement the prediction for users that have not been observed.
                      # in this case the hierachical prior delta should be used for prediction
                      predict = function(x,id) {   
                        y= x %*% .self$beta[[id]] 
                        y
                      },
                      predictBoth = function(x,id) {   
                        y= x %*% .self$beta[[id]] 
                        ret = list()
                        ret$y = y
                        
                        tmp =  1/( (.self$an)/.self$bn  + x%*% .self$cov[[id]]%*%t(x) )
                       #print(sqrt(var))
                        left = y + qt(p=0.025,lower.tail = T,df=2*.self$an) * sqrt(tmp)
                        right = y+ qt(p=0.025,lower.tail = F,df=2*.self$an) * sqrt(tmp)
                        ret$left = left
                        ret$right=right
                        ret$var = tmp
                        varD= 1/( (.self$an)/.self$bn  + x%*% .self$Dcov%*%t(x) )
                        varS = 1/( (.self$an)/.self$bn  + x%*% diag(1/.self$si,length(x))%*%t(x) )
                        ret$varD = varD
                        ret$vars = varS
                        print(tmp)
                        print(varD)
                        print(varS)
                        
                        ret
                      }
                      
                    ))

set.seed(1234)
####test case create some data
hp=c(0.1,0.3)
p1= hp+ rnorm(2,0,0.01)
p2= hp+ rnorm(2,0,0.01)
p3= hp+ rnorm(2,0,0.01)
p4= hp+ rnorm(2,0,0.01)
p5= hp+ rnorm(2,0,0.01)
p6= hp+ rnorm(2,0,0.01)
  
N = 100
Y1 = p1[1]*(1:N) +p1[2] + rnorm(N,0,0.1)
Y2 = p2[1]*(1:N) +p2[2] + rnorm(N,0,0.1)
Y3 = p3[1]*(1:N) +p3[2] + rnorm(N,0,0.1)
Y4 = p4[1]*(1:N) +p4[2] + rnorm(N,0,0.1)
Y5 = p5[1]*(1:N) +p5[2] + rnorm(N,0,0.1)
Y6 = p6[1]*(1:N) +p6[2] + rnorm(N,0,0.1)


X = (1:N)
Xp = rbind(X, rep(1,N))

X = Xp
X = cbind(X, Xp)
X = cbind(X, Xp)
X = cbind(X, Xp)
X = cbind(X, Xp)
X = cbind(X, Xp)
X = t(X)
Y = as.matrix(c(Y1,Y2,Y3,Y4,Y5,Y6))
id= c(rep(1,N),rep(2,N),rep(3,N),rep(4,N),rep(5,N),rep(6,N))

reg = hreg()

reg$train(X,Y,id)
##check delta

cat(c('Variational Model:','\n'))
cat(c('Delta:',hp,'\n'))
cat(c('Delta estimate:',reg$delta,'\n'))

plot(reg$bounds)

plot(Y[id==1])
pred = c()
left = c()
right = c()
for(i in which(id==1)){
  ret = reg$predictBoth(matrix(X[i,],nrow=1) ,1)
  pred = c(pred, ret$y)
  left = c(left, ret$left)
  right = c(right, ret$right)
}


lines(pred,col='red')
lines(left,col='green')
lines(right,col='green')
 
p1 = cbind((1:N), rep(1,N))
data = as.data.frame(cbind(Y1,p1[,1]) )
colnames(data) = c('Y','X1')
model = lm("Y~.", data)
test = as.data.frame(p1[,1])
colnames(test) = c('X1')
lmpred = predict(model, test, interval = "confidence")

lines(lmpred[,1],col='blue')
lines(lmpred[,2],col='blue')
lines(lmpred[,3],col='blue')


##Gibbs sampler implementation of the same model
hregGIBBS<-function(Y,X,a0=0.001,b0=0.001,c0=0.001,d0=0.001,e0=0.001,f0=0.001,n.samples=1000,id){
  
  library("mvtnorm")

  N <- length(Y)
  p <- ncol(X)
  C = max(id)
  
  ##precomputed stuff
  xx = list()
  xy = list()
  keep.bi = list()
  
  ##check if all ids are there
  if(length(unique(id)) != max(id)){
    print('the ids does not seam in order')
    stop()
  }
  ###
  for(i in 1:C){
    x = X[id==i ,]
    y = Y[id==i ]
    x = as.matrix(x)
    y = as.matrix(y)
    
    xx[[i]] = t(x)%*%x
    xy[[i]] = t(x)%*%y
    keep.bi[[i]] <-  matrix(0,n.samples,p)
  }
  
  #intial values:
  sigma <- rgamma(1,1,1)
  si <- rgamma(1,1,1)
  
  
  delta   <-  as.matrix(rnorm(p,0,1),ncol = 1 ) # solve( t(X)%*%X) %*% t(X)%*%Y ##could be used for an initial guess
  w = rgamma(p,1,1)
  
  # beta   <- rnorm(p,0,10)
  
  #Initialize vectors to store the results 
  #keep.bi <-  matrix(0,n.samples,C)
  
  keep.sigma = rep(0,n.samples)
  keep.si = rep(0,n.samples)
  keep.delta   <- matrix(0,n.samples,p)
  keep.w  <- matrix(0,n.samples,p)
  keep.log =rep(0,n.samples)
  
  keep.prediction = list()
  for(i in 1:sum(id==1) ){
    keep.prediction[[i]] =   rep(0,n.samples)
  }

  test = matrix(c(1,1))
  
  # keep.beta   <- matrix(0,n.samples,p)
  
  #Start the MCMC sampler!
  betai =list()
  betaiCov =list()
  
  for(i in 1:n.samples){
    
    sse = 0; 
    ddiv = 0;
    
    for(c in 1:C){
      cxx =  xx[[c]]
      cxy =  xy[[c]]
      y = Y[id==c ]
      x = X[id==c ,]
      x = as.matrix(x)
      y = as.matrix(y)
      cov = diag(p) * si +  sigma *cxx
      Icov = solve(cov)
      mu = (diag(p) * sigma) %*% cxy + (diag(p) * si) %*% delta
      betai[[c]] = rmvnorm(1,Icov %*% mu,Icov)   
      betaiCov[[c]] = Icov
      
      mse = sum( (y- x %*% t(betai[[c]]) )^2 )
      sse= sse+ mse

      ddiv = ddiv +sum( (betai[[c]] - t(delta) )^2 )
      ##keep bi
      tmp =keep.bi[[c]]
      tmp[i,] =  betai[[c]]
      keep.bi[[c]]=tmp
    }
    # update sigma
    sigma <- rgamma(1,a0+N/2,b0+0.5*sse);
    ##updeate si
    si <- rgamma(1,c0+(C*p)/2,d0+0.5*ddiv);
    #update delta
    dcov = diag(w) + (si*diag(p)*C)
    allb = matrix(unlist(betai), ncol = C)

    meanb = matrix(rowSums(allb))
    
    dmu  =    (diag(p) *si) %*%  meanb 
    Idcov = solve(dcov)
     
    delta = matrix(rmvnorm(1,Idcov %*% dmu,Idcov),ncol = 1 )
    
    #update w
    for(wp in 1:p){
      w[wp] = rgamma(1,e0+0.5,f0+0.5*delta[wp]^2 )
    }
    
    ##store results:
    keep.sigma[i] = sigma  
    keep.si[i] = si
    keep.delta[i,]=delta
    keep.w[i,] =w
    
  ##used for prediction of a target
    x = X[id==1 ,]
    for(t in 1:nrow(x)){
      #print(dim(x)) ##betai[[1]] 0.0003596903 / 0.01033604 # delga 0.000246839/ 0.0102771
      test= x[t,]
      keep.prediction[[t]][i] = rnorm(1, keep.bi[[1]][i,]%*% test ,1/(sigma) )
    }
  }
  #return a list with the posterior samples:
  retlist = list(sigma=keep.sigma,si=keep.si,delta = keep.delta, w= keep.w ,loglik=keep.log,prediction=keep.prediction  ) 
  for(i in 1:C){
    retlist[[length(retlist)+1]]= keep.bi[[i]]
    names(retlist)[length(retlist)] =  paste('beta',i,sep='')
  }
  retlist
}

### small testcase
  
n.samples =10000
fit<- hregGIBBS(Y,X,n.samples=n.samples,id=id)


cat(c('Gibbs Sampler:','\n'))
cat(c('Delta:',hp,'\n'))
cat(c('Delta estimate:',colMeans(fit$delta),'\n'))

samplerPred = c()
samplerVar = c()
for(i in 1:length(fit$prediction)){
  samplerPred = c(samplerPred,mean(fit$prediction[[i]]) )
  samplerVar = c(samplerVar,var(fit$prediction[[i]]) )
}
lines(samplerPred, col='orange')
lines(samplerPred+ sqrt(samplerVar)*qnorm(0.025), col='red')
lines(samplerPred- sqrt(samplerVar)*qnorm(0.025), col='red')

plot(fit$delta)
require(mixtools)
ellipse(mu=colMeans(fit$delta), sigma=cov(fit$delta), alpha =  .05, npoints = 250, col="red")
ellipse(mu=as.numeric(reg$delta), sigma=reg$Dcov, alpha = .05, npoints = 250, col="blue")

plot(fit$beta1)
ellipse(mu=colMeans(fit$beta1), sigma=cov(fit$beta1), alpha = .05, npoints = 250, col="red")
ellipse(mu=as.numeric(reg$beta[[1]]), sigma=reg$cov[[1]], alpha = .05, npoints = 250, col="blue")



