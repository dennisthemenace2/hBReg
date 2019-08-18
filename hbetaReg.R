## Hierarchical Beta Regression


## Variational implementation of a Hierarchical Beta Regression model.
# It utilizes Laplace approximation to fit the non-conjugate Beta distribution.

hBetareg <- setRefClass("hBetareg",
                    fields = list(Dcov= 'matrix', cov='list',beta='list', delta = 'matrix',sigma='numeric',si='numeric', a0 = 'numeric',b0='numeric' ,W='matrix',c0 = 'numeric',d0='numeric',e0 = 'numeric',f0='numeric',mse='numeric',iterations='numeric',lowerbound ='numeric',bounds='numeric',an='numeric',bn='numeric' ),
                    methods = list(
                      
                      # X Values
                      # Y Values
                      # id user id
                      train = function(x,y,id,initalGuess=T,deltaStop =1e-3 ) {#1e-1
                        C = max(id)
                        N = nrow(x)##samples ## == length(id)
                        p = ncol(x)##features
                        
                        ###check
                        if(any(y<=0 |y>=1 ) ){
                          print('Data out of valid range.....')
                          if(any(y==0)){
                            y[y==0] = 10^-9
                          }
                          if(any(y==1)){
                            y[y==1] = 1-10^-9
                          }
                        }
                        
                        
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
                        if(initalGuess==T){
                            
                          ##take a guess
                          y_new = log(y/(1-y))
                          x_new = as.matrix(x[,-ncol(x)],ncol=(ncol(x)-1))
                          x_new = cbind(x_new,y_new)
                          x_new = as.data.frame(x_new)
                          colnames(x_new)[ncol(x_new)]= 'Y'
                         
                          model = lm('Y~.',x_new )
                          .self$delta[p] = model$coefficients[1]
                          .self$delta[1:(p-1)] = model$coefficients[2:p]
                          print(.self$delta)
                          
                          olsfittednew = x%*%.self$delta 
                          olserrorvar = sum( (y_new-olsfittednew/(N-(p-1)) )^2 ) 
                          olsfitted = exp(olsfittednew) / (1 + exp(olsfittednew)) 
                          
                          vp1 =mean( olserrorvar*(olsfitted *(1.0-olsfitted))  )
                          .self$sigma = 20#vp1
                          print(vp1)
                          
                          an = as.numeric( (vp1)^2/ 10) # just assume somthing for the variance since its a point estimate
                          bn = as.numeric( vp1/ 10 )
                          print('initial guess')
                          print(.self$delta)
                          print(vp1)
                        }
                        
                        fnbeta = function( beta){
                          eta = X%*%beta 
                          mu = exp(eta) / (1.0+exp(eta)) 
                          diff = beta -.self$delta
                          res = sum ( -lgamma(mu*.self$sigma)  - lgamma((1-mu)*.self$sigma)+ (mu*.self$sigma-1)*log(Y) +  ( (1-mu)*.self$sigma -1 )*log(1-Y) )  -( .self$si/2)*(t(diff) %*%diff) 
                          res
                        }
                        fngrbeta= function(beta){
                          eta = X%*%beta 
                          mu = exp(eta) / (1.0+exp(eta)) 
                          ynew = log( Y / (1.0-Y) );
                          
                          munew = digamma(mu*.self$sigma) - digamma((1.0-mu)*.self$sigma)
                          T_mat = diag( as.numeric(exp(eta) / (1.0+exp(eta))^2) );
                          gradbeta = .self$sigma* t(X)%*%T_mat%*%(ynew-munew);
                          diff = beta -.self$delta
                          gradbeta = gradbeta - ( .self$si*diff)
                        }
                        fnphi = function(phi){
                          res = 0
                         # print(length(length(xclist)))
                          for( i in 1:length(xclist) ){
                            X = xclist[[i]]
                            Y = yclist[[i]]
                            beta = .self$beta[[i]]
                            eta = X%*%beta 
                            mu = exp(eta) / (1.0+exp(eta)) 
                            res =res+ sum (lgamma(phi)-lgamma(phi*mu) -lgamma((1-mu)*phi)+(mu*phi-1)*log(Y) +  ( (1-mu)*phi -1 )*log(1-Y) ) 
                          }
                        #  print(res)
                          res = res + (an-1)*log(phi) - phi*bn
                          #  print(res)
                          res
                        }
                        fnphigr = function(phi){
                          gradphi =0
                          for( i in 1:length(xclist) ){
                            X = xclist[[i]]
                            Y = yclist[[i]]
                            beta = .self$beta[[i]]
                            #print(nrow(X))
                            eta = X%*%beta 
                            mu = exp(eta) / (1.0+exp(eta)) 
                            gradphi = gradphi  + sum( digamma(phi) - mu *digamma(mu*phi) - (1.0-mu) * digamma( (1.0-mu)*phi) + mu * log(Y) + (1.0-mu) * log(1.0-Y) )
                          }
                          gradphi = gradphi+ (an-1)/phi-bn
                        }  
                        
                        fnphigr2 = function(phi){
                          eta = X%*%beta 
                          mu = exp(eta) / (1.0+exp(eta)) 
                          
                          gradphi = sum( trigamma(phi) - mu^2 *trigamma(mu*phi) - (1.0-mu)^2 * trigamma( (1.0-mu)*phi))
                          gradphi = gradphi - (an-1)/(phi^2) 
                        }  
                        fngr2beta= function(beta){
                          #mu = linkfn(beta,X)
                          eta = X%*%beta 
                          mu = exp(eta) / (1.0+exp(eta)) 
                          
                          psi1 = trigamma(mu*.self$sigma)
                          psi2 = trigamma((1.0-mu)*.self$sigma)
                          T_mat = diag( as.numeric(exp(eta) / (1.0+exp(eta))^2) );
                          
                          W = (diag(as.numeric(-.self$sigma^2*(psi1+psi2),ncol(T_mat)) ) %*% T_mat^2); 
                          tempinv = t(X)%*% W %*% X; 
                          tempinv -.self$si
                        }
                        ####prepare pre computed lists
                        xclist = list()
                        yclist = list()
                        xxlist = list()
                      #  xylist = list()
                        
                        .self$si = .self$c0 /.self$d0 
                        
                        ##maybe set up delta and betas here and the rest
                        
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
                      #    xylist[[c]] = t(xc)%*%yc
                        }
                        
                        
                        #### ESTIMATE GAMMA a parameters
                        cn = (.self$c0 + (C * p)/2)
                        #an = .self$a0 + (N)/2
                        en = (.self$e0 + 1/2)
                        an =.self$a0
                        bn =.self$b0
                        
                        ##compute variational bound   
                        lastbound = -Inf
                        .self$bounds= c(lastbound)    
                        for(i in 1:1000){
                          sse = 0
                          deltadiv = 0
                  #        sumVarBi = 0
                          sumVar = 0
                          ##
                          sumlogdetbi = 0
                   #       sumlpY =0
                          ##
                          for (c in 1:C) {
                            xc = xclist[[c]]
                            yc = yclist[[c]]
                           
                            beta = .self$beta[[c]]
                            X = xc
                            Y = yc
                            betaoptimres = optim(beta, fnbeta,method = "BFGS",gr=fngrbeta,
                                                 lower = -Inf, upper = Inf,
                                                 control = list(fnscale = -1), hessian = T)  
                            
                          #  betaoptimres$par
                            
                           # chkGradient(fnbeta,fngrbeta,beta)
                            ##check gradient.
                            
                            beta= betaoptimres$par
                            Icov = solve(-fngr2beta(betaoptimres$par))# 
                          #  Icov = solve(-betaoptimres$hessian) 
                            ##optimize phi
                            
                            .self$cov[[c]] = Icov
                          #  print(Icov)
                            .self$beta[[c]] = beta# matrix( Icov %*% (AX %*% .self$delta + .self$sigma * xy) )
                           # sse = sse + sum((xc %*% .self$beta[[c]] - yc)^2)
                            deltadiv = deltadiv + sum((.self$beta[[c]] - .self$delta)^2)
                            
                    #        sumVarBi = sumVarBi + diag(Icov %*%xxlist[[c]])
                            if(is.nan(sum(diag(Icov))) ){
                              stop()
                            }
                            sumVar = sumVar +sum(diag(Icov) )
                   
                            #if(is.nan(log(det(Icov))) ){
                            #  stop()
                            #}
                            detIcov = det(Icov)
                            if(detIcov>0){
                              sumlogdetbi = sumlogdetbi +log(detIcov)
                            }else{
                              print('det=0!')
                            }
                            ##
                          }
                          mvar = diag(.self$si,p) 
                          allb = matrix(unlist(.self$beta), ncol = C)
                          
                          meanb = matrix(rowSums(allb))
                          
                          mmu = mvar %*% meanb
                          mvar = mvar *C
                          
                          Dcov = .self$W + mvar
                          DIcov = solve(Dcov)
                          
                          .self$delta = DIcov %*% (mmu)
                          
                          dn = (.self$d0 + 0.5 * (deltadiv+ sum(diag(DIcov)) +  sumVar ) )
                          
                          .self$si = as.numeric(cn/dn )
                          
                          
                          phioptimres = optim(.self$sigma, fnphi,method = "BFGS",gr = fnphigr,
                                              #                         lower = 0, upper = Inf,
                                              control = list(fnscale = -1), hessian = T)  
                          #beta_var = -(betaoptimres$hessian)
                          
                          #chkGradient(fnphi,fnphigr,.self$sigma)
                          ##moments matching
                          .self$sigma= as.numeric(phioptimres$par)
                          phi_var =  -1/fnphigr2(phioptimres$par) #-(phioptimres$hessian)#
                       #   cat(c('phi_var: ',phi_var, 'phioptimres$hessian: ',phioptimres$hessian,'\n') )
                          #print(phioptimres)
                          an = as.numeric( (phioptimres$par)^2/ phi_var )
                          bn = as.numeric( phioptimres$par/ phi_var     )
                      #    cat(c('sigma:',.self$sigma,'\n'))
                          #bn = as.numeric(.self$b0 + 0.5 * (sse + sum(sumVarBi)  ))
                         # .self$sigma = an/bn
                          
                          fn  =as.numeric(.self$f0 + 0.5 * (.self$delta^2 +  diag(DIcov)))
                          tw = en/ fn
                          .self$W = diag(p) * tw
                       
                          lpY = 0
                          for(c in 1:C){ 
                            eta = xclist[[c]]%*%.self$beta[[c]] 
                            Y =  yclist[[c]]
                            mu = exp(eta) / (1.0+exp(eta)) 
                            lpY = lpY +  sum(lgamma(.self$sigma)-lgamma(.self$sigma*mu) -lgamma((1-mu)*.self$sigma)+(mu*.self$sigma-1)*log(Y) +  ( (1-mu)*.self$sigma -1 )*log(1-Y) ) 
                          }
                         # print(lpY)
                          ##calculate lower bound.
                          lpsigma = (.self$a0-1)*( digamma(an) - log(bn) ) - .self$b0*(.self$sigma)
                          lpsi = (.self$c0-1)*( digamma(cn) - log(dn) ) - .self$d0*(.self$si)
                          lpw =sum( (.self$e0-1)*( digamma(en) - log(fn) ) - .self$e0*(diag(.self$W) ) )
                          plDelta = sum( (p/2)*(digamma(en) - log(fn))-sum( (diag(.self$W)/2)*( .self$delta^2 +  diag(DIcov)) ) )
                          plbeta = (C/2)*(digamma(cn) - log(dn)) - (.self$si/2)*(deltadiv+ sum(diag(DIcov)) +  sum(sumVar) )
                          #----
                          esigma = an -log(bn)+lgamma(an)+(1-an)*digamma(an)
                          esi = cn -log(dn)+lgamma(cn)+(1-cn)*digamma(cn)
                          eW = sum(en - log(fn)+lgamma(en)+(1-en)*digamma(en) )
                          eDelta = (1/2)*log(det(DIcov))
                          ebeta = (1/2) *sumlogdetbi
                          
                       #   cat(c('an:',an,' bn:',bn,'\n'))
                      #    print('coefs')
                      #    cat(c(
                      #      lpsigma,' ',lpsi,' ',lpw,' ',lpY,' ',plDelta,' ',plbeta,' ',esigma,' ' , esi,' ' , eW,' ' ,eDelta,' ' , ebeta,'\n'
                      #    ))
                          
                          bound = lpsigma+lpsi+lpw+lpY+plDelta+plbeta -esigma - esi - eW -eDelta - ebeta
                          
                          .self$bounds= c(.self$bounds,bound)
                          if(lastbound > bound){
                            print('warning: lower bound should increase') 
                          }
                            print(bound)
                          if( abs(lastbound - bound) < deltaStop ){
                               cat(c('converged:',i ,'\n' ))
                            break;
                          }
                          
                          lastbound = bound
                        }
                        ##some statistics and bound
                        #.self$mse = sse
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
                        y = exp(y)/(1+exp(y))
                        y
                      }#,
                     #predictDelta = function(x) {   
                    #  y= x %*% .self$delta 
                    #  y = exp(y)/(1+exp(y))
                    #  y
                    # }
                      
                    ))




set.seed(1234)
####test case create some data


w = matrix(c(0.1,0.3) ,ncol=1)
X = matrix(rep(seq(0.1,1, length.out = 10),nrow(w)),ncol=nrow(w))
X[,nrow(w)]=1
phi = 1.2
eta = X%*%w #+rnorm( nrow(X),0,0.01)
mu = exp(eta)/(1+exp(eta))
Y = mu +rnorm(10,0,0.001)
#p =mu*phi
#q= (1-mu)*phi
#Y=  p/ (p+q ) 


w2 = matrix(c(0.12,0.28) ,ncol=1)
eta2 = X%*%w2 #+rnorm( nrow(X),0,0.01)
mu2 = exp(eta2)/(1+exp(eta2))
Y2 = mu2 +rnorm(10,0,0.001)

X = rbind(X,X)
Y = rbind(Y,Y2)
id = c(rep(1,10),rep(2,10))
#id = c(rep(1,10) )
#id = c(rep(1,32) )

reg = hBetareg()
reg$train(X,Y,id,initalGuess = F)
plot(reg$bounds)

#reg2 = hreg()
#reg2$train(X,Y,id)
#plot(reg2$bounds)
