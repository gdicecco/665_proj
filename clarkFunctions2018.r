
logit <- function(x) log(x/(1 - x))
invlogit <- function(x) exp(x)/(1 + exp(x))

getJagsPars <- function(fit){
  
  x <- fit$BUGSoutput$summary
  w <- grep('[', rownames(x), fixed=T)
  
  if(length(w) < nrow(x)){
    fixed <- signif(x[-w,], 4)
  }
  
  out <- list(fixed = fixed)
  
  if(length(w) > 0){
    
    z <- x[w,]
    cf  <- matrix(unlist(strsplit(rownames(z),'[', fixed=T)), ncol=2, byrow=T)
    idx <- unlist(strsplit(cf[,2],']', fixed=T))
    idx <- as.numeric(idx)
    cf  <- cf[,1]
    
    ip <- sort(unique(cf))
    id <- sort(unique(idx))
    
    n <- length(id)
    imat <- matrix(NA, n, length(ip))
    imat[ cbind(match(idx,id),match(cf,ip)) ] <- z[,1]
    colnames(imat) <- ip
    
    lo <- hi <- imat*0
    lo[ cbind(match(idx,id),match(cf,ip)) ] <- z[,'2.5%']
    hi[ cbind(match(idx,id),match(cf,ip)) ] <- z[,'97.5%']
    
    out$mean <- signif(imat, 4)
    out$ciLo <- signif(lo, 4)
    out$ciHi <- signif(hi, 4)
  }
  out
}

metRatio <- function(beta, x, y, likelihood, link,
                     priorB=beta*0, priorVB=diag(1000,length(beta))){
  B <- T
  if(likelihood == 'dpois')B <- F
  
  link  <- match.fun(link)
  like  <- match.fun(likelihood)
  theta <- link(x %*% beta)      
  
  if(B){
    p1 <- like(y, 1, theta, log=T)
  }else{
    p1 <- like(y, theta, log=T)
  }
  sum( p1 ) + .dMVN(beta, priorB, priorVB, log=T)
}
updateBetaGLM <- function(bg, cmat, x, y, likelihood='binom',
                          link = 'logit', lo = NULL, hi=NULL, ...){
  
  # metropolis sampler for GLM reg parameters
  # lo, hi - truncation points for MVN proposals
  
  ac <- 0
  if(is.null(lo)){             # proposal
    bs <- t(.rMVN(1,bg,cmat))        
  }else{
    bs <- .tnormMVNmatrix(t(bg), t(bg), smat=cmat, lo = t(lo), hi=t(hi))
    bs <- t(bs)
  }
  
  if(link == 'log')link <- 'exp'        # invert link function
  if(link == 'logit')link <- 'invlogit'
  pnow <- metRatio(bg, x, y, likelihood, link, ...)
  pnew <- metRatio(bs, x, y, likelihood, link, ...)
  a    <- exp(pnew - pnow)
  z    <- runif(1,0,1)
  if(a > z){
    bg <- bs
    ac <- 1
  }
  list(beta = bg, accept = ac)
}
.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}


.sqrtRootMatrix <- function(xmat,sigma,DIVIDE=F){
  
  # xmat from or to correlation scale
  # xmat is n by p
  # sigma is p by p
  
  if(DIVIDE){
    if(length(sigma) == 1)return(xmat/sqrt(sigma))
    return( xmat%*%diag(1/sqrt(diag(sigma))) )
  }
  
  if(length(sigma) == 1)return(xmat*sqrt(sigma))
  
  xmat%*%diag( sqrt(diag(sigma)) )
}

      
################### BBS #################

bbsOxV <- function(data){
  
  # data has columns:
  #        species codes 'Aou'
  #        routes 'Route'
  #        years 'Year'
  # returns plot-years by species counts
  
  allCodes <- sort(unique(data$Aou))              # species codes
  S        <- length(allCodes)                    # no. species
  rt_yr    <- paste(data$Route,data$Year,sep='_') # each route/yr is an obs
  allObs   <- sort(unique(rt_yr))                 # each obs
  n        <- length(allObs)                      # no. obs
  x        <- matrix(0,n,S)                       # OxV
  i        <- match(rt_yr,allObs)                 # row index
  j        <- match(data[,'Aou'],allCodes)        # col index
  ij       <- cbind(i,j)
  x[ij]    <- data[,'count']                      # fill obs matrix
  colnames(x) <- allCodes                         # name everything
  rownames(x) <- allObs
  rtYr     <- matrix( unlist(strsplit(allObs,'_')),n,2,byrow=T) # rownames to rt, yr
  route  <- as.numeric(rtYr[,1])
  year   <- as.numeric(rtYr[,2])
  cbind(route, year, x)
}

bbsSpecs <- function(path="dataBBS/"){
  
  f      <- paste(path,'SpeciesList.txt',sep='')
  tt     <- replaceNonAscii(f)$x
  widths <- c(0,6,6,51,51,51,51,51,51,51)
  specs  <- vector('list', length = length(widths)-1)
  start  <- cumsum(widths)
  for(i in 1:(length(widths)-1)){
    ti <- substr(tt,start[i]+1,start[i+1])
    names(specs)[i] <- gsub(" ", "", ti[1], fixed = TRUE)
    t1 <- gsub(" ", "", ti[-1], fixed = TRUE)
    if(i < 3)t1 <- as.numeric(t1)
    specs[[i]] <- t1
  }
  s <- as.data.frame(specs)[,c(1:3,6:9)]
  colnames(s)[1:3] <- c('seq','Aou','name')
  s
}

##################### Bayesian regression, with Tobit ##############

bayesReg <- function(formula, data, ng = 3000, burnin = 100, TOBIT=NULL){
  
  fc <- as.character(formula)
  
  if(missing(data)){
    data <- environment(formula)
    yx   <- match.call()
    m    <- match(c("formula", "data"), names(yx), 0L)
    yx   <- yx[c(1L, m)]
    yx[[1L]] <- quote(stats::model.frame)
    yx   <- eval(yx, parent.frame())
    y    <- yx[,1]
  } else {
    yy <- unlist( strsplit( fc, '~' ) )
    yy <- yy[ nchar(yy) > 0]
    y  <- data[,yy[1]]
  }
  
  yzero <- which (y == 0)
  nzero <- length(yzero)
  
  if(is.null(TOBIT)){
    TOBIT <- F
#    if(nzero > 0)TOBIT <- T
  } 
  
  if(TOBIT)message('fitted as Tobit model')
  
  tmp <- model.frame(formula, data, na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames    <- colnames(x)
  snames    <- colnames(y)
  Q         <- ncol(x)
  n         <- nrow(x)
  predXcols <- 2:Q
  
  ymiss <- which(is.na(y))
  
  yy <- y
  xx <- x
  if(length(ymiss) > 0){
    yy <- yy[-ymiss]
    xx <- xx[-ymiss,]
  }
  
  XX  <- crossprod(xx)
  IXX <- solve(XX)
  bg  <- IXX%*%crossprod(xx,yy)
  py  <- x%*%bg
  w   <- y
  w[ymiss] <- mean(y,na.rm=T)
  if(TOBIT){
    p0 <- py[y == 0]
    p0[p0 > 0] <- -p0[p0 > 0]
    py[y == 0] <- p0
    w[w == 0] <- py[y == 0]
  }
  wi  <- c(which(y == 0),ymiss)   
  
  priorB   <- matrix(0,Q,1)         #prior mean regression parameters
  priorIVB <- solve(diag(1000,Q))   #inverse prior covariance matrix
  s1       <- .1                    #variance prior values
  s2       <- .1
  
  bchains <- matrix(NA,ng,Q)
  schains <- rep(0,ng)  #store parameters
  colnames(bchains) <- xnames
  ychains <- matrix(0,ng,n)
  
  
  for(g in 1:ng){
    
    sg <- .updateSigma(x, w, bg, s1, s2)
    bg <- .updateBeta(x, w, sg, priorIVB, priorB, XX)
    mu <- x%*%bg
    if(TOBIT)w[yzero] <- .tnorm(nzero,-Inf,0,mu[yzero],sqrt(sg))
    w[ymiss] <- py[ymiss]
    py <- rnorm(n,mu,sqrt(sg))
    
    py[py <= 0] <- 0
    
    bchains[g,] <- bg   #store estimates
    schains[g]  <- sg
    ychains[g,] <- py
  }
  
  beta <- signif( t( apply(bchains,2,quantile,c(.5,.025,.975)) ), 4)
  bhat <- beta[drop=F,,1]
  bsd  <- signif( apply(bchains,2,sd ), 4)
  beta <- cbind(beta[,1],bsd,beta[,2:3])
  
  shat <- mean(schains)
  sse  <- sd(schains)
  
  zero <- which(beta[,3] < 0 & beta[,4] > 0)
  notZero <- rep('*',Q)
  notZero[ zero ] <- ' '
  beta <- data.frame( cbind(beta,notZero) )
  
  colnames(beta) <- c('median','std error','0.025','0.975','not zero')
  
  py <- signif( t( apply( ychains, 2, quantile,c(.5,.025,.975) ) ), 3)
  py <- cbind(y,py)
  
  rmspe <- signif( sqrt( mean((y - py[,2])^2, na.rm=T) ), 4)
  sig   <- signif(sqrt(shat),4)
  
  ssr   <- crossprod(y - x%*%bhat)/(n - Q)
  rsq   <- signif(1 - ssr/var(y), 3)
  
  cat("\nCoefficients:\n")
  print( beta )
  
  cat("\n * indicates that 95% predictive interval does not include zero\n")
  
  out <- paste("\nResidual standard error ", sig, ", with ",n - Q, 
               " degrees of freedom, \n root mean sq prediction error ",
               rmspe, ".", sep='')
  cat(out,"\n")

  
  list(beta = beta, predictY = py, sigma = median(schains), 
       rmspe = rmspe)
}

pmake <- function(pars){  
  
  p   <- pars[1]                   #frequency of a
  f   <- pars[2]                   #inbreeding coefficient 
  paa <- p^2 + f*p*(1 - p)	  #frequency of p.aa
  pab <- 2*p*(1 - p)*(1 - f)	  #frequency of p.ab
  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb
  c(paa,pab,pbb)
}

minf <- function(p){   #minimum inbreeding coefficient
  
  lof <- rep(-1,length(p))
  lop <- -(1 - p)/p
  hip <- -p/(1 - p)
  lof[p > (1 - p)] <- lop[p > (1 - p)]
  lof[p < (1 - p)] <- hip[p < (1 - p)]
  lof
}

pf.update <- function(p,f){  #log posterior
  
  multlik(c(p,f)) + dbeta(p,g1,g2,log=T) + 
    dnorm(f,priorF,priorFSD,log=T)
}

multlik <- function(pars){  #multinom likelihood
  
  p    <- pmake(pars)
  pmat <- matrix(rep(p,npop),3,npop)
  sum(y*log(pmat))
  
}

update_pf <- function(){
  
  ac     <- 0
  propp  <- .tnorm(1,.02,.98,pg,.002) #propose pa
  fl     <- minf(propp)              #lower f limit for pa
  propf  <- .tnorm(1,fl,1,fg,.01)     #propose f > fl
  pnew   <- pf.update(propp,propf)
  pnow   <- pf.update(pg,fg)
  a      <- exp(pnew - pnow)
  z      <- runif(1,0,1)
  if(z < a){
    pg <- propp
    fg <- propf
    ac <- 1
  }
  
  list(pg = pg, fg = fg, accept = ac)
}



###################### general functions ###################

rbetaBinom <- function(n, m, a=NULL, b=NULL, mu=NULL){
  
  # need (a, b), (mu, a), or (mu, b)
  
  if(!is.null(mu)){
    if(is.null(a) & is.null(b))stop('if mu is supplied, need also a or b')
    if(is.null(a))a <- b/(mu - 1)
    if(is.null(b))b <- a*(mu - 1)
  }
  rbinom(n, m, rbeta(n, a, b))
}
dbetaBinom <- function(x, m, a=NULL, b=NULL, mu=NULL, log=F){
  
  # need (a, b), (mu, a), or (mu, b)
  
  if(!is.null(mu)){
    if(is.null(a) & is.null(b))stop('if mu is supplied, need also a or b')
    if(is.null(a))a <- b/(1/mu - 1)
    if(is.null(b))b <- a*(1/mu - 1)
  }
  xx <- lchoose(m,x) + lbeta(x + a, m - x + b) - lbeta(a, b)
  if(!log)xx <- exp(xx)
  xx
}
.colorSegment <- function(xx, yy, colors, nc = 20, lwd=1, lty=1,
                          clear=T){
  
  # add line to plot, colored by magnitude
  # colors - colors to interpolate, see colorRampPalette
  #          interpolate for range(y)
  # x, y   - vectors to plot
  # clear  - clear background of line
  
  yr    <- range(yy)
  yr[1] <- yr[1] 
  yr[2] <- yr[2] 
  sc    <- seq(yr[1],yr[2],length=nc)
  yc    <- c(1:nc)
  y1    <- 1:(nc-1)
  y2    <- y1 + 1
  sc    <- round( yr[1] + yc/nc*diff(yr) ,0)
  
  cf   <- colorRampPalette(colors)(nc)
  ycol <- findInterval(yy,sc,all.inside=T)
  
  if(clear) lines( new$density, new$mids, lwd=lwd*3, col='white')
  segments(xx[y1], yy[y1], xx[y2], yy[y2], col=cf[ycol],lwd=lwd,lty=lty)
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

.combineFacLevels <- function(xfactor,fname=NULL, 
                              aname = 'reference', vminF=1){
  tmp <- as.character(xfactor)
  tmp[tmp %in% fname] <- aname
  tab <- table(tmp)
  wm  <- names(tab)[tab < vminF]
  tmp[tmp %in% wm] <- aname
  as.factor(tmp)
}

covarGrid <- function(lonlat, path='../dataFiles/'){
  
  require(RANN)
  file <- paste(path,'UScovarGrid.rdata',sep='')
  load(file)
  
  ii <- nn2(gridded.data[,c('lon','lat')], lonlat,k=1)[[1]]
  xx <- gridded.data[ii,c('soil','temp','prec','therm','def','nlcd')]
  colnames(xx)[colnames(xx) == 'temp'] <- 'winterTemp'  # two temp variables
  xx
}

replaceNonAscii <- function(filename){
  l   <- 'latin1'
  a   <- 'ASCII'
  x   <- readLines(filename)
  nar <- grep("notASCII", iconv(x, l, a, sub="notASCII"))
  if(length(nar) > 0) x[nar] <- iconv(x[nar], l, a, sub="-")
  list(notAsciRows = nar, x = x)
}


.processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                         sigOnly = F,burnin=1,xlimits = NULL){  
  
  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  NOPARS <- F
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
      btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }
    
    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    
    if(length(wq) == ncol(btmp))NOPARS <- T
    if(NOPARS) warning('no significant pars to plot')
    if(length(wq) > 0 & !NOPARS){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
  }
  
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
    if(burnin > (nrow(xgb) + 100))stop("burnin too large")
    xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)
  
  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }
  
  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(3,2,1,1))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(3,2,1,1))
  
  oseq <- 1:nrow(xgb)
  if(nrow(xgb) > 1000){
    oseq <- round(seq(1, nrow(xgb), length=1000))
    oseq <- unique(oseq)
  }
  
  if(CPLOT & !NOPARS){
    for(j in 1:nc){
      plot(xgb[oseq,j],type='l')
      abline(h=out[j,],lty=2)
      if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
      abline(h = 0, col='grey',lwd=2)
      title(colnames(xgb)[j])
    }
  }
  xlims <- xlimits
  if(DPLOT & !NOPARS){
    for(j in 1:nc){
      xj <- density(xgb[,j])
      if(is.null(xlimits))xlims <- range(xj$x)
      plot(xj$x,xj$y,type='l',xlim=xlims)
      abline(v=out[j,],lty=2)
      if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
      title(colnames(xgb)[j])
    }
  }
  list(summary = signif(out,4)
  )
  
}

.dMVN <- function(xx,mu,smat=NULL,sinv=NULL,log=F){ #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  if(!is.matrix(mu))mu <- matrix(mu,1)
  
  xx <- xx - mu
  if(ncol(xx) != nrow(smat))xx <- t(xx)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logd    <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
      tiny  <- min(abs(xx))/100 + 1e-5
      smat  <- smat + diag(diag(smat + tiny))
      testv <- try(chol(smat),T)
    }
    covm    <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev      <- eigen(smat, only.values = TRUE)$values 
    logd    <- sum(log( ev ))
  }
  
  z <- -(ncol(xx) * log(2 * pi) + logd + distval)/2
  if(!log)z <- exp(z)
  z
}


rnormCond <- function(a, b, sigma, xx, ww, column=1){
  
  mu  <- b%*%xx + a%*%ww       # correlation scale
  tmp <- conditionalMVN(mu,mu,sigma,column)
  rnorm(1, tmp$mu, sqrt(tmp$vr) )
}

conditionalMVN <- function(xx, mu, sigma, cindex){  
  
  # x and mu are vectors, cindex is vector index for conditional
  
  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(xx) != nm)stop('x and mu different length in conditionalMVN')
  
  xx <- matrix(xx,nrow=1)
  mu <- matrix(mu,nrow=1)
  
  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
    return( list(mu = numeric(0), vr = numeric(0)) )
  }
  
  sin <- chol2inv(testv)
  p1  <- sigma[cindex,-cindex]%*%sin
  
  mu1 <- mu[cindex] +  p1%*%(xx[-cindex] - mu[-cindex]) 
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[-cindex,cindex]
  
  list(mu = mu1, vr = vr1)
}

.rMVN <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}

.rwish <- function(df,SS){
  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

.riwish <- function(df,S){
  solve(.rwish(df,solve(S)))
}



.tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  tiny <- 10e-6
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}

.updateBeta <- function(x, y, sigma, priorIVB, priorB, 
                        XX=NULL){  # random vector of coefficients
  
  if(is.null(XX))XX <- crossprod(x)
  
  V  <- solve( XX/sigma + priorIVB ) 
  v  <- crossprod(x,y)/sigma + priorIVB%*%priorB
  t( .rMVN(1,t(V%*%v),V) )                     
}

.updateSigma <- function(x, y, beta, s1, s2){ # random value for residual variance
  n  <- length(y)
  u1 <- s1 + n/2
  u2 <- s2 + 0.5*crossprod(y - x%*%beta)
  1/rgamma(1,u1,u2)                          
}


smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind =T)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}
shadeInterval <- function(xvalues,loHi,col='grey',ylim=range(loHi),
                          PLOT=T,add=T,xlab=' ',ylab=' '){
  
  #draw shaded interval
  
  tmp <- smooth.na(xvalues,loHi)
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(!add)plot(xvalues,loHi[,1]*0,cex=.01,ylim=c(range(loHi,na.rm=T)),
               xlab=xlab,ylab=ylab)
  if(PLOT)polygon(xbound,ybound,border=NA,col=col)
  
  invisible(cbind(xbound,ybound))
  
}




tnorm.mvt <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),
                      hi=rep(Inf,length(avec)),
                      whichSample=c(1:length(avec)),times=1){   
  
  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample
  
  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))
  
  for(j in 1:times){
    for(k in whichSample){
      
      tmp <- conditionalMVN(avec,muvec,smat,k)
      muk <- tmp$mu
      sgk <- tmp$vr
      
      if(length(muk) == 0)next
      
      avec[k] <- .tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
    }
  }
  avec
}
.updateBetaMVN <- function(X,Y,sg){
  
  XX  <- crossprod(X)
  IXX <- chol2inv(chol( XX ) )
  WX  <- crossprod(X,Y)
  WIX <- IXX%*%WX
  matrix( .rMVN(1,as.vector(WIX),
                kronecker(sg,IXX)),nrow(IXX),ncol(WIX) )
}

varPrior <- function(mu, wt){#mean and wt of prior IG distribution
  c(wt, mu*(wt - 1))
}

initialStatesSS <- function(y){
  
  require(stats)
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)
  
  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y
  
  if(length(wm) > 0){
    notMiss <- notMiss[-wm]
    x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }
  
  list(x = as.vector(x), notMiss = notMiss, miss = wm)
}


updateSSRW <- function(states,yy,zb=rep(0,length(y)),missing,tg,sg){        
  #state-space random walk 
  #update continuous states, random walk
  #missing times, obs y, obs error tg, process error sg
  #zb = z%*%beta
  
  nn <- length(states)
  
  for(t in 1:nn){
    
    VI <- 0
    v  <- 0
    
    if(!t %in% missing){          #observations
      v  <- yy[t]/tg
      VI <- 1/tg
    }
    
    if(t < nn){              #t+1 term excluded for last 
      v  <- v + (states[t+1] - zb[t])/sg 
      VI <- VI + 1/sg
    }
    
    if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + (states[t-1] + zb[t-1])/sg
      VI <- VI + 1/sg
    }
    
    V     <- 1/VI
    states[t] <- rnorm(1,V*v,sqrt(V))
  }
  states
}


updateVariance <- function(yy,mu,s1=1,s2=1){
  
  #y and mu are n by 1
  
  u1 <- s1 + length(yy)/2
  u2 <- s2 + .5*crossprod(yy - mu)
  
  1/rgamma(1,u1,u2) 
}

updateSSbeta <- function(states,z,sg,priorB,priorIVB,addStates=F){  
  #covariates z
  #if addStates = T:  states[t+1] <- states[t] + zbeta[t]
  
  require(mvtnorm)
  
  nt <- length(states)
  xx <- states[-1]
  zz <- z
  if(nrow(zz) == nt)zz <- zz[-nt,]     # last value already clipped?
  
  V  <- solve( crossprod(z[-nt,])/sg + priorIVB )
  
  if(addStates)xx <- xx - states[-nt]
  
  v  <- 1/sg*crossprod(zz,xx)
  t( rmvnorm(1,V%*%v,V) )
}



getMissX <- function(x, y, first=NULL, last=NULL){
  
  prows <- rrows <- numeric(0)
  missX <- missY <- numeric(0)
  missx <- which(is.na(x),arr.ind=T)
  missy <- which(is.na(y))
  
  if(!is.null(first)){        # first, last bound groups
    ngroup <- length(first)
    for(j in 1:ngroup){           #rows to predict and respond
      prows <- c(prows,c( first[j]:(last[j]-1)) )
      rrows <- c(rrows,c( (first[j]+1):last[j] ) )
      jr <- first[j]:last[j]
      missX <- append(missX, list(which(is.na(x[jr,]),arr.ind=T)) )
      missY <- append(missY, list(which(is.na(y[jr])) ))
    }
    missXbyGroup <- which(is.na(x[prows,]),arr.ind=T)
    missYbyGroup <- which(is.na(y[rrows]))
  }
  
  list(prows = prows, rrows = rrows, missx = missx, missy = missy,
       missX = missX, missY = missY,
       missXbyGroup = missXbyGroup, missYbyGroup = missYbyGroup)
}

imputeX <- function(missx, xx, yy, bg, sg, priorx ){
  
  ix <- sort( unique(missx[,2]) )
  IV <- 1/1000
  
  for(i in ix){
    
    wmi <- missx[missx[,2] == i,1]
    
    V <- 1/( bg[i]^2/sg + IV )
    v <- bg[i]*(yy[wmi+1] - xx[wmi,-i,drop=F]%*%bg[drop=F,-i])/sg + 
      priorx[i]*IV
    xx[wmi,i] <- v*V
  }
  xx
}
.updateBeta <- function(x, y, sigma, priorIVB, priorB, 
                        XX=NULL){  # random vector of coefficients
  
  if(is.null(XX))XX <- crossprod(x)
  
  V  <- solve( XX/sigma + priorIVB ) 
  v  <- crossprod(x,y)/sigma + priorIVB%*%priorB
  t( .rMVN(1,t(V%*%v),V) )                     
}


updateZ <- function(ag, bg, x, z, sg, sinv, low, hiw){
  
  nt   <- nrow(z)
  csig <- cov2cor(sg)
  cinv <- solve( csig )
  yp   <- z*0
  
  tmp  <- .sqrtRootMatrix(t(cbind(bg,ag)),sg,DIVIDE=T)
  bcor <- t(tmp[1:p,])
  acor <- t(tmp[c(p,p+1),])
  
  V2cor <- t(acor)%*%cinv%*%acor
  V2sig <- t(ag)%*%sinv%*%ag
  
  for(t in 1:nt){
    
    # first variable on correlation scale
    
    v <- matrix(0,2)
    VI <- sg*0
    
    if(t > 1){
      zcor <-  t(.sqrtRootMatrix(z[t-1,],sg,DIVIDE=T))
      v  <- v + cinv%*%(bcor%*%x[t,] + acor%*%zcor)
      VI <- VI + cinv
    }
    if(t < nt){
      zcor <-  t(.sqrtRootMatrix(z[t+1,],sg,DIVIDE=T))
      v  <- v + t(acor)%*%cinv%*%(zcor - bcor%*%x[t+1,])
      VI <- VI + V2cor
    }
    V <- solve(VI)
    
    mu <- V%*%v
    zcor <-  t(.sqrtRootMatrix(z[t,],sg,DIVIDE=T))
    z[t,1] <- tnorm.mvt(zcor,V%*%v,V,low[t,],
                        hiw[t,],whichSample=1)[1]
    
    v <- matrix(0,2)
    VI <- sg*0
    
    # second variable on variance scale
    if(y[t,2] == 0){
      
      if(t > 1){
        v  <- v + sinv%*%(bg%*%x[t,] + ag%*%z[t-1,])
        VI <- VI + sinv
      }
      if(t < nt){
        v  <- v + t(ag)%*%sinv%*%(z[t+1,] - bg%*%x[t+1,])
        VI <- VI + V2sig
      }
      V  <- solve(VI)
      mu <- V%*%v
      z[t,2] <- tnorm.mvt(z[t,],V%*%v,V,low[t,],
                          hiw[t,],whichSample=2)[2]
    }
    if(t > 1){
      zcor <-  t(.sqrtRootMatrix(z[t-1,],sg,DIVIDE=T))
      yp[t,1] <- rnormCond(acor, bcor, csig, t(x[drop=F,t,]), 
                           zcor, 1) # correlation 
      yp[t,2] <- rnormCond(ag, bg, sg, t(x[drop=F,t,]), 
                           t(z[drop=F,t-1,]), 2)  
    }
  }
  list( z = z, yp = yp )
}

.updateWishart <- function(xx,yy,df,beta=NULL,IXX=NULL,WX=NULL,WIX=NULL,
                                  TRYPRIOR=F){
  #df  <- n - Q + S - 1
  S     <- ncol(yy)
  index <- 0
  XX  <- crossprod(xx)
  IXX <- solve(XX)
  
  D  <- diag(1,nrow(xx)) - xx%*%IXX%*%t(xx)
  SS  <-  t(yy)%*%D%*%yy
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){
    tiny <- 1e-8
    SS[SS < tiny] <- tiny
    message('warning: updateWishartNoPrior')
    SS <- crossprod(yy - xx%*%beta) +  diag(diag(SS)*.001)#*nrow(SS)
    SS <- SS + diag(diag(SS)*.1)
    testv <- try(chol(SS),T)
    index <- 1
  }
  SI <- chol2inv(testv)
  
  z  <- matrix(rnorm(df*S),df,S)%*%chol(SI)
  
  sinv  <- crossprod(z)
  sigma <- solve(sinv)
  list( sigma = sigma, sinv = sinv, indicator = index )
}

####### spatial 


subETOPO5 <- function(lon,lat){
  
  #lon is negative degrees W, positive degrees E
  #lat is positive N
  
  require(geomapdata)
  
  data(ETOPO5)
  
  ilon <- nrow(ETOPO5)/360
  ilat <- ncol(ETOPO5)/180
  
  x <- lon
  y <- lat
  
  
  xx <- 360 + lon
  
  lat[lat > 0] <- 90 - lat[lat > 0]
  lat[lat < 0] <- -lat[lat < 0] + 90
  
  rlon <- xx*ilon
  rlat <- lat*ilat
  
  z <- ETOPO5[rlon[1]:rlon[2],]
  z <- z[,rlat[1]:rlat[2]]
  
  x <- seq(x[1],x[2],length=nrow(z))
  y <- seq(y[1],y[2],length=ncol(z))
  
  list(x = x, y = y, z = z)
  
}


mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  
  if(is.null(scale))scale <- 1
  
  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
  
  par(pin=c(px,py))
  
  invisible( c(px,py) )
  
}


speciesNAmerMap <- function(lon, lat, opt){
  
  ngrid <- 60; trim <- c(.001,.999)
  mfrow <- c(1,1)
  POINTS <- add <- STATES <- F
  newPage <- T 
  yy <- zlevs <- cRamp <- titles <- NULL
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  INT <- F
  
  maplon <- xm <- range(lon)
  maplat <- ym <- range(lat)
  zm <- NULL
  
  dscale <- min( diff(xm),diff(ym) )
  
  mdata <- 'world'
  if(dscale < 10){
    INT <- T
    mdata <- 'state'
  }
  if(!STATES)mdata <- 'world'
  
  if(!is.matrix(yy)){
    yy <- matrix(yy,ncol=1)
    colnames(yy) <- ' '
  }
  scaleCoords <- rbind(maplon[1] + diff(maplon)*c(.9,.97),
                       maplat[1] + diff(maplat)*c(.05,.25))
  
  n <- nrow(yy)
  S <- ncol(yy)
  
  ncol   <- 50
  
  if(is.null(cRamp))cRamp <- c('white','yellow','orange',
                               'red','brown','black')
  
  colF   <- colorRampPalette(cRamp)
  
  if(!add & newPage)par(mfrow=mfrow, bty='n', mar=c(.1,.1,.1,.1) )
  
  for(j in 1:S){
    
    pj <- colnames(yy)[j]
    ssj <- yy[,j]              
    
    maps::map(boundary=T,col='grey',lwd=2,xlim=range(xm),ylim=range(ym))
    
    sc <- zl <- zlevs
    
    if(is.null(zlevs)){
      sc <- seq(min(ssj,na.rm=T),max(ssj,na.rm=T),length=ncol)
      ss <- range(sc)
      zl <- round(seq(ss[1],ss[2],by=diff(ss)/50),4)
    }
    ss <- range(sc)
    zint <- findInterval(sc,zl)
    
    colq   <- colF(length(zl))
    ncc  <- length(colq)
    
    colz  <- colq[zint]
    
    values2contour(xx=lon,yy=lat,
                   z=ssj,nx=ngrid,ny=ngrid,col=colz,
                   zlevs=zl,lwd=4,add=T,fill=T)
    
    zcol <- colq[findInterval(ssj,sc)]
    if(POINTS)symbols(lon,lat,circles=ssj/max(ssj)/2,inches=F,fg=zcol,bg=zcol,add=T)
    
    mapMask(lon,lat,dx=.2*dscale/ngrid,dy=.2*dscale/ngrid,
            whiteOutDist=4*dscale/ngrid,col='white')
    
    map(mdata,add=T,interior=INT,col='white',lwd=5)
    
    map(mdata,add=T,interior=INT,col='black',lwd=1)
    
    if(max(ss) > 1000)ss <- round(ss,-2)
    if(max(ss) > 100)ss <- round(ss,-1)
    if(max(ss) > 10)ss <- round(ss,0)
    if(max(ss) > 1)ss <- round(ss,1)
    if(max(ss) < 1)ss <- signif(ss,2)
    
    if(!is.null(titles))title(titles[j])
    
    colorLegend(scaleCoords[1,],scaleCoords[2,],
                scale=c(1:ncc),cols=colq,text.col='black',
                labside='left',endLabels=ss)
  }
  
}




regMap <- function(xx, yy, opt){
  
  # if IMAGE z is a matrix with length(xx) rows, length(yy) columns)
  # if !IMAGE, xx, yy, and z are all same length
  require(maps)
  require(maptools)
  require(mapdata)
  # require(sp)
  # require(rworldmap)
  
  z <- zminmax <- NULL
  fg <- rep(1,length(x)); bg <- rep(1,length(x)); scaleSym <- 1; mapscale <- 1
  xl <- range(xx); yl <- range(yy); lineCol <- 'black'; lwd  <- 1
  axes <- T
  POINTS <- IMAGE <- ADD <- county <- F
  stateVec <- c('north carolina','south carolina','tennessee',
                'georgia','virginia')
  region <- 'dummy'
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  
  if(!ADD){
    mapSetup(xl,yl,scale=mapscale)
  }
  
  colSea   <- colorRampPalette(c('black','darkblue','blue','lightblue','white'))
  cols <- c(colSea(5),terrain.colors(100))
  
  zr <- max(z)
  bb <- c( -10000,seq(-300,-5,length=5), seq(-2,160,length=50), seq(170,zr+1,length=50))
  
  if(!is.null(zminmax)){
    bb <- seq(zminmax[1],zminmax[2],length=106)
  }
  
  if(IMAGE){
    image(xx,yy,z,col=cols,breaks=bb,add=ADD,xlab=' ',ylab=' ',xlim=xl,ylim=yl) 
    ADD <- T
  }
  
  if(!is.null(stateVec) & county & region != 'andes'){
    map('county',stateVec,boundary=T,col='grey',lwd=lwd,
        add=ADD,xlim=range(xx),ylim=range(y))
    ADD <- T
  }
  if(axes)map.axes()
  
  if(IMAGE)return()
  
  if(POINTS){
    if(is.null(z)) points(xx,yy,col=bg)
    if(!IMAGE & !is.null(z) & length(xx) == length(y)){
      symbols(xx,y,circles=z*scaleSym,fg=fg,bg=bg,inches=F,add=T)
    }
  }
}


values2contour <- function(xx,yy,z,nx=100,ny=100,lty=1,labcex=.7,
                           col='black',lwd=1,zlevs=NULL,add=T,fill=F,
                           drawlabels=F){    
  
  # contours where x,y is not a uniform grid
  
  require(MBA)
  
  xyzmat <- cbind(xx,yy,z)
  
  colnames(xyzmat) <- c('x','y','z')
  
  surf  <- mba.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  # surf  <- myMBA.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  if(is.null(zlevs)){
    zlevs <- signif(seq(min(z), max(z), length=3),1)
  }
  contour(surf, levels=zlevs,lwd=lwd,lty=lty,col=col,add=add,labcex=labcex,
          drawlabels=drawlabels)
  
  if(fill){
    zl <- zlevs
    if(length(zl) == 1)stop('fill.contour() needs at least 2 contour lines')
    .filled.contour(surf$x,surf$y,surf$z,levels=zl,col=col)
  }
  invisible(surf)
}


mapMask <- function(xx,yy,dx=NULL,dy=NULL,whiteOutDist=1,col='white'){    
  
  #mask parts of map with no data, grid density (dx,dy)
  #add to exising map
  #xy is an expand.grid of unmasked pixels
  
  require(RANN)
  
  rx <- range(xx,na.rm=T)
  ry <- range(yy,na.rm=T)
  
  if(is.null(dx)){
    dx <- diff(rx)/50
    dy <- diff(rx)/50
  }
  
  
  xnew <- seq(rx[1],rx[2],by=dx)
  ynew <- seq(ry[1],ry[2],length=length(xnew))
  xxy   <- as.matrix( expand.grid(xnew,ynew) )
  nr   <- nrow(xxy)
  
  nx <- 0
  wd <- whiteOutDist*1.5
  
  while(nx < min(c(5000,nr/4))){
    tmp <- nn2(cbind(xx,yy),xxy,k = 1)$nn.dists
    xy  <- xxy[tmp > wd,]
    nx  <- nrow(xy)
    if(is.null(nx))return()
    wd  <- wd/1.3
  }
  
  symbols(xy[,1],xy[,2],squares=rep(dx,nx),inches=F,fg=col,bg=col,add=T)
  
}


############################################################
colorLegend <- function(xx,yy,ytick=NULL,scale=seq(yy[1],yy[2],length=length(cols)),
                        cols,labside='right',
                        text.col=NULL,
                        bg=NULL,endLabels=NULL,endCols=NULL){  
  
  nn <- length(scale)
  ys <- seq(yy[1],yy[2],length=nn)
  
  for(j in 1:(length(scale)-1)){
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=3)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    cx <- 'black'
    if(!is.null(endCols))cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels,col=cx)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2,col=cx)
  }
}

finleyChains <- function(betaS, beta.tuning, priors){
  
  betaS <- rnorm(length(betaS),betaS,abs(betaS)*.1)
  phiS  <- .tnorm(1,.1,2,.05,1)
  sigS  <- .tnorm(1,0,3,3,1)
  
  out <- spGLM(y[,spec]~x, family="poisson", coords=coords, knots=knots,
               starting=list("beta"=betaS, "phi"= phiS,"sigma.sq"= sigS, "w"=0),
               tuning=list("beta"=beta.tuning, "phi"=0.4,"sigma.sq"=1, "w"=0.1),
               priors=priors,
               amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,
                          "accept.rate"=0.43),
               cov.model="exponential", verbose=F, n.report=500)
  out$DIC <- unlist( spDiag(out, verbose=F) )['DIC4']
  out
}