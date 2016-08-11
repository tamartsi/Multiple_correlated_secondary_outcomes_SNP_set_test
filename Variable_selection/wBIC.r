#-------------------------------------------------------------------------------
#          wBIC.r 
#
# Version: v0.1 (7-12-2016)
# Description: 
#   + Scaled marginal models for multiple continuous (secondary) outcomes. 
#   + Allows for unpenalized covariates and penalized "shared" effects (of interest)
#-------------------------------------------------------------------------------
# y       = response  (n x m)
# x       = covariate (n x p); be sure to include a column of 1's if you want intercepts
# w       = shared exposures  (n x q)
# d       = weight    (n x 1); set to 1 if sample is a random sample from a cohort
# working = "ind","exch","unstr" choices for working correlation structure                                                 
# conv    = covergence threshold
# maxiter = maximum number of iterations
# penalty.fun = penalty function to use on "shared" effects
# lambda  = default is NULL; sequence of penalty parameters to use
# amcp    = default is 3.7; additional parameter for mcp penalty
# npts    = default is NULL; will be length(lambda)
# wt      = pm+q vector specifying which columns of CapX will be penalized (1=penalize, 0=no penalization)
# eps     = default is 10^(-6); used in penalty function
# accuracy= default is 1e-4; set all values below this number to 0
#use.sample.init= logical; if FALSE, fixed starting values (not based on data) are used             
# center  = logical; center design matrix?
# std     = logical; standardize design matrix?
# gamhat0 = default is NULL; can provide initial value of gamma=(beta_1',...,beta_m',alpha)'
# warm    = logical; should we use warm start in path?
#-------------------------------------------------------------------------------


# indices for lower trianglar matrix
lower.tri <- function (x, diag=F) {
	if(!is.matrix(x)) x = as.matrix(x)
    if (diag) row(x) >= col(x)
    else row(x) > col(x)
  }

# get lower trianglar submatrix
lowerTriangle <- function(x, diag=F) {
	x[lower.tri(x,diag=diag)]
  }

# construct block diagonal matrix from vectors
vdiag <- function(v) {
	nr = nrow(v)
	nc = ncol(v)
	mat = matrix(0,nr=nr,nc=nc*nr)
	for(i in 1:nrow(v)) {
		mat[i,(1:nc)+(i-1)*nc] = v[i,]
	  }
	return(mat)
  }

# for MM
pen.mat = function(abs.b,abs.b0, eps, penalty.fun, wt, lambda, a=3.7){
	   if (penalty.fun == "alasso") q.b = q.alasso(abs.b0,lambda, wt) else
	   if (penalty.fun == "mcp")    q.b = q.mcp(abs.b,lambda,wt,a) else
  	 if (penalty.fun=="scad")     q.b = q.scad(abs.b,lambda,wt,a) 
			
		eps.b = eps + abs.b
		return(diag(drop(q.b/eps.b)))
		}

		
q.alasso = function(abs.b0,lambda,wt){
		return(wt*lambda/abs.b0)
  }
		
q.mcp = function(abs.b,lambda,wt,a){
  out= wt*pmax(0,(lambda - abs.b/a))
  return(out)
  }		

q.scad = function(abs.b,lambda,wt,a){
  abeta=abs.b
  pbeta0=(abeta<=lambda)+pmax(a*lambda-abeta,0)*(abeta>lambda)/((a-1)*lambda)

  return(wt*lambda*pbeta0)
  }


 
# main function; one lambda            
wscamar_pen <- function(yy,yx,yw,xx,xw,xy,ww,wx,wy,x,y,w, 
                        sigma2, R,d,working,conv=1e-4,maxiter=5e2,
                        penalty.fun, lambda, amcp=3.7, wt, eps = 10^(-6),
                        accuracy = 1e-10,use.sample.init=TRUE,  
                        center=F, std=F, gamhat0=NULL, sum.d, lM, n, est.sigma2) {

 
	# big while loop for estimation ---------------------------
	max.diff = 1e6
	iter = 1
	
	p <- nrow(xx)
	M <- nrow(yy)
	q <- nrow(ww)
	
  gamhat = rep(0,p*M+q)
  psi = diag(sigma2) 
  lower = rep(1,M)
		
	while(max.diff>conv & iter<maxiter) {

		invR = solve(R)
    
		#------------------
		#  update gamhat 
		#------------------
    gamhat.old = gamhat

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))

		if (is.null(gamhat0)) {
		gamhat0 = try(solve(XinvRX,XinvRy),silent=T); if (is.character(gamhat0)) gamhat0=ginv(XinvRX)%*%XinvRy      }
    if (iter==1) {gamhat = gamhat0; beta0 = matrix(gamhat[1:(M*p)],nr=p); alpha0 = matrix(gamhat[(M*p+1):(M*p+q)],nr=q)}
     
    A = -XinvRX
    
    mm=0
	  b.old = rep(Inf, length(gamhat))
	  gamhat.old = gamhat
	  b.new = gamhat
	  while(mm<maxiter & max(abs(b.old-b.new))>accuracy){
		mm=mm+1
		b.old = b.new
		q.mat = pen.mat(abs(b.old),abs(gamhat0), eps, penalty.fun, wt=wt, lambda, a=amcp)

		inv = try(solve(A - n*q.mat),silent=T); if (is.character(inv)) inv=ginv(A-n*q.mat)
		b.new = b.old -inv %*% (XinvRy - XinvRX%*%b.old - n*q.mat%*%b.old)
		}
		
	  b.new[which(abs(b.new)<accuracy)] = 0
    gamhat = b.new
    beta = matrix(gamhat[1:(M*p)],nr=p)
		alpha = matrix(gamhat[(M*p+1):(M*p+q)],nr=q)
    
		#------------------
		#  update sigma2
		#------------------
    
		sigma2.old = sigma2
		if (est.sigma2){
		for(j in 1:M) {
 			old2=sigma2[j]
			incre=old2
			while(abs(incre)/old2>=conv) {
				old2  = sigma2[j]
				old   = sqrt(sigma2[j])
				bj    = beta[,j]

				upper = yy[j,j]-old*yx[j,]%*%bj-old*yw[j,]%*%alpha-sum.d*old2
				lower[j] = sum.d+t(bj)%*%xx%*%bj/2+t(alpha)%*%wx%*%bj+t(alpha)%*%ww%*%alpha/2
				incre = upper/lower[j]
				sigma2[j] = old2 + incre
			}
		}
    
		#------------------
		#  update R
		#------------------
		if(working=="exch" | working=="unstr") {
			phi = diag(1/sqrt(c(sigma2))) # not psi, updated by new sigma2
  			a = alpha%*%t(lM)
			rr = phi%*%yy%*%phi - phi%*%yx%*%beta - t(beta)%*%t(yx)%*%phi - phi%*%yw%*%a - t(a)%*%t(yw)%*%phi + t(beta)%*%xx%*%beta + t(a)%*%ww%*%a + t(beta)%*%xw%*%a + t(a)%*%wx%*%beta
			sum.rr = sum(lowerTriangle(rr))
			if(working=="exch") {
				theta = sum.rr / (sum.d*M*(M-1)/2) 
				R = matrix(theta,M,M)
				diag(R) = 1
			} else {
				rr = rr/sum.d                                    
				R = rr
				diag(R) = 1
			}
		}
		}
		# compute differences
		max.diff = max(abs(gamhat-gamhat.old)/(abs(gamhat.old)+eps),abs(sigma2-sigma2.old)/sigma2.old)
		iter = iter+1
		
	} # end of big while loop
	# big while loop for estimation completed -----------------
	  
	# invR
	invR = solve(R)
  
	out = vector("list",8)

  names(out) = c("beta","alpha","sigma2","R",
                 "wBIC","wedf","L","flag.w")
	
	w.out = get.wbics(y,x,w,d, invR, sigma2,alpha,beta,beta0,alpha0,n,p,M,lM, eps, penalty.fun, lambda, amcp)
	
  
  out$beta = beta
	out$alpha = alpha
	out$sigma2 = sigma2
	out$R = R
  out$wBIC  = w.out$wBIC
  out$wedf  = w.out$edf
  out$L     = w.out$L
  out$flag.w   = w.out$flag

  
	if(iter>=maxiter) warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.")

  
  return(out)
}
            


get.wbics<- function(y,x,w,d, invR, sigma2,alpha,beta,beta0,alpha0,n,p,M,lM, 
                        eps, penalty.fun, lambda, amcp){

  yy=crossprod(y,y)
  yx=crossprod(y,x)
  yw=crossprod(y,w)
  xx=crossprod(x,x)
  xw=crossprod(x,w)
  xy=t(yx)
  ww=crossprod(w,w)
  wx=t(xw)
  wy=t(yw)


# y.star
y.star = y%*%diag(1/sqrt(sigma2))
# residual matrix at estimated values (nxM)
resid = y.star - x%*%beta - w%*%alpha%*%t(lM)

U.beta = matrix(0,nrow=n,ncol=p*M)

n0dx = which(alpha!=0)
L    = length(n0dx)

resid.invR = resid %*% invR # (nxM)	  
UU = matrix(0,nrow=p*M+L,ncol=p*M+L)
if (L==0) {  # all alpha =0
  for (i in 1:n) {
    Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
    U.beta[i,] = Ui.beta
    Ui = c(Ui.beta)
    UU = UU + tcrossprod(Ui)
  } } else {
    U.alpha = matrix(0,nrow=n,ncol=L)
    for(i in 1:n) {
      Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
      U.beta[i,] = Ui.beta
      Ui.alpha = w[i,n0dx] %*% t(resid.invR[i,]) %*% lM
      U.alpha[i,] = Ui.alpha
      Ui = c(Ui.beta, Ui.alpha)
      UU = UU + tcrossprod(Ui)
    } }

if (L!=0){
  yw=crossprod(y,w[,n0dx])
  xw=crossprod(x,w[,n0dx])
  ww=crossprod(w[,n0dx],w[,n0dx])
  wx=t(xw)
  wy=t(yw)
  
  # hessian matrix on 2-2 corner	
  phi = diag(1/sqrt(c(sigma2))) # not psi
  invR.phi = invR %*% phi
  l.invR.phi = t(lM) %*% invR.phi
  l.invR = t(lM) %*% invR
  invR.l = t(l.invR)
  l.invR.l = l.invR %*% lM
  XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                 cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
  H.gam.gam=XinvRX
  } else { # no w's	
    XinvRX = kronecker(invR,xx)
    H.gam.gam = XinvRX 
    }
H =H.gam.gam

Hinv = try(solve(H),silent=T); if (is.character(Hinv)) Hinv=ginv(H)
tH = t(H)
Htinv = try(solve(tH),silent=T); if (is.character(Htinv)) Htinv=ginv(tH)
cov.mat = Hinv%*%UU%*%Htinv
logdet = log(det(cov.mat))
flag=0
if (logdet=="NaN") flag=1

  WW = H.gam.gam 
  WWinv = try(solve(WW),silent=T) 
  if (is.character(WWinv)) WWinv = ginv(WW)
  edf = sum(diag(WWinv%*%UU))
  wBIC=calc.wBIC(resid, p, n, M, L, edf,invR,sigma2)
  
  out = vector("list",4)
  names(out) = c("wBIC","edf","L","flag")
  
  
  out$wBIC      = wBIC
  out$edf       = edf
  out$L         = L
  out$flag      = flag
return(out)
}


get.wbics.null<- function(y,x,d, invR, sigma2,beta,n,p,M,lM, eps, penalty.fun, lambda, amcp){
  
  
    yy=crossprod(y,y)
    yx=crossprod(y,x)
    xx=crossprod(x,x)
    xy=t(yx)
  
  
  # y.star
  y.star = y%*%diag(1/sqrt(sigma2))
  # residual matrix at estimated values (nxM)
  resid = y.star - x%*%beta 
  
  
  U.beta = matrix(0,nrow=n,ncol=p*M)
  resid.invR = resid %*% invR # (nxM)    
  UU = matrix(0,nrow=p*M,ncol=p*M)
  for (i in 1:n) {
    Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
    U.beta[i,] = Ui.beta
    Ui = c(Ui.beta)
    UU = UU + tcrossprod(Ui)
  } 
  
  # no w's	
  XinvRX = kronecker(invR,xx)
  H.gam.gam = XinvRX 
  H = H.gam.gam
    
  # covariance matrix
  Hinv = try(solve(H),silent=T); if (is.character(Hinv)) Hinv=ginv(H)
  tH = t(H)
  Htinv = try(solve(tH),silent=T); if (is.character(Htinv)) Htinv=ginv(tH)
  cov.mat = Hinv%*%UU%*%Htinv

  logdet = log(det(cov.mat))
  flag=0
  if (logdet=="NaN") flag=1
  
  L=0
  WW = H.gam.gam 
  WWinv = try(solve(WW),silent=T) 
  if (is.character(WWinv)) WWinv = ginv(WW)
  edf = sum(diag(WWinv%*%UU))
  wBIC=calc.wBIC(resid, p, n, M, L, edf,invR,sigma2)
  
  out = vector("list",4)
  names(out) = c("wBIC","edf","L","flag")
  

  out$wBIC       = wBIC
  out$edf        = edf
  out$L           = L
  out$flag        = flag
  
  return(out)
}



wscamar_pen_null <- function(yy,yx,xx,xy,x,y, 
                        sigma2, R,d,working,conv=1e-4,maxiter=5e2, eps=1e-6,
                        use.sample.init=TRUE, center=F, std=F, gamhat0=NULL, 
                        sum.d, lM, n, est.sigma2) {
  
  
  # big while loop for estimation ---------------------------
  max.diff = 1e6
  iter = 1
  
  p <- nrow(xx)
  M <- nrow(yy)
  
  gamhat = rep(0,p*M)
  psi = diag(sigma2) 
  lower = rep(1,M)
  
  while(max.diff>conv & iter<maxiter) {
    
    invR = solve(R)
    #------------------
    #  update gamhat 
    #------------------
    gamhat.old = gamhat
    
    phi = diag(1/sqrt(c(sigma2))) # not psi
    invR.phi = invR %*% phi
    l.invR.phi = t(lM) %*% invR.phi
    l.invR = t(lM) %*% invR
    invR.l = t(l.invR)
    l.invR.l = l.invR %*% lM
    XinvRX = kronecker(invR,xx)
                   
    XinvRy = matrix(t(invR.phi%*%yx),ncol=1)
    
    gamhat = try(solve(XinvRX,XinvRy),silent=T); if (is.character(gamhat)) gamhat=ginv(XinvRX)%*%XinvRy      
    beta = matrix(gamhat,nr=p)    
    #------------------
    #  update sigma2
    #------------------
    
    sigma2.old = sigma2
    if (est.sigma2){
      for(j in 1:M) {
        old2=sigma2[j]
        incre=old2
        while(abs(incre)/old2>=conv) {
          old2  = sigma2[j]
          old   = sqrt(sigma2[j])
          bj    = beta[,j]
          
          upper = yy[j,j]-old*yx[j,]%*%bj-sum.d*old2
          lower[j] = sum.d+t(bj)%*%xx%*%bj/2
          incre = upper/lower[j]
          sigma2[j] = old2 + incre
        }
      }
      
      #------------------
      #  update R
      #------------------
      if(working=="exch" | working=="unstr") {
        phi = diag(1/sqrt(c(sigma2))) # not psi, updated by new sigma2
        rr = phi%*%yy%*%phi - phi%*%yx%*%beta - t(beta)%*%t(yx)%*%phi + t(beta)%*%xx%*%beta 
        sum.rr = sum(lowerTriangle(rr))
        if(working=="exch") {
          theta = sum.rr / (sum.d*M*(M-1)/2) # double check this
          R = matrix(theta,M,M)
          diag(R) = 1
        } else {
          rr = rr/sum.d
          R = rr
          diag(R) = 1
        }
      }
    }
    # compute differences
    max.diff = max(abs(gamhat-gamhat.old)/(abs(gamhat.old)+eps),abs(sigma2-sigma2.old)/sigma2.old)
    iter = iter+1
    
  } # end of big while loop
  # big while loop for estimation completed -----------------
  
  # invR
  invR = solve(R)
  
  
  out = vector("list",16)
  names(out) = c("beta","sigma2","R",
                 "wBIC","L","wedf","flag.w")
  
  w.out = get.wbics.null(y,x,d, invR, sigma2,beta,n,p,M,lM,eps, penalty.fun, lambda, amcp)
    
  out$beta = beta
  out$sigma2 = sigma2
  out$R = R
  out$wBIC  = w.out$wBIC
  out$wedf  = w.out$edf2
  out$L      = w.out$L
  out$flag.w = w.out$flag
 
  
  if(iter>=maxiter) warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.")
  
  
  return(out)
}

	

calc.wBIC <- function(resid, p, n, M, L,edf,invR,sigma2){
  P <- edf
  arg1 <- 0
  for (i in 1:n){
    resi <- resid[i,]
    arg1 <- arg1 + t(resi) %*% invR %*% resi
  }
  covar <- solve(invR)
  arg0 <- log(det(covar))
  
  arg1 <- arg1/n
  arg2 <- P*log(n)/n
  return(arg0+ arg1 + arg2)                                                            
}		

                                                                                                            
                           
get.lambda.seq <- function(yy,yx,yw,xx,xw,xy,ww,wx,wy, sigma2, R ,penalty.fun, n.pts, gamhat0,lM, n) {


		invR = solve(R)

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
    l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
	
    XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))                   
                       
 
    lambda.max = max(abs(rev(XinvRy)[1:q]))/n+1
    if (penalty.fun=="alasso"){
    alambda.max= lambda.max*max(abs(rev(gamhat0)[1:q]))
    lambda.seq <- exp(seq(-log(1e5),log(alambda.max+2),length=n.pts))  }else {
    lambda.seq <- exp(seq(-log(1e5),log(lambda.max+20),length=n.pts)) }

    return(lambda.seq)
    }

   
                        
get_gamhat0 <- function(yy,yx,yw,xx,xw,xy,ww,wx,wy, sigma2, R,lM) {
  
	
		invR = solve(R)
		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))

		gamhat0 = try(solve(XinvRX,XinvRy),silent=T); if (is.character(gamhat0)) gamhat0=ginv(XinvRX)%*%XinvRy  
		return(gamhat0)
		}


## path function for multiple lambdas
fit.path <- function(y,x,w,d=NULL,working,conv=1e-4,maxiter=5e2,penalty.fun, 
                             lambda=NULL, amcp=3.7, npts=NULL, wt, eps = 10^(-6),
                             accuracy = 1e-4,use.sample.init=TRUE, 
                             center=FALSE, std=FALSE, gamhat0=NULL, warm=TRUE) {

  x.orig <- x
  w.orig <- w
  y.orig <- y
  
	# check input types and dimensions
	n = nrow(y)
	M = ncol(y)
	p = ncol(x)
  q = ncol(w)

	is.wt = F # consider weighting?
	if(!is.null(d)) { is.wt = T; sum.d = sum(d) } else { d = rep(1,n); sum.d = n }

  if (center){
    if (sd(x[,1]) == 0) {x[,-1] <- scale(x[,-1],T,F)} else {x <- scale(x,T,F)}    
		w <- scale(w,T,F)
	}
	
	if (std){
    sds.x <- apply(x, 2, function(xx) sd(xx)/sqrt(nrow(x)))
		sds.w <- apply(w, 2, function(x) sd(x)/sqrt(nrow(w)))
       if (sds.x[1] != 0)  x[,1] <- x[,1]/sds.x[1]  
       for (jj in 1:ncol(w)) w[, jj] <- w[, jj]/sds.w[jj]
       for (jj in 2:ncol(x)) x[, jj] <- x[, jj]/sds.x[jj]
	}
		
	# precomputation
	if(is.wt) {
		y = diag(sqrt(d))%*%y
		x = diag(sqrt(d))%*%x
		w = diag(sqrt(d))%*%w	
	}
	yy=crossprod(y,y)
	yx=crossprod(y,x)
	yw=crossprod(y,w)
	xx=crossprod(x,x)
	xw=crossprod(x,w)
	xy=t(yx)
	ww=crossprod(w,w)
	wx=t(xw)
	wy=t(yw)
	
	
	# initial values of parameters
  if (use.sample.init) {sigma2 = apply(y,2,sd)^2} else {sigma2 = rep(1,M)}
	psi = diag(sigma2)
	eyeM = diag(M)
	lM = matrix(1,ncol=1,nrow=M)
	if(working=="ind") { 
		R=eyeM 
	} else if(working=="exch") { 
		if (use.sample.init) {theta=mean(lowerTriangle(cor(y)))} else {theta=.3}
		R=matrix(theta,M,M)
		diag(R)=1
	} else { 
		if (use.sample.init) {R=cor(y)} else {theta=.3; R=matrix(theta,M,M); diag(R)=1}  
	} # unstructured, sample correlation

  sig2 = sigma2
	Rmat = R
	

  if (is.null(gamhat0)) gamhat0= get_gamhat0(yy,yx,yw,xx,xw,xy,ww,wx,wy, sigma2, R, lM)
  if (is.null(npts)) npts=25                 
	if (is.null(lambda)) {  lambda_path=get.lambda.seq(yy,yx,yw,xx,xw,xy,ww,wx,wy,sigma2, R ,penalty.fun, npts, gamhat0,lM, n) }
    else {lambda_path=lambda; npts = length(lambda_path)}    
     
	
	# initialize output                                                               
	beta = array(0,dim=c(npts,p,M))                                                       
	alpha = matrix(0, npts, q)
	sigma2 = matrix(0,npts, M)
	R = array(0,dim=c(npts,M,M))
	
	for (j in 1:npts){
  if (j==1){
  out = wscamar_pen(yy,yx,yw,xx,xw,xy,ww,wx,wy,x,y,w, 
                    sigma2=sig2, R=Rmat, d=d,working=working,conv=conv,maxiter=maxiter,
                    penalty.fun=penalty.fun, lambda=lambda_path[j], amcp=amcp, wt=wt, eps=eps,
                    accuracy=accuracy, gamhat0=gamhat0, sum.d=sum.d,lM=lM, 
                    n = n, est.sigma2=TRUE)
  } else {
  out = wscamar_pen(yy,yx,yw,xx,xw,xy,ww,wx,wy,x,y,w, 
                    sigma2=sig2, R=Rmat, d=d,working=working,conv=conv,maxiter=maxiter,
                    penalty.fun=penalty.fun, lambda=lambda_path[j], amcp=amcp, wt=wt, eps=eps,
                    accuracy=accuracy, gamhat0=gamhat0, sum.d=sum.d,lM=lM, 
                    n = n, est.sigma2=FALSE)
  }
  beta[j,,]  = out$beta                                                          
	alpha[j,] = out$alpha
	sigma2[j,] = out$sigma2
	R[j,,] = out$R
  if (warm){
  gamhat0 = c(as.vector(beta[j,,]),alpha[j,])
  sig2    = as.vector(sigma2[j,])
  Rmat    = as.matrix(R[j,,])
  } 
  }


  if (std) {
  for (j in 1:npts){ 
			 if (sds.x[1] != 0) {
			 	beta[j,1,] <- beta[j,1,]/sds.x[1]
			 	}
  			  if (p>2){
          beta[j,2:p,] <- diag(1/sds.x[2:ncol(x)]) %*% beta[j,2:p,] }  else {
          beta[j,2,]   <- beta[j,2,]/sds.x[2]} 
           
  		    alpha[j,] <- alpha[j,]/sds.w
			}
	}
  	
		if (center) {  
                    
			if (sd(x.orig[,1]) == 0)				
				{
        for (j in 1:npts){ 
        intercepts <- beta[j,1,] - t(colMeans(x.orig[,-1]) %*% beta[j,2:ncol(x),]) - rep(colMeans(w.orig)%*%alpha[j,],M)
        beta[j,1,] <- intercepts} 
			}	
 }
 
out = vector("list",5)
names(out) = c("beta","alpha","sigma2","R","lambda")

  out$beta = beta
	out$alpha = alpha
	out$sigma2 = sigma2
	out$R = R
  out$lambda     =lambda_path
  
  return(out)
}


# called after fit.path; wind is output from fit.path
get.BIC= function(y,x,w,d,wind,working,conv=1e-4,maxiter=5e2,use.sample.init=TRUE, 
                   center=FALSE, std=FALSE){
  x.orig <- x
  w.orig <- w
  y.orig <- y
  
  # check input types and dimensions
  n = nrow(y)
  M = ncol(y)
  p = ncol(x)
  q = ncol(w)
  
  is.wt = F # consider weighting?
  if(!is.null(d)) { is.wt = T; sum.d = sum(d) } else { d = rep(1,n); sum.d = n }
  
  if (center){
    if (sd(x[,1]) == 0) {x[,-1] <- scale(x[,-1],T,F)} else {x <- scale(x,T,F)}    
    w <- scale(w,T,F)
  }
  
  if (std){
    sds.x <- apply(x, 2, function(xx) sd(xx)/sqrt(nrow(x)))
    sds.w <- apply(w, 2, function(x) sd(x)/sqrt(nrow(w)))
    if (sds.x[1] != 0)  x[,1] <- x[,1]/sds.x[1]  
    for (jj in 1:ncol(w)) w[, jj] <- w[, jj]/sds.w[jj]
    for (jj in 2:ncol(x)) x[, jj] <- x[, jj]/sds.x[jj]
  }
  
  
  # precomputation
  if(is.wt) {
    y = diag(sqrt(d))%*%y
    x = diag(sqrt(d))%*%x
    w = diag(sqrt(d))%*%w  
  }
  yy=crossprod(y,y)
  yx=crossprod(y,x)
  xx=crossprod(x,x)
  xy=t(yx)  
  
  # initial values of parameters
  if (use.sample.init) {sigma2 = apply(y,2,sd)^2} else {sigma2 = wind$sigma2[1,]}
  psi = diag(sigma2)
  eyeM = diag(M)
  lM = matrix(1,ncol=1,nrow=M)
  if(working=="ind") { 
    R=eyeM
  } else if(working=="exch") { 
    if (use.sample.init) {theta=mean(lowerTriangle(cor(y)))} else {theta=mean(lowerTriangle(wind$R[1,,]))}
    R=matrix(theta,M,M)
    diag(R)=1 
  } else { 
    if (use.sample.init) {R=cor(y)} else {R=wind$R[1,,]}  
  } # unstructured, sample correlation
  
  sig2 = sigma2
  Rmat = R
  gamhat0 = c(as.vector(wind$beta[1,,]),wind$alpha[1,])
  
  # compute BIC on original scale 
  npts <- length(wind$lambda)
  wBIC <- rep(Inf,npts)
  flag.w <- rep(0,npts)
  beta = array(0,dim=c(npts,p,M))                                                       
  alpha = matrix(0, npts, q)
  sigma2 = matrix(0,npts, M)
  R = array(0,dim=c(npts,M,M))
  
  a.n0 = apply(wind$alpha,1,function(x) which(x!=0))
  wdx.prev = 0
  for (j in 1:npts){
    wdx = a.n0[[j]] 
    wdx.cur = wdx
    if (!identical(wdx.cur,wdx.prev)){
      if(length(wdx)>0) {    
        yw=crossprod(y,w[,wdx])
        xw=crossprod(x,w[,wdx])
        ww=crossprod(w[,wdx],w[,wdx])
        wx=t(xw)
        wy=t(yw)
        gm0 = gamhat0[c(1:(M*p),wdx)]
        nowt = rep(0,length(gm0))
        out = wscamar_pen(yy,yx,yw,xx,xw,xy,ww,wx,wy,x,y,as.matrix(w[,wdx]), 
                          sigma2=sig2, R=Rmat, d=d,working=working,
                          conv=conv,maxiter=maxiter,penalty.fun="mcp", 
                          lambda=0, amcp=3.7, wt=nowt, eps=1e-6,
                          accuracy=1e-4, 
                          gamhat0=gm0, sum.d=sum.d,lM=lM, n = n, est.sigma2=use.sample.init)
        alpha[j,wdx] = out$alpha
        } else { # wdx=0
          out = wscamar_pen_null(yy,yx,xx,xy,x,y, 
                                 sigma2=sig2, R=Rmat, d=d,working=working,
                                 conv=conv,maxiter=maxiter,
                                 sum.d=sum.d,lM=lM, n = n, est.sigma2=use.sample.init)
          alpha[j,] = rep(0,q)
          } # end wdx=0
        beta[j,,]  = out$beta                                                              
        sigma2[j,] = out$sigma2
        R[j,,] = out$R
        wBIC[j] = out$wBIC
        flag.w[j] = out$flag.w
        wdx.prev = wdx.cur
    } else { #wdx.cur=wdx.prev
      wBIC[j] = wBIC[j-1]
      beta[j,,]  = beta[j-1,,]                                                         
      alpha[j,] = alpha[j-1,]
      sigma2[j,] = sigma2[j-1,]
      R[j,,] = R[j-1,,]
      wdx.prev = wdx.cur  
    }
  }
  
  if (std) {
    for (j in 1:npts){ 
      if (sds.x[1] != 0) {
        beta[j,1,] <- beta[j,1,]/sds.x[1]
      }
      if (p>2){
        beta[j,2:p,] <- diag(1/sds.x[2:ncol(x)]) %*% beta[j,2:p,] }  else {
          beta[j,2,]   <- beta[j,2,]/sds.x[2]} 
      alpha[j,] <- alpha[j,]/sds.w
    }
  }
  
  if (center) {  
    
    if (sd(x.orig[,1]) == 0)				
    {
      for (j in 1:npts){ 
        intercepts <- beta[j,1,] - t(colMeans(x.orig[,-1]) %*% beta[j,2:ncol(x),]) - rep(colMeans(w.orig)%*%alpha[j,],M)
        beta[j,1,] <- intercepts} 
    }	
  }
  
   
  out = vector("list",6)
  names(out) = c("beta","alpha","sigma2","R",
                 "wBIC","flag.w")
  
  out$beta = beta
  out$alpha = alpha
  out$sigma2 = sigma2
  out$R = R
  out$wBIC       = wBIC
  out$flag.w     = flag.w
 
   return(out)
}






