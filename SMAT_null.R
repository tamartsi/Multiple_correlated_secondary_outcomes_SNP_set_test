#-------------------------------------------------------------------------------
#   Liz's code modified - under the null, no testing. 
#          WSCAMAR.R 
#
# Authors: Lin Li (linli@hsph.harvard.edu, ll323@cornell.edu)
#          Elizabeth D. Schifano (eschifan@hsph.harvard.edu)
# Version: v0.9 (07-09-2011)
# Description: 
#   + Scaled marginal models for multiple continuous outcomes. 
#   + Features efficient implementation for large-scale problems.
#   + Based on Dr. Elizabeth D. Schifano's work and R scripts as well as Dr. 
#     Jason Roy's paper (Roy, Lin and Ryan, 2003) and SAS scripts.
#-------------------------------------------------------------------------------
# y       = response  (n x M)
# x       = covariate (n x p)
# d       = weight    (n x 1)
# working = "ind","exch","unstr"
# conv    = covergence threshold
# maxiter = maximum number of iterations
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

# main function
SMAT.null <- function(y,x,d=NULL,working,conv=1e-4,maxiter=5e2,timing=T) {

  # timing
  if(timing) time.0=proc.time()
  
	# check input types and dimensions
	##TODO##
	n = nrow(y)
	M = ncol(y)
	p = ncol(x)

	is.wt = F # consider weighting?
	if(!is.null(d)) { is.wt = T; sum.d = sum(d) } else { d = rep(1,n); sum.d = n }

	# precomputation
	if(is.wt) {
		y = diag(sqrt(d))%*%y
		x = diag(sqrt(d))%*%x
	}
	yy=crossprod(y,y)
	yx=crossprod(y,x)
	xx=crossprod(x,x)
	xy=t(yx)

	# initial values of parameters
	betahat = rep(0,p*M)
	sigma2 = apply(y,2,sd)^2 # sample variance
	psi = diag(sigma2)
	eyeM = diag(M)
	lM = matrix(1,ncol=1,nrow=M)
	if(working=="ind") { 
		R=eyeM 
	} else if(working=="exch") { 
		theta=mean(lowerTriangle(cor(y)))
		R=matrix(theta,M,M)
		diag(R)=1 
	} else { 
		R=cor(y) 
	} # unstructured, sample correlation

	lower = rep(1,M)

	# big while loop for estimation ---------------------------
	max.diff = 1
	iter = 1
	while(max.diff>conv & iter<maxiter) {

		invR = solve(R)

		#------------------
		#  update betahat 
		#------------------
		betahat.old = betahat

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = kronecker(invR,xx)
		XinvRy = matrix(t(invR.phi%*%yx),ncol=1)

		 betahat = solve(XinvRX,XinvRy)
		beta.mat = matrix(betahat, nrow = p)
		
		#------------------
		#  update sigma2
		#------------------
		sigma2.old = sigma2

		for(j in 1:M) {
 			old2=sigma2[j]
			incre=old2
			while(abs(incre)/old2>=conv) {
				old2  = sigma2[j]
				old   = sqrt(sigma2[j])
				bj    = beta.mat[,j]
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
			rr = phi %*% yy %*% phi - phi %*% yx %*% beta.mat - t(beta.mat) %*% t(yx) %*% phi + t(beta.mat) %*% xx %*% beta.mat
			sum.rr = sum(lowerTriangle(rr))
			if(working=="exch") {
				theta = sum.rr / (n*M*(M-1)/2 - p*M - M) # double check this
				R = matrix(theta,M,M)
				diag(R) = 1
			} else {
				rr = rr/n
				R = rr
				diag(R) = 1
			}
		}
		
		# compute differences
		max.diff = max(abs(betahat-betahat.old)/abs(betahat.old),abs(sigma2-sigma2.old)/sigma2.old)
		iter = iter+1
		
	} # end of big while loop
	# big while loop for estimation completed -----------------

  if(timing) time.fit = proc.time()
  

  
   
	# output a list
	out = vector("list",4)
	names(out) = c("beta","sigma2","R", "resids")
	out$beta = betahat
	out$sigma2 = sigma2
	out$R = R
	y.star <- t(apply(y, 1, function(x) x/sqrt(sigma2)))
	out$resids <- y.star - x%*% beta.mat

	if(iter>=maxiter) warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.")

  if(timing) {
    running.time = time.fit-time.0
    print(running.time)
  }
  
  return(out)
}


