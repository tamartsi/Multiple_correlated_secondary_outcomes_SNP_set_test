
require(CompQuadForm)




test.set.mult.out.cc.null <-  function(y, D, X, p.D, weights = NULL, working = "unstr", common = T, sandwich = F, timing = F){
	## y is an n x m matrix of outcomes
	## D is case control status
	## X is an n x p matrix of covariates
	## p.D is disease prevalence in the population\
	## working - type of  working correlation structure
	## common: test using a common effect assumption of the set on the outcomes? 
	## sandwich - the sandwich form of the correlation structure should be used?
	## timing: for the wscamar null fit
	
	if (is.null(weights)) {
		p.d1 <- p.D/(sum(D)/length(D))
		p.d0 <- (1-p.D)/((length(D) - sum(D))/length(D))
		weights <- rep(p.d1, length(D))
		weights[which(D == 0)] <- p.d0
	
		}
	
	null.fit <- SMAT.null(y, X, d = weights, working = working, timing = timing)
	y.star <- y %*% diag(null.fit$sigma2^(-1/2))
	inv.cor <- solve(null.fit$R)
	
	weight.mat <- diag(weights)
	
	 P <-  weight.mat -   tcrossprod(crossprod(weight.mat, tcrossprod(tcrossprod(X , chol2inv(chol(crossprod(X,crossprod(weight.mat ,X))))) ,X)), weight.mat)
	Py.1 <- P %*% y.star

    Py <- P %*% y.star %*% inv.cor

	return(list(P= P, Py = Py, inv.cor = inv.cor))




}



test.set.mult.out.cc.from.null <-  function(null.model, G,  common = T){
	## n.out is the number of outcomes
	## G is an n x d matrix of things to test (i.e. genotyping data)
	## M.sqInvR and Py are outputs of the function test.set.mult.out.cc.prepare
	## common: test using a common effect assumption of the set on the outcomes? 
	
	Py <- null.model$Py
	P <- null.model$P
	inv.cor <- null.model$inv.cor
	
	m <- ncol(Py)
	
	
	 if (common == T) {
	 	Q.a <- crossprod(Py , crossprod(tcrossprod(G), Py ))
		Q <- sum(Q.a)	
		K <- t(G) %*% P %*% P %*% G*sum(inv.cor)
	

	 	}
		
		## Q is the test statistics
		## K is the matrix of which we take the eigenvalues for the coefficients of the chi square variables in the sum given to davis method. 
	
	 pval <- Get_PValue(K, Q)$p.value
		
	return(list(p.value = pval, Q= Q))
}











Get_PValue<-function(K,Q){
	
	lambda<-Get_Lambda(K)
	#print(lambda)
	n1<-length(Q)

	p.val<-rep(0,n1)
	p.val.liu<-rep(0,n1)
	is_converge<-rep(0,n1)

	for(i in 1:n1){
		out<-davies(Q[i],lambda = lambda, acc=10^(-6))

		p.val[i]<-out$Qq
		p.val.liu[i]<-liu(Q[i],lambda)

		is_converge[i]<-1
		
		# check convergence
		if(length(lambda) == 1){
			p.val[i]<-p.val.liu[i]
		} else if(out$ifault != 0){
			is_converge[i]<-0
		}
	
		# check p-value
		if(p.val[i] > 1 || p.val[i] < 0 ){
			is_converge[i]<-0
			p.val[i]<-p.val.liu[i]
		}
	}

	return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge))

}



Get_Lambda<-function(K){
	
	if (sum(K != 0 ) == dim(K)[1]) {
		lambda1 <- diag(K)
		
	} else{
	out.s<-eigen(K,symmetric=TRUE)
	
	lambda1<-out.s$values}
	IDX1<-which(lambda1 >= 0)

	# eigenvalue bigger than sum(eigenvalues)/1000
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)

	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	}
	lambda<-lambda1[IDX2]
	lambda

}


