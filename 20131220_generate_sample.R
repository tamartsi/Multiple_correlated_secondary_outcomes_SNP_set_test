
expit <- function(x) {
	exp(x)/(1 + exp(x))
}

g.cause <- "rs6869352"

generate.sample <- function(G, g.cause, sig.gy = NULL, sig.yd = NULL, sig.gd, n.y = 4, pop.rate.D = 0.01 , type = "ind", sd.phen = c(1,1.5, 1.4, 1.2), cor.phen  = diag(n.y), under.null = T){
	
	
	stopifnot(is.element(type, c("ind",  "YG-D")))
	n.samp <- nrow(G)
	l.or <- log(pop.rate.D/(1- pop.rate.D)) 

	ind.cause <- which(colnames(G) == g.cause)
	beta.y <- rep(0, dim(G)[2])
	if (!is.null(sig.gy)) beta.y[ind.cause] <- sig.gy
	
	if (type == "ind"){
		
		prob.D <- expit(l.or)
		D <- rbinom(dim(G)[1], 1, prob.D)
		
		
		if (under.null) pheno <- rmvnorm(dim(G)[1], sigma = cor.phen) else {
			pheno <- t(apply(G, 1, function(x) rmvnorm(1, mean = rep(sum(x*beta.y), n.y), sigma = cor.phen)))
		}
		pheno <- scale(pheno)
		Y <- pheno %*% diag(sd.phen)
			
		samp <- cbind(G, Y, D)
		colnames(samp) <- c(colnames(G), paste("Y", 1:4, sep = "."), "D")
		return(samp)
	} else { # type = "YG-D
		
		if (under.null) pheno <- rmvnorm(n.samp, sigma = cor.phen) else {
			pheno <- t(apply(G, 1, function(x) rmvnorm(1, mean = rep(sum(x*beta.y), n.y), sigma = cor.phen)))
		}
		pheno <- scale(pheno)
		Y <- pheno %*% diag(sd.phen)

		beta.D.y <- rep(0, n.y)
		beta.D.y[1] <- sig.yd
		
		beta.D.g <- rep(0, dim(G)[2])
		beta.D.g[ind.cause] <- sig.gd
		
		G.list <- lapply(data.frame(t(G)), function(x) x)
		Y.list <- lapply(data.frame(t(Y)), function(x) x)
		
		D <- unlist(mapply(function(A,B){
			p <- expit(l.or + sum(A*beta.D.y) + sum(B*beta.D.g))
			return(rbinom(1,1,p))		
			}, Y.list, G.list))

		samp <- cbind(G, Y, D)
		colnames(samp) <- c(colnames(G), paste("Y", 1:4, sep = "."), "D")
		return(samp)	

		
	}
}
