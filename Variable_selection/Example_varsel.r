# ================================================ #
# Example of variable selection using weighted BIC #
# ================================================ #
# setwd("")  ## if need to set working directory 
source("wBIC.r")


# This example dataset contains 1000 cases and 1000 controls, 
# where the true disease prevalence, pi.D=0.08, with
# m = 4 secondary outcomes (first 4 columns)
# p = 3 = 2 covariates + 1 intercept, (next three columns)
# q = 39 SNPs, (next 39 columns).
# The final column indicates case=1/control=0 status; 
# controls are in rows 1-1000, cases are rows 1001-2000

dat<- read.csv("Example_varsel.csv")
y <- as.matrix(dat[,1:4])
x <- as.matrix(dat[,5:7])
g <- as.matrix(dat[,8:46])
status <- as.vector(dat[,47])

m <- ncol(y)
p <- ncol(x)
q <- ncol(g)

# define weight based on disease prevalence, pi.D, and proportion of cases, p.d
n <- nrow(dat)
n1 <- sum(status)
n0 <- n-n1
p.d <- n1/(n1+n0)
pi.D <- 0.08
wco <- (1-pi.D)/(1-p.d)
wca <- pi.D/p.d
ipw <- as.vector(c(rep(wco,n0),rep(wca,n1)))

# indicate with wt that only SNP effects should be penalized (not covariate effects)
wt <- rep(0,p*m+q)
wt[(p*m+1):(p*m+q)] <- 1

# fit penalized estimation path with mcp along a lambda sequence;
# details on argument values are given in wBIC.r
lambda.seq <- c(0,exp(seq(log(1e-4),log(15),length=99))) 
path.out=fit.path(y=y,x=x,w=g,d=ipw,working="unstr",maxiter=5000,penalty.fun="mcp", 
                lambda=lambda.seq, wt=wt, 
                accuracy = 1e-4,use.sample.init=TRUE, 
                center=TRUE, std=TRUE, gamhat0=NULL) 

# use get.BIC in conjunction with fit.path output to get BICs on re-evaluated data
bics.out = get.BIC(y,x,g,ipw,path.out,working="unstr",maxiter=5000,use.sample.init=FALSE)

# find index with lowest BIC
idx <- which.min(bics.out$wBIC)

# SNP effects for select index
bics.out$alpha[idx,]
