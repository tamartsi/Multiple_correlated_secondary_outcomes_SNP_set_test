##################################################################################################################################
##############  An example of generating simulating data and testing using the test proposed by Sofer T., Shifano L. et al. (2016). 
############# For questions and complaints contact Tamar Sofer at tsofer@uw.edu. 
##################################################################################################################################

require(CompQuadForm)
require(mvtnorm)

### source various function needed to generate and analyze the data
### function to generate a study population
work.dir <- getwd()  ### set to whatever!
source(file.path(work.dir, "20131220_generate_sample.R" ))
#### the test!
source(file.path(work.dir, "20160609_test_set_mult_out_cc_efficient_Py.R"))
### SMAT function under the null, which we use to estimate the correlation matrix between the outcomes (and not for testing associations with genotypes).
source(file.path(work.dir, "SMAT_null.R"))
### load correlation matrix between the outcomes, and simulation genetic data (simulated using Hapgen)
load(file.path(work.dir,"data_cor.dat"))
load(file.path(work.dir, "simulated_gene_6k.dat"))



G <- as.matrix(t(adtv.dat[,3:60002]))
colnames(G) <- adtv.dat[,1]



n.sim <- 100 # say ;)

### set up simulation parameters
sim.type <-  "YG-D"
sig.yd <- log(2)
sig.yd.string <- "log2"
sig.gd <- log(1.7)
sig.gd.string <- "log_1.7"

#### the effect we want to test:
sig.gy <- 0.14
sig.gy.string <- "0.14"
under.null = FALSE


n.snp <- ncol(G)
p.val <- matrix(NA, nrow = n.sim	, ncol = 3) ## matrix of estimated p-values from simulations
colnames(p.val) <- c("weighted", "unweighted", "controls") 

pval_file <- paste0(work.dir, "/pvals_", sim.type,  "_sigGY_", sig.gy.string, "_sigYD_", sig.yd.string , "_sigGD_", sig.gd.string, "_common.txt")


#### set the baseline probability of disease for people without the causal variant or outcomes of value 0, so that the final population 
#### prevalence of the disease is 8%:

        if ((abs(sig.gd - log(1.7)) < 1e-7)  & (abs(sig.yd - log(2)) < 1e-7)) pop.rate.D <- 0.044
        if ((abs(sig.gd - log(1.7)) < 1e-7)  & abs(sig.yd - log(2)/2) < 1e-7) pop.rate.D <- 0.05
        if (abs(sig.gd - log(1.7)/2)< 1e-7  & abs(sig.yd - log(2)) < 1e-7) pop.rate.D <- 0.055
        if (abs(sig.gd - log(1.7)/2)<1e-7 & abs(sig.yd - log(2)/2) < 1e-7) pop.rate.D <- 0.063        
        if (sig.gd == 0 & sig.yd == 0) pop.rate.D <- 0.08
		if (sig.gd != 0 & sig.yd == 0) quit("n")
		if (sig.gd == 0 & sig.yd != 0) quit("n")
		if (sim.type == "ind" & (sig.yd != 0 | sig.gd != 0)) quit("n")

seed <- 0

for (k in 1:n.sim){
	set.seed(seed + k)


	samp <- generate.sample(G, g.cause = g.cause, sig.gy = sig.gy, sig.yd = sig.yd, sig.gd = sig.gd, cor.phen = cor.mat,  type = 	sim.type, under.null = under.null,  pop.rate.D = pop.rate.D)


	col.D.stat <- ncol(samp)
	
	p.D <- sum(samp[,col.D.stat])/nrow(samp)

	## sampling case-control data:
	n = 500

	if (sum(samp[,col.D.stat]) >= n) {
		
	
	cases <- samp[sample(which(samp[,col.D.stat] == 1), n),]
	controls <- samp[sample(which(samp[,col.D.stat] == 0), n),]

	D <- c(rep(1,n), rep(0,n))
	
	we.cases <- p.D/(nrow(cases)/(nrow(cases) + nrow(controls)))
	we.controls <- (1-p.D)/(nrow(controls)/(nrow(cases) + nrow(controls)))
	we <- c(rep(we.cases, nrow(cases)), rep(we.controls, nrow(controls) ))
	y <- rbind(cases[,88:91], controls[,88:91])
	G.t <- rbind(cases[,1:87], controls[,1:87])   ### test 87 variants, using 1000 individuals. 

	X <- as.matrix(rep(1, nrow(cases) + nrow(controls)))
	

	set.test.null.weighted <- test.set.mult.out.cc.null(y, D, X, pop.rate.D, weights = we)
	set.test.null.unweighted <- test.set.mult.out.cc.null(y, D, X, pop.rate.D, weights = rep(1, length(D)))
	set.test.null.cont <- test.set.mult.out.cc.null(y[which(D == 0), , drop = F], D[which(D==0)], X[which(D==0), , drop = F ], pop.rate.D, weights = rep(1, sum(D == 0)))
	
	p.val[k, "weighted"] <- test.set.mult.out.cc.from.null(set.test.null.weighted, G.t)$p.value
	p.val[k, "unweighted"] <- test.set.mult.out.cc.from.null(set.test.null.unweighted, G.t)$p.value
	p.val[k, "controls"] <- test.set.mult.out.cc.from.null(set.test.null.cont, G.t[which(D == 0), ,drop = F])$p.value
	
	
	}
	

	}
	
	write.table(p.val, col.names = T, row.names = F, file = pval_file, append = T)
	
############   example computing times:
#### for fitting the null model (i.e. prior to testing):
system.time(set.test.null.weighted <- test.set.mult.out.cc.null(y, D, X, pop.rate.D, weights = we))
   # user  system elapsed
  # 0.266   0.000   0.266	
	
system.time( test.set.mult.out.cc.from.null(set.test.null.weighted, G.t)$p.value)
   # user  system elapsed
  # 0.057   0.000   0.058