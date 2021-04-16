#!/usr/bin/env Rscript

dir.create(file.path(Sys.getenv("R_LIBS_USER")),recursive = TRUE, showWarnings = FALSE)
if(!require(optparse, quietly=T)){install.packages('optparse',repos="http://cran.us.r-project.org"); require(optparse, quietly=T);} 

maf1 = 0.3
maf2 = 0.3
# main effects of SNPs
eff1 = c(0,0,1)
eff2 = c(0,0,1)

n_case = 1000
n_control = 1000

# coefficients for determining P(case)
# this is baseline probability, main effect snp1, main effect snp2, interaction effect
coeff = c(0,0,0,0.5)


opt_list <- list(
 	make_option("--maf1", type="numeric", default=0.3, help="MAF of SNP 1 (default = %default)"),
 	make_option("--maf2", type="numeric", default=0.3, help="MAF of SNP 2 (default = %default)"),
 	make_option("--eff1", type="character", action="store", default="0,0,1", help="Biological action of SNP 1 (default = %default)"),
 	make_option("--eff2", type="character", action="store", default="0,0,1", help="Biological action of SNP 2 (default = %default)"),
 	make_option("--case", type="integer", default=1000, help="Number of cases to generate (default = %default)"),
 	make_option("--control", type="integer", default=1000, help="Number of controls to generate (default = %default)"),
 	make_option("--model", type="character", default="0.5,0,0,0", help="Model to generate (default = %default)"),
 	make_option("--penetrance", type="character", default=NULL, help="If given, specify the penetrance table directly (9 numbers separated by commas)"),
 	make_option("--penbase", type="numeric", default=0.25, help="Baseline to use for penetrance tables"),
 	make_option("--pendiff", type="numeric", default=NA, help="Difference between min and max probabilities in the penetrance table (default = 1-2*penbase)"),
 	make_option("--quant", type="integer", default=NULL, help="If given, generate quantitative traits (integer sample size)"),
 	make_option("--snr", type="numeric", default=NULL, help="Percent of variance explained between null + full model (quant only)"),
	make_option("--seed", type="numeric", default=NULL, help="RNG seed")
)


opts <- parse_args(OptionParser(option_list=opt_list))

if(!is.null(opts$seed)){
set.seed(opts$seed)
}

# parse the string values
eff1 <- as.numeric(strsplit(opts$eff1, ",")[[1]])
eff2 <- as.numeric(strsplit(opts$eff2, ",")[[1]])
coeff <- as.numeric(strsplit(opts$model, ",")[[1]])
maf1 <- opts$maf1
maf2 <- opts$maf2
n_case <- opts$case
n_control <- opts$control
# n_samp is the number of quantitative samples
n_samp <- n_case + n_control

penbl <- opts$penbase
pendf <- 1-2*penbl
if(!is.na(opts$pendiff)){
	pendf <- opts$pendiff
}
pentable <- matrix(0,nrow=3,ncol=3)

CC <- TRUE
pen <- FALSE

if(length(coeff) == 3){
# auto-scale the coefficients
	coeff=(coeff-min(coeff))/sum(coeff)
	coeff=c(penbl, pendf * coeff)

} else if(length(coeff) != 4){
	stop("Coefficients in model must be 3 or 4 elements")
}

# Do quantitative traits
if(!is.null(opts$quant)){
	CC <- FALSE
	n_samp <- opts$quant
}

# Specify the penetrance table directly
if(!is.null(opts$penetrance)){
	pen <- TRUE
	pentable <- t(matrix(as.numeric(strsplit(opts$penetrance, ",")[[1]]),nrow=3,ncol=3))
	
	if(min(pentable) != 0 || max(pentable) != 1){
		warning("Penetrance table not scaled from 0 to 1, scaling")
		pentable = (pentable - min(pentable)) / (max(pentable) - min(pentable))
	}
	
	if(CC){
		if(opts$penbase < 0 || opts$penbase >= 1){
			stop("Penetrance baseline not in range [0,1), exiting")
		}
		
		if(pendf < 0 || penbl + pendf > 1){
			stop("Penetrance baseline and difference have nonsenical values (diff < 0 or baseline + diff > 1)")
		}

	}
	
	pentable = penbl + pendf * pentable;

}

# sanity check some inputs
if(opts$maf1 < 0 || opts$maf1 > 1){
	stop("SNP 1 MAF out of [0,1]")
}
if(opts$maf2 < 0 || opts$maf2 > 1){
	stop("SNP 2 MAF out of [0,1]")
}

if(length(eff1) != 3 || length(eff2) != 3){
	stop("Effect lengths must be exactly 3 numbers")
}

if(length(coeff) != 4){
	stop("Coefficients for model must be exactly 4 numbers")
}

if(CC && (opts$case <= 0 || opts$control <= 0)){
	stop("You must have at least one case and one control")
}

if(CC && any(coeff < 0)){
	stop("Cannot have negative coefficients in the model")
}

if(CC && sum(coeff) && is.null(opts$snr) > 1){
	stop("Coefficients must sum to <= 1")
}

if(min(eff1) != 0 || max(eff1) != 1){
	warning("Effect for SNP 1 not scaled from 0 to 1, scaling")
	eff1 = (eff1 - min(eff1)) / (max(eff1) - min(eff1))
}

if(min(eff2) != 0 || max(eff2) != 1){
	warning("Effect for SNP 2 not scaled from 0 to 1, scaling")
	#nprint(eff2)
	eff2 = (eff2 - min(eff2)) / (max(eff2) - min(eff2))
}

if(!is.null(opts$snr)){
	# Auto-set the pct. explained variance here
	S1 <- as.factor(1:3)
	S2 <- as.factor(1:3)
	d <- expand.grid(S1,S2)
	d <- data.frame(S1=d$Var1, S2=d$Var2)
	wt_vec <- as.vector(outer(dbinom(0:2,2,maf1), dbinom(0:2,2,maf2)))
	# Variance of bernoulli is p*(1-p), want weight to be 1/variance
	wt_vec <- wt_vec * (1-wt_vec)
	
	# now, set the phenotype matrix
	if(!pen){
		pentable <- coeff[1] + coeff[2]*outer(eff1, rep(1,length(eff2)) ) + coeff[3]*outer(rep(1,length(eff1)), eff2) + coeff[4]*outer(eff1, eff2)
		pen <- TRUE
	}
	
	# get the penetrance table to between 0 and 1
	pentable <- (pentable - min(pentable)) / (max(pentable) - min(pentable))
	
	d$Y <- pentable[3*(as.numeric(d$S2)-1) + as.numeric(d$S1)]
	
	mod <- lm(Y~S1+S2, data=d, weights=wt_vec)
	
	noise_var = 1
	smod = summary(mod)
	#print(smod)
	#print(smod$sigma)
	# in this case, there is no interaction
	if(smod$r.squared == 1){
		mod <- lm(Y~1, data=d, weights=wt_vec)
		smod <- summary(mod)
		#print(mod)
		#print(smod)
		#print(smod$sigma)
		rm <- weighted.mean(mod$residuals, mod$weights)
	}
	
	
	
	# scale the penetrance by the amount of unexplained variance times the signal-to-noise ratio (snr)

	# In case/control data, set the odds ratio based on the amount of signal
	# i.e. OR = exp(SNR/sigma) 
	if(CC){
		minP = 1/(1+exp(1/smod$sigma * opts$snr))
		penDiff = 1-2*minP
		pentable <- minP + pentable * penDiff
	}else{
		pentable <- pentable / smod$sigma * opts$snr
	}
	
	#print(pentable)
			
}

id_num=1

cat(c("#IID", "Pheno", "SNP1_1", "SNP1_2", "SNP2_1", "SNP2_2"), sep=" ")
cat("\n")

while (n_case + n_control > 0 && n_samp > 0){
	snp1 = sum(runif(2) < maf1)
	snp2 = sum(runif(2) < maf2)
	
	snpvec = as.numeric(c(snp1 > 1, snp1 > 0, snp2 > 1, snp2 > 0)) + 1
	
	if(pen){
		prob = pentable[snp1+1, snp2+1]
	}else{
		prob = coeff[1] + coeff[2]*eff1[snp1+1] + coeff[3]*eff2[snp2+1] + coeff[4]*eff1[snp1+1]*eff2[snp2+1]
	}
	
	print = F
	if(CC){
		status = runif(1) < prob
		if (status == 1 && n_case > 0){
			n_case = n_case - 1
			print=T
		} else if( status == 0 && n_control > 0){
			n_control = n_control - 1
			print = T
		}
	} else {
		status = prob + rnorm(1)
		n_samp = n_samp - 1
		print = T
	}
	
	
	if(print){
		cat(c(id_num,status, snpvec), sep=" ")
		id_num = id_num + 1
		cat("\n")		
	}
	
}

