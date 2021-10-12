#' gmte_aalen
#'
#' Performs analyses for the Triangulation WIthin A STudy (TWIST) framework to calcuate the ‘genetically moderated treatment effect’ (GMTE) for a time-to-event Aalen additive hazards model. The \code{gmte_aalen} function returns an object of class \code{twistR_GMTE}, containing effect estimates from the individual tests (such as RGMTE) and the results when combinations are performed (such as RGMTE+MR).
#'
#' @param Y_t0 Variable name (string) for when participants "enter" the model, which appears in data.frame \code{D}. When participants enter the model (can be all 0s if \code{Y_t1} is time since start of exposure). Variable can be in date format from the \code{as.Date()} function, or numeric.
#' @param Y_t1 Variable name (string) for when participants "exit" the model, which appears in data.frame \code{D}. Either days since start of exposure period (numeric) or in date format from the \code{as.Date()} function.
#' @param Y_d Variable name for the binary "event" variable (string) which appears in data.frame \code{D}.
#' @param T The treatment variable name (string) which appears in data.frame \code{D}. Assumed to be binary.
#' @param G The genotype variable name (string) which appears in data.frame \code{D}. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
#' @param Z A string containing the model covariates to appear in the \code{glm()} models (for example "age+sex"). All need to be in data.frame \code{D}.
#' @param D A data.frame containing the above variables.
#' @param Nsim Number of simulations to perform in the \code{aalen()} models. Default is 100.
#' @param alpha The p-value threshold for the chi-square test, estimating whether two estimates should be combined. Default is 0.05.
#' @return An object of class \code{twistR_GMTE} containing the following components:\describe{
#' \item{\code{CAT}}{The summary statistics from the Corrected As Treated (CAT) analysis.}
#' \item{\code{GMTE1}}{The summary statistics from the GMTE(1) analysis (i.e. in the treated individuals).}
#' \item{\code{GMTE0}}{The summary statistics from the GMTE(0) analysis (i.e. in the untreated individuals).}
#' \item{\code{MR}}{The summary statistics from the MR analysis.}
#' \item{\code{RGMTE}}{The summary statistics from the Robust GMTE analysis (GMTE1 corrected for GMTE0).}
#' \item{\code{FullCombined}}{The combined summary statistics from all analyses performed, inclduing combinations.}
#'}
#' @author Jack Bowden; Luke Pilling.
#' @references Bowden, J., et al., The Triangulation WIthin A STudy (TWIST) framework for causal inference within Pharmacogenetic research. PLoS Genetics. https://doi.org/10.1371/journal.pgen.1009783
#' @export
#' @examples
#' # Example using a time-to-event outcome (mortality), binary treatment (statins), and binary genotype (SLCO1B1*5 homozygotes) variables
#' Y_t0="date_first_statin"
#' Y_t1="date_of_death_or_censor"
#' Y_d="dead"
#' T="statin"
#' G="slco1b1_5_hmz"
#' Z="age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
#' results=gmte_aalen(Y_t0,Y_t1,Y_d,T,G,Z,D)

gmte_aalen = function(Y_t0,Y_t1,Y_d,T,G,Z,D,Nsim=100,alpha=0.05)
{
	start_time = Sys.time()
	cat("TWIST (Triangulation WIthin A STudy) analysis in R - Aalen additive hazards (time-to-event) model\n")
	require(timereg)

	## check inputs
	if (class(Y_t0) != "character")  stop("Outcome Y_t0 needs to be a variable name i.e. a string (class `character`)")
	if (class(Y_t1) != "character")  stop("Outcome Y_t1 needs to be a variable name i.e. a string (class `character`)")
	if (class(Y_d) != "character")  stop("Outcome Y_d needs to be a variable name i.e. a string (class `character`)")
	if (class(T) != "character")  stop("Treatment T needs to be a variable name i.e. a string (class `character`)")
	if (class(G) != "character")  stop("Genotype G needs to be a variable name i.e. a string (class `character`)")
	if (class(Z) != "character")  stop("Covariates Z needs to be a formula i.e. a string (class `character`) of variable(s) in data.frame D e.g. \"age+sex\"")

	if (class(D) != "data.frame")  stop("D needs to be a data.frame")
	if (! Y_t0 %in% colnames(D))  stop(paste0("Outcome Y_t0 [", Y_t0, "] needs to be in data.frame D"))
	if (! Y_t1 %in% colnames(D))  stop(paste0("Outcome Y_t1 [", Y_t1, "] needs to be in data.frame D"))
	if (! Y_d %in% colnames(D))  stop(paste0("Outcome Y_d [", Y_d, "] needs to be in data.frame D"))
	if (! T %in% colnames(D))  stop(paste0("Treatment T [", T, "] needs to be in data.frame D"))
	if (! G %in% colnames(D))  stop(paste0("Genotype G [", G, "] needs to be in data.frame D"))
	
	## check covariates are in data - make new formula with const() wrapper if required
	Zs=strsplit(Z,"[+]|[*]")[[1]]
	Zwrapped=""
	Zunwrapped=""
	for (ii in 1:length(Zs))  
	{
		Zx=Zs[ii]
		if (ii>1) Zwrapped=paste0(Zwrapped,"+")
		if (ii>1) Zunwrapped=paste0(Zunwrapped,"+")
		## if wrapper is supplied for variable in formula
		if (grepl("(",Zx,fixed=TRUE))
		{
			## strip this to check data is in D
			Zx=strsplit(strsplit(Zx,"(",fixed=TRUE)[[1]][2],")",fixed=TRUE)[[1]][1]
			## add to new formula
			Zwrapped=paste0(Zwrapped,Zs[ii])
			Zs[ii]=Zx
			Zunwrapped=paste0(Zunwrapped,Zx)
		} else {
			## no wrapper provided. Assume constant (time invarying)
			Zwrapped=paste0(Zwrapped,"const(",Zx,")")
			Zunwrapped=paste0(Zunwrapped,Zx)
		}
		if (! Zx %in% colnames(D))  stop(paste0("Covariate [", Zx, "] needs to be in data.frame D"))
	}

	cat(paste0("- Outcome Y_t0 [", Y_t0, "] i.e. when participants enter model\n"))
	cat(paste0("- Outcome Y_t1 [", Y_t1, "] i.e. when participants exit model\n"))
	cat(paste0("- Outcome Y_t0 [", Y_d, "] i.e. binary variable indicating event\n"))
	cat(paste0("- Treatment T [", T, "]\n"))
	cat(paste0("- Genotype G [", G, "]\n"))
	cat(paste0("- Covariates [", Zwrapped, "]\n"))

	## check variable formats
	if (class(D[,Y_t0]) != "numeric" & class(D[,Y_t0]) != "Date") stop(paste0("Outcome Y_t0 [", Y_t0, "] needs to be in `numeric` or `Date` format"))
	if (class(D[,Y_t1]) != "numeric" & class(D[,Y_t1]) != "Date") stop(paste0("Outcome Y_t1 [", Y_t1, "] needs to be in `numeric` or `Date` format"))
	if (class(D[,Y_d]) != "numeric" & class(D[,Y_d]) != "integer") stop(paste0("Outcome Y_d [", Y_d, "] needs to be `numeric`"))
	
	## if dates provided, convert to numeric 
	if (class(D[,Y_t0]) == "Date")  D[,Y_t0] = as.numeric(D[,Y_t0])
	if (class(D[,Y_t1]) == "Date")  D[,Y_t1] = as.numeric(D[,Y_t1])
	
	## subset dataset to columns specified and remove NAs
	D=D[,colnames(D) %in% c(Y_t0,Y_t1,Y_d,T,G,Zs)]
	D=as.data.frame(na.omit(D))

	## create variables named Y, T and G for formulas, and compute T* (interaction between T and G)
	D[,"Y_t0"]=D[,Y_t0]
	D[,"Y_t1"]=D[,Y_t1]
	D[,"Y_d"]=D[,Y_d]
	D[,"T"]=D[,T]
	D[,"G"]=D[,G]
	D[,"Tstar"] = D[,"T"]*D[,"G"]

	## enough data for people on treatment with genotype?
	cat(paste0("- N with complete data [", nrow(D), "]\n"))

	n_treated=length(D[ D[,"T"] == 1 , "T"])
	cat(paste0("- N on treatment (T=1) [", n_treated, "]\n"))
	
	n_treated_geno=length(D[ D[,"T"] == 1 & D[,"G"] != 0 , "T"])
	cat(paste0("- N on treatment (T=1) & carrying genotype (G!=0) [", n_treated_geno, "]\n"))
	
	n_treated_geno_outcome=length(D[ D[,"T"] == 1 & D[,"G"] != 0  & D[,"Y_d"] == 1 , "T"])
	cat(paste0("- N on treatment (T=1) & carrying genotype (G!=0) & experience event (Y_d=1) [", n_treated_geno_outcome, "]\n"))
	
	if (n_treated_geno == 0) stop("Not enough observations for analysis")	
	if (n_treated_geno < 100) cat("Warning: Low numbers of Treated individuals carrying the Genotype - model may not converge\n")	

	if (n_treated_geno_outcome == 0) stop("Not enough participants for analysis")	
	if (n_treated_geno_outcome < 100) cat("Warning: Low numbers of Treated individuals carrying the Genotype experience Outcome - model may not converge\n")	

	## any exit dates before enter dates?
	t1_before_t0=D[ D[,"Y_t1"] < D[,"Y_t0"],"Y_t1"]
	if (length(t1_before_t0)>0) warning("Some exit dates are before enter dates - check your input data")
	
	cat("\n")

	####################
	## begin analysis ##
	####################

	## survival object
	survival_object = Surv(D[,"Y_t0"],D[,"Y_t1"],D[,"Y_d"])

	# MR
	cat("Run MR model\n")
	D[,"Tshat"] = glm(as.formula(paste0("Tstar~G+",Zunwrapped)),data=D)$fitted
	MRfit       = invisible(aalen(as.formula(paste0("survival_object~const(Tshat)+",Zwrapped)),data=D,n.sim=Nsim))
	MR          = coef(MRfit)["const(Tshat)",1]
	sMR         = coef(MRfit)["const(Tshat)",2]
	pMR         = coef(MRfit)["const(Tshat)",5]

	# Corrected-As treated (CAT)
	cat("Run CAT model\n")
	D[,"Tcat"] = D[,"T"]*mean(D[,"G"][D[,"T"]==1])
	CATfit     = invisible(aalen(as.formula(paste0("survival_object~const(Tcat)+",Zwrapped)),data=D,n.sim=Nsim))
	CAT        = coef(CATfit)["const(Tcat)",1]
	sCAT       = coef(CATfit)["const(Tcat)",2]
	pCAT       = coef(CATfit)["const(Tcat)",5]

	# GMTE(1)
	cat("Run GMTE(1) model\n")
	GMTE1fit = invisible(aalen(as.formula(paste0("survival_object~const(T)+const(Tstar)+",Zwrapped)),data=D,n.sim=Nsim))
	GMTE1    = coef(GMTE1fit)["const(Tstar)",1]
	sGMTE1   = coef(GMTE1fit)["const(Tstar)",2]
	pGMTE1   = coef(GMTE1fit)["const(Tstar)",5]

	# GMTE(0)
	cat("Run GMTE(0) model\n")
	D[,"tt"] = (1-D[,"T"])
	D[,"ts"] = D[,"tt"]*D[,"G"]
	GMTE0fit = invisible(aalen(as.formula(paste0("survival_object~const(tt)+const(ts)+",Zwrapped)),data=D,n.sim=Nsim))
	GMTE0    = coef(GMTE0fit)["const(ts)",1]
	sGMTE0   = coef(GMTE0fit)["const(ts)",2]
	pGMTE0   = coef(GMTE0fit)["const(ts)",5]

	# RGMTE
	cat("Run RGMTE model\n")
	RGMTEfit = invisible(aalen(as.formula(paste0("survival_object~const(T)+const(Tstar)+const(Tshat)+",Zwrapped)),data=D,n.sim=Nsim))
	RGMTE    = coef(RGMTEfit)["const(Tstar)",1]
	sRGMTE   = coef(RGMTEfit)["const(Tstar)",2]
	pRGMTE   = coef(RGMTEfit)["const(Tstar)",5]

	# Combined methods
	cat("Combined methods\n")
	Ests         = c(MR,RGMTE); SEs = c(sMR,sRGMTE)
	RGMTE_MR     = gmte_combine(Ests,SEs,alpha)
	Ests         = c(CAT,RGMTE); SEs = c(sCAT,sRGMTE)
	RGMTE_CAT    = gmte_combine(Ests,SEs,alpha)
	Ests         = c(CAT,MR); SEs = c(sCAT,sMR)
	MR_CAT       = gmte_combine(Ests,SEs,alpha)
	Ests         = c(CAT,GMTE1); SEs = c(sCAT,sGMTE1)
	GMTE1_CAT    = gmte_combine(Ests,SEs,alpha)
	Ests         = c(MR,RGMTE,CAT); SEs = c(sMR,sRGMTE,sCAT)
	RGMTE_MR_CAT = gmte_combine(Ests,SEs,alpha)

	# Final output
	FullCombined      = matrix(nrow=10,ncol=6)
	FullCombined[1,]  = c(CAT,sCAT,pCAT,NA,NA,NA)
	FullCombined[2,]  = c(GMTE0,sGMTE0,pGMTE0,NA,NA,NA)
	FullCombined[3,]  = c(GMTE1,sGMTE1,pGMTE1,NA,NA,NA)
	FullCombined[4,]  = c(RGMTE,sRGMTE,pRGMTE,NA,NA,NA)
	FullCombined[5,]  = c(MR,sMR,pMR,NA,NA,NA)
	FullCombined[6,]  = RGMTE_MR
	FullCombined[7,]  = RGMTE_CAT
	FullCombined[8,]  = MR_CAT
	FullCombined[9,]  = GMTE1_CAT
	FullCombined[10,] = RGMTE_MR_CAT

	colnames(FullCombined) = c("Est","SE","EstP","Qstat","Qp","Combine?")
	rownames(FullCombined) = c("CAT","GMTE0","GMTE1","RGMTE","MR",
	                           "RGMTE_MR","RGMTE_CAT","MR_CAT","GMTE1_CAT","RGMTE_MR_CAT")

	cat("Results:\n")
	print(FullCombined)

	end_time = Sys.time()
	time_taken = round(as.numeric(end_time-start_time),1)
	cat("Analysis completed in ", time_taken, " seconds\n")

	output_list=list(model="gmte_aalen",CAT=CATfit,GMTE0=GMTE0fit,GMTE1=GMTE1fit,RGMTE=RGMTEfit,MR=MRfit,FullCombined=FullCombined)
	class(output_list)="twistR_GMTE"
	return(output_list)

}
