#' gmte_continuous
#'
#' Performs analyses for the Triangulation WIthin A STudy (TWIST) framework to calcuate the ‘genetically moderated treatment effect’ (GMTE) for a continuous outcome. The \code{gmte_continuous} function returns an object of class \code{twistR_GMTE}, containing effect estimates from the individual tests (such as RGMTE) and the results when combinations are performed (such as RGMTE+MR).
#'
#' @param Y The continuous outcome variable name (string) which appears in data.frame \code{D}.
#' @param T The treatment variable name (string) which appears in data.frame \code{D}. Assumed to be binary.
#' @param G The genotype variable name (string) which appears in data.frame \code{D}. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population) but can be additive (0, 1, 2) or a score of multiple variants if desired.
#' @param Z A string containing the model covariates to appear in the \code{glm()} models (for example "age+sex"). All need to be in data.frame \code{D}.
#' @param D A data.frame containing the above variables.
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
#' # Example using a continuous outcome (LDL), binary treatment (statins), and binary genotype (SLCO1B1*5 homozygotes) variables
#' Y="ldl"
#' T="statin"
#' G="slco1b1_5_hmz"
#' Z="age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
#' results=gmte_continuous(Y,T,G,Z,D)

gmte_continuous = function(Y,T,G,Z,D,alpha=0.05)
{
	start_time = Sys.time()
	cat("TWIST (Triangulation WIthin A STudy) analysis in R - continuous outcome\n")

	## check inputs
	if (class(Y) != "character")  stop("Outcome Y needs to be a variable name i.e. a string (class `character`)")
	if (class(T) != "character")  stop("Treatment T needs to be a variable name i.e. a string (class `character`)")
	if (class(G) != "character")  stop("Genotype G needs to be a variable name i.e. a string (class `character`)")
	if (class(Z) != "character")  stop("Covariates Z needs to be a formula i.e. a string (class `character`) of variable(s) in data.frame D e.g. \"age+sex\"")

	if (class(D) != "data.frame")  stop("D needs to be a data.frame")
	if (! Y %in% colnames(D))  stop(paste0("Outcome Y [", Y, "] needs to be in data.frame D"))
	if (! T %in% colnames(D))  stop(paste0("Treatment T [", T, "] needs to be in data.frame D"))
	if (! G %in% colnames(D))  stop(paste0("Genotype G [", G, "] needs to be in data.frame D"))
	Zs=strsplit(Z,"[+]|[*]")[[1]]
	for (Zx in Zs)  if (! Zx %in% colnames(D))  stop(paste0("Covariate [", Zx, "] needs to be in data.frame D"))

	cat(paste0("- Outcome Y [", Y, "]\n"))
	cat(paste0("- Treatment T [", T, "]\n"))
	cat(paste0("- Genotype G [", G, "]\n"))
	cat(paste0("- Covariates [", Z, "]\n"))

	## subset dataset to columns specified and remove NAs
	D=D[,colnames(D) %in% c(Y,T,G,Zs)]
	D=as.data.frame(na.omit(D))

	## create variables named Y, T and G for formulas, and compute T* (interaction between T and G)
	D[,"Y"]=D[,Y]
	D[,"T"]=D[,T]
	D[,"G"]=D[,G]
	D[,"Tstar"] = D[,"T"]*D[,"G"]

	## enough data for people on treatment with genotype?
	cat(paste0("- N with complete data [", nrow(D), "]\n"))

	n_treated=length(D[ D[,"T"] == 1 , "T"])
	cat(paste0("- N on treatment (T=1) [", n_treated, "]\n"))
	
	n_treated_geno=length(D[ D[,"T"] == 1 & D[,"G"] != 0 , "T"])
	
	if (n_treated_geno == 0) stop("Not enough observations for analysis")	
	if (n_treated_geno < 100) cat("Warning: Low numbers of Treated individuals carrying the Genotype - model may not converge\n")	

	cat("\n")

	####################
	## begin analysis ##
	####################

	# MR
	cat("Run MR model\n")
	D[,"Tshat"] = lm(as.formula(paste0("Tstar~G+",Z)),data=D)$fitted
	MRfit       = lm(as.formula(paste0("Y~Tshat+",Z)),data=D)
	MR          = summary(MRfit)$coef["Tshat",1]
	sMR         = summary(MRfit)$coef["Tshat",2]
	pMR         = summary(MRfit)$coef["Tshat",4]

	# Corrected-As treated (CAT)
	cat("Run CAT model\n")
	D[,"Tcat"] = D[,"T"]*mean(D[,"G"][D[,"T"]==1])
	CATfit     = lm(as.formula(paste0("Y~Tcat+",Z)),data=D)
	CAT        = summary(CATfit)$coef["Tcat",1]
	sCAT       = summary(CATfit)$coef["Tcat",2]
	pCAT       = summary(CATfit)$coef["Tcat",4]
	
	# GMTE(1)
	cat("Run GMTE(1) model\n")
	GMTE1fit = lm(as.formula(paste0("Y~T+Tstar+",Z)),data=D)
	GMTE1    = summary(GMTE1fit)$coef["Tstar",1]
	sGMTE1   = summary(GMTE1fit)$coef["Tstar",2]
	pGMTE1   = summary(GMTE1fit)$coef["Tstar",4]

	# GMTE(0)
	cat("Run GMTE(0) model\n")
	D[,"tt"] = (1-D[,"T"])
	D[,"ts"] = D[,"tt"]*D[,"G"]
	GMTE0fit = lm(as.formula(paste0("Y~tt+ts+",Z)),data=D)
	GMTE0    = summary(GMTE0fit)$coef["ts",1]
	sGMTE0   = summary(GMTE0fit)$coef["ts",2]
	pGMTE0   = summary(GMTE0fit)$coef["ts",4]

	# RGMTE
	cat("Run RGMTE model\n")
	RGMTEfit = lm(as.formula(paste0("Y~T+Tstar+Tshat+",Z)),data=D) 
	RGMTE    = summary(RGMTEfit)$coef["Tstar",1]
	sRGMTE   = summary(RGMTEfit)$coef["Tstar",2]
	pRGMTE   = summary(RGMTEfit)$coef["Tstar",4]

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

	cat("\nResults:\n")
	print(FullCombined)

	end_time = Sys.time()
	time_taken = round(as.numeric(end_time)-as.numeric(start_time),1)
	cat(paste0("\nAnalysis completed in ", time_taken, " seconds\n"))

	output_list=list(model="gmte_continuous",CAT=CATfit,GMTE0=GMTE0fit,GMTE1=GMTE1fit,RGMTE=RGMTEfit,MR=MRfit,FullCombined=FullCombined)
	class(output_list)="twistR_GMTE"
	return(output_list)

}
