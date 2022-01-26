#' gmte_binary
#'
#' Performs analyses for the Triangulation WIthin A STudy (TWIST) framework to calcuate the ‘genetically moderated treatment effect’ (GMTE) for a binary outcome. The \code{gmte_binary} function returns an object of class \code{twistR_GMTE}, containing effect estimates from the individual tests (such as RGMTE) and the results when combinations are performed (such as RGMTE+MR).
#'
#' @param Y The binary outcome variable name (string) which appears in data.frame \code{D}.
#' @param T The treatment variable name (string) which appears in data.frame \code{D}. Assumed to be binary.
#' @param G The genotype variable name (string) which appears in data.frame \code{D}. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population) but can be additive (0, 1, 2) or a score of multiple variants if desired.
#' @param Z A string containing the model covariates to appear in the \code{glm()} models (for example "age+sex"). All need to be in data.frame \code{D}.
#' @param D A data.frame containing the above variables.
#' @param Link Link function for the \code{glm()} - needs to be one of "logit","probit" or "identity". If unspecified the default is "logit".
#' @param alpha The p-value threshold for the chi-square test, estimating whether two estimates should be combined. Default is 0.05.
#' @param verbose Return lots of output - useful for error checking. Default is FALSE.
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
#' # Example using a binary outcome (high LDL), binary treatment (statins), and binary genotype (SLCO1B1*5 homozygotes) variables
#' Y="ldl_high"
#' T="statin"
#' G="slco1b1_5_hmz"
#' Z="age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
#' Link="logit"
#' results=gmte_binary(Y,T,G,Z,D,Link)

gmte_binary = function(Y,T,G,Z,D,Link="logit",alpha=0.05,verbose=FALSE)
{
	start_time = Sys.time()
	cat("TWIST (Triangulation WIthin A STudy) analysis in R - binary outcome\n")
	require(margins)

	## check inputs
	if (class(Y) != "character")  stop("Outcome Y needs to be a variable name i.e. a string (class `character`)")
	if (class(T) != "character")  stop("Treatment T needs to be a variable name i.e. a string (class `character`)")
	if (class(G) != "character")  stop("Genotype G needs to be a variable name i.e. a string (class `character`)")
	if (class(Z) != "character")  stop("Covariates Z needs to be a formula i.e. a string (class `character`) of variable(s) in data.frame D e.g. \"age+sex\"")
	if (class(Link) != "character")  stop("Link needs to be a string (class `character`) and needs to be one of c(\"logit\",\"probit\",\"identity\")")

	if (! Link %in% c("logit","probit","identity"))  stop(paste0("Link [", Link, "] needs to be one of c(\"logit\",\"probit\",\"identity\")"))
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
	D[,"Tshat"] = glm(as.formula(paste0("Tstar~G+",Z)),data=D)$fitted

	## check outcome is a binary variable [0]/[1]
	Y_bin=unique(D[,"Y"])
	if (length(Y_bin) != 2) stop(paste0("Outcome Y [", Y, "] needs to be a binary variable with only two values: 0 or 1"))
	if (! Y_bin[1] %in% c(0,1) | ! Y_bin[2] %in% c(0,1)) stop(paste0("Outcome Y [", Y, "] needs to be a binary variable with only two values: 0 or 1"))

	## enough data for people on treatment with genotype?
	n_total=nrow(D)
	cat(paste0("- N with complete data [", n_total, "]\n"))

	n_treated=length(D[ D[,"T"] == 1 , "T"])
	cat(paste0("- N on treatment (T=1) [", n_treated, "]\n"))
	
	if (n_treated == 0) stop("Need to include some treated individuals")	
	if (n_treated < 100) cat("Warning: Low numbers of treated individuals - model may not converge\n")	

	n_untreated=n_total-n_treated
	if (n_untreated == 0) stop("Need to include untreated individuals for control group")	
	if (n_untreated < 100) cat("Warning: Low numbers of untreated individuals - model may not converge\n")	
	
	n_treated_geno=length(D[ D[,"T"] == 1 & D[,"G"] != 0 , "T"])
	n_treated_geno_outcome=length(D[ D[,"T"] == 1 & D[,"G"] != 0  & D[,"Y"] == 1 , "T"])
	
	if (n_treated_geno == 0) stop("Not enough observations for analysis")	
	if (n_treated_geno < 100) cat("Warning: Low numbers of Treated individuals carrying the Genotype - model may not converge\n")	
	
	if (n_treated_geno_outcome == 0) stop("Not enough participants for analysis")	
	if (n_treated_geno_outcome < 100) cat("Warning: Low numbers of Treated individuals carrying the Genotype experience Outcome - model may not converge\n")	

	cat("\n")

	####################
	## begin analysis ##
	####################

	# Corrected-As treated (CAT)
	cat("Run CAT model\n")
	D[,"Tcat"] = D[,"T"]*mean(D[,"G"][D[,"T"]==1])
	CATfit     = glm(as.formula(paste0("Y~Tcat+",Z)),family=binomial(link=Link),data=D)
	MarCAT     = summary(margins(CATfit))
	CAT        = MarCAT[MarCAT[,"factor"]=="Tcat",2]
	sCAT       = MarCAT[MarCAT[,"factor"]=="Tcat",3]
	pCAT       = MarCAT[MarCAT[,"factor"]=="Tcat",5]
	if (verbose)  print(CATfit)

	# GMTE(0)
	cat("Run GMTE(0) model\n")
	D[,"tt"] = (1-D[,"T"])
	D[,"ts"] = D[,"tt"]*D[,"G"]
	GMTE0fit = glm(as.formula(paste0("Y~tt+ts+",Z)),family=binomial(link=Link),data=D)
	MarGMTE0 = summary(margins(GMTE0fit))
	GMTE0    = MarGMTE0[MarGMTE0[,"factor"]=="ts",2]
	sGMTE0   = MarGMTE0[MarGMTE0[,"factor"]=="ts",3]
	pGMTE0   = MarGMTE0[MarGMTE0[,"factor"]=="ts",5]
	if (verbose)  print(GMTE0fit)

	# GMTE(1)
	cat("Run GMTE(1) model\n")
	GMTE1fit = glm(as.formula(paste0("Y~T+Tstar+",Z)),family=binomial(link=Link),data=D)
	MarGMTE1 = summary(margins(GMTE1fit))
	GMTE1    = MarGMTE1[MarGMTE1[,"factor"]=="Tstar",2]
	sGMTE1   = MarGMTE1[MarGMTE1[,"factor"]=="Tstar",3]
	pGMTE1   = MarGMTE1[MarGMTE1[,"factor"]=="Tstar",5]
	if (verbose)  print(MarGMTE1)

	# RGMTE
	cat("Run RGMTE model\n")
	RGMTEfit = glm(as.formula(paste0("Y~T+Tstar+Tshat+",Z)),family=binomial(link=Link),data=D) 
	MarRGMTE = summary(margins(RGMTEfit))
	RGMTE    = MarRGMTE[MarRGMTE[,"factor"]=="Tstar",2]
	sRGMTE   = MarRGMTE[MarRGMTE[,"factor"]=="Tstar",3]
	pRGMTE   = MarRGMTE[MarRGMTE[,"factor"]=="Tstar",5]
	if (verbose)  print(MarRGMTE)

	# MR
	cat("Run MR model\n")
	MRfit       = glm(as.formula(paste0("Y~Tshat+",Z)),family=binomial(link=Link),data=D)
	MarMR       = summary(margins(MRfit))
	MR          = MarMR[MarMR[,"factor"]=="Tshat",2]
	sMR         = MarMR[MarMR[,"factor"]=="Tshat",3]
	pMR         = MarMR[MarMR[,"factor"]=="Tshat",5]
	if (verbose)  print(MRfit)

	# Partial results DF
	if (verbose)
	{
		FullCombined      = matrix(nrow=5,ncol=3)
		FullCombined[1,]  = c(CAT,sCAT,pCAT)
		FullCombined[2,]  = c(GMTE0,sGMTE0,pGMTE0)
		FullCombined[3,]  = c(GMTE1,sGMTE1,pGMTE1)
		FullCombined[4,]  = c(RGMTE,sRGMTE,pRGMTE)
		FullCombined[5,]  = c(MR,sMR,pMR)
		colnames(FullCombined) = c("Est","SE","EstP")
		rownames(FullCombined) = c("CAT","GMTE0","GMTE1","RGMTE","MR")
		cat("\nResults (initial):\n")
		print(FullCombined)
	}

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

	output_list=list(model="gmte_binary",CAT=MarCAT,GMTE0=MarGMTE0,GMTE1=MarGMTE1,RGMTE=MarRGMTE,MR=MarMR,FullCombined=FullCombined)
	class(output_list)="twistR_GMTE"
	return(output_list)

}
