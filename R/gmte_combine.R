#' gmte_combine
#'
#' Internal function to combine estimates from multiple models and provide summary statistics
#'
#' @param Ests The estimates.
#' @param SEs Standard Errors for the estimates.
#' @param alpha Significance threshold (default 0.05)
#' @return An object containing the combined estimate (plus SE, P, and Q statistic)
#'
#' @author Jack Bowden; Luke Pilling.
#' @references Bowden, J., et al., The Triangulation WIthin A STudy (TWIST) framework for causal inference within Pharmacogenetic research. medrxiv https://doi.org/10.1101/2021.05.04.21256612
#' @export

gmte_combine = function(Ests=Ests,SEs=SEs,alpha=0.05){
	w  = 1/SEs^2
	ws = w/sum(w)
	E  = sum(ws*Ests)
	SE = sqrt(sum((ws*SEs)^2))
	
	E2  = (Ests-E)^2
	Q   = sum(w*E2)
	df  = length(Ests)-1
	Qp  = 1-pchisq(Q,df)
	Dec = 0 ; if(Qp>=alpha){Dec=1}
	
	P.Est    = 2*(1-pnorm(abs(E/SE)))
	R        = c(E,SE,P.Est,Q,Qp,Dec)
	ResultsX = matrix(R,nrow=1)
	return(ResultsX)
}
