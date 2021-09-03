
GMTE_combine = function(Ests=Ests,SEs=SEs,alpha=0.05){
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
