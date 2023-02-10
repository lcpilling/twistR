#' Forest plot individual and combined TWIST/GMTE estimates
#'
#' Plot the individual and combined estimates from a Triangulation WIthin A STudy (TWIST) analysis. The \code{gmte_plot} function takes an object of class \code{twistR_GMTE}, containing effect estimates from the individual tests (such as RGMTE) and the results when combinations are performed (such as RGMTE+MR), and creates a forest plot, highlighting the individual and combined estimates, and indicating with a "*" when the combined estimate may be valid.
#'
#' @param x An object of class \code{twistR_GMTE} e.g., the output from \code{gmte_continuous}
#' @param plot_title A string to print as the plot title
#' @param plot_cat Logical. Plot the CAT estimates? (Default=TRUE)
#' @param cols Three colours to indiciate the three model types (GMTE0, individual estimates, combined estimates)
#' @param pchs Three point types to indiciate the three model types (GMTE0, individual estimates, combined estimates)
#'
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
#'
#' gmte_plot(results, plot_title = "SLCO1B1*5 effect on LDL during statin treatment")
#'
#' # If desired, remove CAT estimates for "cleaner" plot, as these are often orders of magnititude larger than the other estimates
#' gmte_plot(results, plot_title = "SLCO1B1*5 effect on LDL during statin treatment", plot_cat=FALSE)


gmte_plot = function(x, 
                     plot_title = "", 
                     plot_cat = TRUE,
                     cols = c("#f46036","#2e294e","#1b998b"),
                     pchs = c(15,16,23))  {

	## check input
	if (class(x)!="twistR_GMTE")  stop("`x` needs to be a twistR_GMTE object")

	## make sure objects only visible in plotting space
	par(xpd=FALSE)

	## get bits we need for twistR_GMTE object
	model = x$model
	res   = x$FullCombined

	## move CAT below other estimates
	rownames(res)=c(6,1:5,7:10)
	res=res[order(as.numeric(rownames(res))),]

	## add colours & point types
	res[,"col"] = c(cols[1],cols[2],cols[2],cols[2],cols[3],cols[2],cols[3],cols[3],cols[3],cols[3])
	res[,"pch"] = c(pchs[1],pchs[2],pchs[2],pchs[2],pchs[3],pchs[2],pchs[3],pchs[3],pchs[3],pchs[3])

	## if plotting CAT but it is 10* larger than other estimates give a warning.
	if (plot_cat & abs(res[res[,"Model"]=="CAT","Est"]) > 10*max(abs(res[res[,"Model"] %in% c("GMTE0","GMTE1","RGMTE","MR","RGMTE_MR"),"Est"])))  warning("The CAT estimate is substantially larger than other estimates. Consider re-running with the flag 'plot_cat=FALSE'")

	## remove CAT estimate if specified by user
	if (!plot_cat) res = res[ ! grepl("CAT",res[,"Model"]) , ]

	## add column for "combine" to use a star "*" rather than 1 in plot
	res[,"combine"] = ""
	res[res[,"Combine?"] %in% "1","combine"] = "*"

	## reverse order so that GMTE0 is at top of plot
	rownames(res)=1:nrow(res)
	res=res[rev(rownames(res)),]

	## calculate 95% CIs 
	res[,"Est_lci"] = res[,"Est"] - 1.96*(res[,"SE"])
	res[,"Est_uci"] = res[,"Est"] + 1.96*(res[,"SE"])

	## determine axis limits
	x_min = min(c(res[,"Est_lci"],res[,"Est_uci"]))
	x_max = max(c(res[,"Est_lci"],res[,"Est_uci"]))
	x_dist = x_max - x_min
	x_min = x_min - (x_dist*0.05)
	x_max = x_max + (x_dist*0.05)

	y_min = 0.8
	y_max = nrow(res)+0.5

	## determine x axis label
	x_lab = ""
	if (model == "gmte_continuous")  x_lab = "Mean difference in outcome"
	if (model == "gmte_binary")      x_lab = "Risk difference in outcome"
	if (model == "gmte_aalen")       x_lab = "Hazard difference per unit time"

	## create empty plot
	##   extend xlim to give space for labels
	##   extend ylim to give space for labels
	plot(0 , type = "n" , xlim=c(x_min-(x_dist/1.5),x_max) , ylim=c(y_min,y_max), cex=1, xaxt='n', yaxt='n', bty='n', xlab=x_lab, ylab="", main=plot_title)

	## add vertical reference line
	segments(x0=0, y0=0-1, x1=0, y1=y_max, col="#BBBBBB", lty=2, lwd=1.5)

	## add x axis
	axis(1, at=c(round(x_min,2),0,round(x_max,2)))

	## allow objects to run outside of plot
	par(xpd=TRUE)

	## add y-axis
	segments(x0=x_min-(x_dist/20), y0=1, x1=x_min-(x_dist/20), y1=nrow(res))
	for (ii in 1:nrow(res))  segments(x0=x_min-(x_dist/10), y0=ii, x1=x_min-(x_dist/20), y1=ii)
	text(x=x_min-(x_dist/8), y=1:nrow(res), labels=res[,"Model"], adj=0, cex=1, font=4, pos=2)

	## use "segments" to draw the CI lines
	segments(res[,"Est_lci"], 1:nrow(res), res[,"Est_uci"], 1:nrow(res),
		   lwd=1, lty=1, col=res[,"col"])

	## plot central estimate
	points(res[,"Est"], 1:nrow(res), 
		 pch=res[,"pch"], col=res[,"col"], bg=res[,"col"], cex=1.5)

	## indicate if combine should be used
	text(x_max, 1:nrow(res), pos=4, res[,"combine"], cex=1)

	## legend
	legend(x = x_min-(x_dist/1.2), 
		 y = y_min-0.2, 
		 c("GMTE0 (untreated)","Individual model","Combined model"),
		 bty = "n",
		 col=cols,
		 pt.bg=cols,
		 pch=pchs,
		 lty=1,
		 cex=0.8)

	## restore margins settings
	par(xpd=FALSE)

}

