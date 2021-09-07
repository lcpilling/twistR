# twistR
TWIST (Triangulation WIthin A STudy) analysis in R. 

The TWIST framework estimates the 'genetically moderated treatment effect' (GMTE). This is akin to the population attributable fraction (PAF) in pharmaco-epidemiology: the excess disease or outcomes in individuals prescribed a specific medication that are attributed to genetic variant(s). TWIST extends the PAF by incorporating information on both treated and untreated individuals. With this larger set of information we show that four analysis approaches for estimating the GMTE are possible. Each one relies on a different set of assumptions to work correctly and provides estimates that are largely uncorrelated with one another. A decision framework is used to decide when a particular estimation strategy is most appropriate and how specific estimators can be combined to further improve efficiency. Triangulation of evidence from different data sources, each with their inherent biases and limitations, is becoming a well established principle for strengthening causal analysis. Hence we called our framework 'Triangulation WIthin a STudy' (TWIST) in order to emphasise that an analysis in this spirit is also possible within a single data set, using causal estimates that are approximately uncorrelated, but reliant on different sets of assumptions.

TWIST therefore robustly estimates how many events (outcome can be mortality, disease, symptoms, etc) in the genotype group could be avoided if this group could experience the same treatment effect as non-carriers.

Preprint on medrxiv https://doi.org/10.1101/2021.05.04.21256612

## Installation
To install `twistR` directly from GitHub use the `devtools` package (to install: `install.packages("devtools")`)

`devtools::install_github("lukepilling/twistR")`

To update the package just run the above command again.

