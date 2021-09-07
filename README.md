# twistR
TWIST (Triangulation WIthin A STudy) analysis in R. 

The TWIST framework estimates the 'genetically moderated treatment effect' (GMTE). This is akin to the population attributable fraction (PAF) in pharmaco-epidemiology: the excess disease or outcomes in individuals prescribed a specific medication that are attributed to genetic variant(s). TWIST extends the PAF by incorporating information on both treated and untreated individuals. With this larger set of information we show that four analysis approaches for estimating the GMTE are possible. Each one relies on a different set of assumptions to work correctly and provides estimates that are largely uncorrelated with one another. A decision framework is used to decide when a particular estimation strategy is most appropriate and how specific estimators can be combined to further improve efficiency. Triangulation of evidence from different data sources, each with their inherent biases and limitations, is becoming a well established principle for strengthening causal analysis. Hence we called our framework 'Triangulation WIthin a STudy' (TWIST) in order to emphasise that an analysis in this spirit is also possible within a single data set, using causal estimates that are approximately uncorrelated, but reliant on different sets of assumptions.

TWIST therefore robustly estimates how many events (outcome can be mortality, disease, symptoms, etc) in the genotype group could be avoided if this group could experience the same treatment effect as non-carriers.

Preprint on medrxiv https://doi.org/10.1101/2021.05.04.21256612

## Installation
To install `twistR` from GitHub use the `devtools` package:

`devtools::install_github("lukepilling/twistR")`

To update the package just run the above command again.

## Models
`twistR` can perform GMTE analysis on three types of outcome:
1. `gmte_binary()` performs logistic regression analysis of a binary outcome (yes/no [1]/[0]). Requires `margins` package.
2. `gmte_continuous()` performs linear regression analysis of a continuous outcome (quantitative measure, such as LDL cholesterol).
3. `gmte_aalen()` performs Aalen additive hazards model i.e. a time-to-event analysis (such as mortality). Requires `timereg` package.

## Function arguments
`gmte_binary()`

Argument | Description
-------- | -----------
Y | The binary outcome variable name (string) which appears in data.frame `D`.
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`.
D | A data.frame containing the above variables.
Link | Link function for the `glm()` - needs to be one of "logit","probit" or "identity". If unspecified the default is "logit".

`gmte_continuous()`

Argument | Description
-------- | -----------
Y | The continuous outcome variable name (string) which appears in data.frame `D`.
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`.
D | A data.frame containing the above variables.

`gmte_aalen()`

Argument | Description
-------- | -----------
Y_t0 | Variable name (string) for when participants "enter" the model, which appears in data.frame `D`. When participants enter the model (can be all 0s if Y_t1 is time since start of exposure). Variable can be in date format from the `as.Date()` function, or numeric.
Y_t1 | Variable name (string) for when participants "exit" the model, which appears in data.frame `D`. Either days since start of exposure period (numeric) or in date format from the `as.Date()` function.
Y_d | Variable name for the binary "event" variable (string) which appears in data.frame `D`. 
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`. Unless otherwise specified covariates will be assumed to be time invarying i.e. the `const()` wrapper will be added  See `aalen()` documentation in `timereg` package.
D | A data.frame containing the above variables.

