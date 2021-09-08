# twistR
TWIST (Triangulation WIthin A STudy) analysis in R. 

The TWIST framework estimates the 'genetically moderated treatment effect' (GMTE). This is akin to the population attributable fraction (PAF) in pharmaco-epidemiology: the excess disease or outcomes in individuals prescribed a specific medication that are attributed to genetic variant(s). TWIST extends the PAF by incorporating information on both treated and untreated individuals. With this larger set of information we show that four analysis approaches for estimating the GMTE are possible. Each one relies on a different set of assumptions to work correctly and provides estimates that are largely uncorrelated with one another. A decision framework is used to decide when a particular estimation strategy is most appropriate and how specific estimators can be combined to further improve efficiency. Triangulation of evidence from different data sources, each with their inherent biases and limitations, is becoming a well established principle for strengthening causal analysis. Hence we called our framework 'Triangulation WIthin a STudy' (TWIST) in order to emphasise that an analysis in this spirit is also possible within a single data set, using causal estimates that are approximately uncorrelated, but reliant on different sets of assumptions.

TWIST therefore robustly estimates how many events (outcome can be mortality, disease, symptoms, etc) in the genotype group could be avoided if this group could experience the same treatment effect as non-carriers.

Preprint on medrxiv https://doi.org/10.1101/2021.05.04.21256612

## Table of Contents
  - [Installation](#installation)
  - [Model types performed](#model-types-performed)
  - [Function arguments](#function-arguments)
      - [gmte_binary()](#gmte_binary)
      - [gmte_continuous()](#gmte_continuous)
      - [gmte_aalen()](#gmte_aalen)
  - [Output and example](#output-and-example)


## Installation
To install `twistR` from GitHub use the `devtools` package:

`devtools::install_github("lukepilling/twistR")`

To update the package just run the above command again.

## Model types performed
`twistR` can perform GMTE analysis on three types of outcome:
1. `gmte_binary()` performs logistic regression analysis of a binary outcome (yes/no [1]/[0]). Requires the [`margins`](https://cran.r-project.org/web/packages/margins/) package.
2. `gmte_continuous()` performs linear regression analysis of a continuous outcome (quantitative measure, such as LDL cholesterol).
3. `gmte_aalen()` performs Aalen additive hazards model i.e. a time-to-event analysis (such as mortality). Requires the [`timereg`](https://cran.r-project.org/web/packages/timereg/) package.

## Function arguments

#### gmte_binary()

Argument | Description
-------- | -----------
Y | The binary outcome variable name (string) which appears in data.frame `D`.
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`.
D | A data.frame containing the above variables.
Link | Link function for the `glm()` - needs to be one of "logit","probit" or "identity". If unspecified the default is "logit".

#### gmte_continuous()

Argument | Description
-------- | -----------
Y | The continuous outcome variable name (string) which appears in data.frame `D`.
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`.
D | A data.frame containing the above variables.

#### gmte_aalen()

Argument | Description
-------- | -----------
Y_t0 | Variable name (string) for when participants "enter" the model, which appears in data.frame `D`. When participants enter the model (can be all 0s if Y_t1 is time since start of exposure). Variable can be in date format from the `as.Date()` function, or numeric.
Y_t1 | Variable name (string) for when participants "exit" the model, which appears in data.frame `D`. Either days since start of exposure period (numeric) or in date format from the `as.Date()` function.
Y_d | Variable name for the binary "event" variable (string) which appears in data.frame `D`. 
T | The treatment variable name (string) which appears in data.frame `D`. Assumed to be binary.
G | The genotype variable name (string) which appears in data.frame `D`. Normally binary (e.g. comparing homozygous rare individuals to the rest of the population).
Z | A string containing the model covariates to appear in the `glm()` models (for example "age+sex"). All need to be in data.frame `D`. Unless otherwise specified covariates will be assumed to be time invarying i.e. the `const()` wrapper will be added  See `aalen()` documentation in the [`timereg`](https://cran.r-project.org/web/packages/timereg/) package.
D | A data.frame containing the above variables.

## Output and example

For each GMTE function a object of class `twistR_GMTE` is returned. This contains the full model outputs from each individual analysis performed (GMTE1, GMTE0, MR, RGMTE, and CAT --  see paper) in addition to a summary of all the models performed and whether the models can be combined to improve estimation in `$FullCombined`. For example:

```R
# Example using: 
#   - a binary outcome (high LDL) - in data.frame D
#   - a binary treatment variable (on statins) - in data.frame D
#   - a binary genotype (SLCO1B1*5 homozygotes) - in data.frame D
#   - adjustment for age and genetic principal components of ancestry 1 to 10 - all in data.frame D
Y="ldl_high"
T="statin"
G="slco1b1_5_hmz"
Z="age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
Link="logit"
results=gmte_binary(Y,T,G,Z,D,Link)
```

Prints the following results table:

 Model | Est | SE | P.Est | Qstat | Qp | Combine?
------ | --- | -- | ----- | ----- | -- | --------
CAT | -1.315816e+01 | 0.07672673 | 0.000000e+00 | NA | NA | NA
GMTE1 | 5.823592e-02 | 0.01373771 | 2.243879e-05 | NA | NA | NA
GMTE0 | -8.747968e-04 | 0.00596922 | 8.834862e-01 | NA | NA | NA
RGMTE | 7.525366e-02 | 0.01630214 | 3.908640e-06 | NA | NA | NA
MR | -8.593626e-04 | 0.04927101 | 9.860844e-01 | NA | NA | NA
RGMTE_MR | 6.774350e-02 | 0.01547698 | 1.202974e-05 | 2.150891 | 0.1424872 | 1
RGMTE_CAT | -4.963459e-01 | 0.01594618 | 0.000000e+00 | 28462.588988 | 0.0000000 | 0
MR_CAT | -3.842414e+00 | 0.04145881 | 0.000000e+00 | 20820.492332 | 0.0000000 | 0
GMTE1_CAT | -3.522933e-01 | 0.01352266 | 0.000000e+00 | 28749.387485 | 0.0000000 | 0
RGMTE_MR_CAT | -4.493672e-01 | 0.01517140 | 0.000000e+00 | 28554.130703 | 0.0000000 | 0

* For a binary analysis, the `Est` is the difference in risk between the genotype groups.
* For a continuous analysis, the `Est` is the difference in outcome (units) between the genotype groups.
* For a time-to-event (Aalen) analysis, the `Est` is the difference in number of events per person year between the genotype groups.

To understand which estimate is best to use, the user needs to consider the following:
* the context of the specific treatment and population,
* the assumptions tested by each model,
* whether a combination of estimates (such as the MR and RGMTE estimates) is optimum.

It is not so simple as to just use the Robust GMTE (RGMTE) estimate, for example. For this reason we do not automatically give a recommendation when the functions are executed. For the combined estimates the `Combine?` column simply reports whether the p-value for the Q-statistic (`Qp`) is >0.05 i.e. the estimates do not significantly differ. This does not necesssarily mean it is the best choice. For futher information on the assumptions tested and decision framework please see the published manuscript (open access in PLOS Genetics). 
