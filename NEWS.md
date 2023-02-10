# twistR 0.1.4

* `FullCombined` object is now a data.frame with a column for "Model" rather than using rownames of a matrix
* New function `gmte_plot()` creates a forest plot of estimates when provided with a `twistR_GMTE` object

# twistR 0.1.3

* Fix bug with gmte_combine where function fails if CAT estimate doesn't converge - allows returning of GMTE (etc) even if failure of one estimate

# twistR 0.1.2

* Add 'verbose' option to help with error checking/debugging, or if the user simply wants lots of output.

# twistR 0.1.1

* Update calculation of Tshat to use linear regression model throughout. Allows use of multi-allelic genetic variants.

# twistR 0.1.0

* First release of package.
