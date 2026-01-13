# GS_f  
This repository provides R scripts for evaluating the goodness-of-fit of the proposed model, the *f*-divergence-based Gaussian Symmetry (GS[*f*]) model, for multi-way contingency tables with ordinal categories.
The methodology builds on the Gaussian Symmetry framework and introduces an *f*-divergence-based extension that allows flexible model assessment.  

This repository accompanies the following paper:  
> Okahara, H. and Tahata, K. (2026). *Modeling asymmetry in multi-way contingency tables with ordinal categories via f-divergence*. *Statistical Papers*.

## CONTENTS  
1. `libraries.R`    : Loads the required R packages for the project.
2. `database.R`     : Provides an example dataset used in the manuscript.
3. `f_divergence.R` : Defines the various *f*-functions and derivatives.
4. `functions.R`    : Contains utility functions for model fitting and evaluation.
5. `constraints.R`  : Specifies linear constraints for the ME, VE, CE and ME_2 models.
6. `mph.Rcode.R`    : Placeholder for Professor Lang's `mph.fit` function (not publicly available). 
7. `main.R`         : Main script to fit the GS[*f*] model and nested models, conduct simulation analysis, and compute plug-in estimates of the potential parameters.



## SOFTWARE AVAILABILITY  
The `mph.fit` function and its supporting code are **not included** in this repository.  
According to the note on [Professor Lang’s website](https://homepage.divms.uiowa.edu/~jblang/mph.fitting/mph.fit.documentation.htm):  
> "To receive a free copy of the R programs referred to in this document, contact me (Joseph B. Lang) at joseph-lang@uiowa.edu.  
> [*Note: The statistical package R is freeware that is available at:* [http://cran.r-project.org](http://cran.r-project.org); *its syntax is very similar to that of the commercial software S-plus.*]"

If you intend to use the 'GS_f' implemented here, you must request the `mph.Rcode.R` directly from Professor Lang using the contact information above.
