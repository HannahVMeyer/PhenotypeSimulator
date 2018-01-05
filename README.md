# PhenotypeSimulator
*PhenotypeSimulator* allows for the flexible simulation of phenotypes from different genetic and non-genetic (noise) components. 

In quantitative genetics, genotype to phenotype mapping is commonly realised by fitting a linear model to the genotype as the explanatory variable and the phenotype as the response variable. Other explanatory variable such as additional sample measures (e.g. age, height, weight) or batch effects can also be included. For linear mixed models, in addition to the fixed effects of the genotype and the covariates, different random effect components can be included, accounting for population structure in the study cohort or environmental effects. The application of linear and linear mixed models in quantitative genetics ranges from genetic studies in model organism such as yeast and *Arabidopsis thaliana* to human molecular, morphological or imaging derived traits. Developing new methods for increasing numbers of sample cohorts, phenotypic measurements or complexity of phenotypes to analyse, often requires the simulation of datasets with a specific underlying phenotype structure. 

*PhenotypeSimulator* allows for the simulation of complex phenotypes under different models, including genetic variant effects and infinitesimal genetic effects (reflecting population structure) as well as correlated, non-genetic covariates and observational noise effects. Different phenotypic effects can be combined into a final phenotype while controlling for the proportion of variance explained by each of the components. For each component, the number of variables, their distribution and the design of their effect across traits can be customised. 

The current CRAN version of *PhenotypeSimulator* is: 0.1.3

An update of version 0.1.3 has been submitted to CRAN (awaiting approval) and can already be downloaded from this github repository via
```{r}
install.packages("devtools")
devtools::install_github("HannahVMeyer/PhenotypeSimulator")
```

## Major changes from version 0.1.3 to the current github version 0.2.0:

**Input**
1. *PhenotypeSimulator* now includes readStandardGenotypes which can read externally simulated or user-provided genotypes in plink, genome, oxgen (hapgen/impute2), bimbam or simple delimited format.
1. A user-specified correlation matrix can be provided for the simulation of the correlatedBdEffects.
1. Short option flags for command-line use of *PhenotypeSimulator* were removed.

**Output**
1. *PhenotypeSimulator* provides the option to save the simulated phenotypes and genotypes in formats compatible with a number of commonly used genetic association software (gemma, bimbam, plink, snptest) via writeStandardOutput.
1. Intermediate phenotype components are now saved per default.
1. Saving additional subsets of the simulated data has been removed.

**Variance components**
1. Genotype simulation and kinship estimation: functions for genotype simulation and kinship estimation have been rewritten for
    significant speed-ups of the computation time [benchmarking](https://github.com/HannahVMeyer/PhenotypeSimulator-profiling).
    
1. geneticFixedEffects and noiseFixedEffects:
    1. The effect size distributions of the shared effects are now modelled as the product of two exponential distributions
        (to yield an approximately uniform distributions) or the product of a normal distribution with user-specified
        parameters and a standard normal distribution.
        
    1. The independent effects can now be specified to affect the same subset or different subsets of traits (via 
        keepSameIndependent).
        
    1. The overall number of traits affected by the effects can now be specified via pTraitsAffected.
    

1. correlatedBgEffects: the additional correlation between the traits can be specified by the user by providing an external
    correlation matrix. 




