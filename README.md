eagle: Environment-ASE through Generalized LinEar modeling
====

EAGLE is an R package which uses Allele Specific Expression (ASE) quantified via RNA-seq to detect gene by environment (GxE) interactions. The underlying idea is that observing an association between the level of allelic imbalance at an exonic SNP, and an environmental factor, is evidence of a *cis*-regulatory element modulating the transcriptomic response to the environmental perturbation. The underlying statistical model is a binomial [Generalized Linear Mixed Model](https://en.wikipedia.org/wiki/Generalized_linear_mixed_model), with a random effect term used to account for count overdisperion. Model fitting is achieved using [Non-conjugate Variational Message Passing](https://papers.nips.cc/paper/4407-non-conjugate-variational-message-passing-for-multinomial-and-binary-regression). 

The paper describing EAGLE has now been published:

**Allele-specific expression reveals interactions between genetic variation and environment.**
*Knowles, D. A; Davis, J. R; Edgington, H.; Raj, A.; Favé, M.; Zhu, X.; Potash, J. B; Weissman, M. M; Shi, J.; Levinson, D.; Awadalla, P.; Mostafavi, S.; Montgomery, S. B; and Battle, A.*
[Nature Methods](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4298.html) 2017.

An early [bioRxiv preprint](http://biorxiv.org/content/early/2015/09/13/025874) is available. 

The code is on [github](https://github.com/davidaknowles/eagle). 

## System requirements

eagle has been tested under the following architectures: 
* Mac OS X Yosemite 10.10.2, R 3.1.2
* Ubuntu 14.04 LTS, R 3.2.2
* Red Hat 4.4.7, R 3.1.2  

I have not tried compiling under Windows. 

## Installation

You will need the RcppEigen R package installed. In R run

`install.packages("RcppEigen")`

To download the code

`git clone git@github.com:davidaknowles/eagle.git`

To compile+install, in this directory run

`R CMD INSTALL --build .`

Alternatively run 
```
# require(devtools)
devtools::install_github("davidaknowles/eagle")
```

## Usage

For a simple example script on synthetic data look at `test.R`. For a slightly more involved/realistic example, including eQTLs and realistic filtering options look at `big_test.R`. 

