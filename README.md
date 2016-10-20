eagle
====

Environment-ASE through Generalized LinEar modeling

The paper is on bioRxiv:

Allele-specific expression reveals interactions between genetic variation and environment

David A Knowles, Joe R Davis, Anil Raj, Xiaowei Zhu, James B Potash, Myrna M Weissman, Jianxin Shi, Douglas F Levinson, Sara Mostafavi, Stephen B Montgomery, Alexis Battle

http://biorxiv.org/content/early/2015/09/13/025874

Note you will need the RcppEigen R package installed. In R run

`install.packages("RcppEigen")`

To compile+install, in this directory run

`R CMD INSTALL --build .`

For a simple example script on synthetic data look at 
test.R. For a slightly more involved/realistic example, including eQTLs and realistic filtering options look at big_test.R. _

eagle has been tested under the following architectures: 
* Mac OS X Yosemite 10.10.2, R 3.1.2
* Ubuntu 14.04 LTS, R 3.2.2
* Red Hat 4.4.7, R 3.1.2  

I have not tried compiling under Windows. 
