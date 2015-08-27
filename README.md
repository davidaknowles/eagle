eagle
====

Environment-ASE through Generalized LinEar modeling

Note you will need the RcppEigen R package installed. In R run

install.packages("RcppEigen")

To compile+install, in this directory run

R CMD INSTALL --build .

For a simple example script on synthetic data look at 
test.R

eagle has been tested under the following architectures: 
Mac OS X Yosemite 10.10.2, R 3.1.2
Ubuntu 14.04 LTS, R 3.2.2
Red Hat 4.4.7, R 3.1.2  

I have not tried compiling under Windows. 