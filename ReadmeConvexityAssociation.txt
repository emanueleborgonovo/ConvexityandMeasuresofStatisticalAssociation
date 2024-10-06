The subroutine TestsBivariatePromitEDW2W22Sobol.m is used to generate Figure 4.1
The subroutine calls:
sobolpoints.m to generate points with Sobol' random generator
bwsi.m to calculate the Square Wasserstein-2 measure of statistical association in Equation 57 of the paper.
mmd2si.m to calculate the energy distance sensitivity measure in Equation 58 using the appropriate kernel 

The subroutine MultivariateGaussianNumericEnergyDist.m is used to generate Figure 4.2.
This subroutine calls:
sobolpoints.m to generate points with Sobol' random generator
mmd2si.m to calculate the energy distance sensitivity measure using the appropriate kernel
edistmim.m which is an alternative to estimate the same quantity but without relying on a kernel
plotpp2.m that generates the probability-probability plot.


The file Multivariate_Wass2ED.mcdx is the Mathcad PCT Prime 8 file used by the authors to calculate the analytical expressions of the sensitivity measure $\xi ^{ED}(Y,X_{i})$ in Equation 59, for the case study. We have also printed the file in the corresponding .pdf.
