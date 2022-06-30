# Master-1---Economics---Thesis
#
# Supplementary material for the thesis: "Modeling the relation between adaptation investment and public debt"
#
# These files are here to reproduce the results of my thesis. There are divided in two files. 

**********************************************************************************************************
→ One folder contains the mod.file to run on Dynare. 

Millner_no_staggering.mod implements my model with no staggering of expenditures in the time-to-build. Most of the graphs in the thesis stem from this model. 
It is also used to implement progressively all items in the model, and it includes the possibility to run the model with other specifications for the utility function.

Millner_with_staggering.mod implements my model with the staggering of public expenditures as in Leeper(2010). However, as mentioned in the paper, this model is not complete enough to produce realistic impulse-response-functions.

Millner_optimal_vs_alternative.mod provides the model with different adaptation capital reaction function, whether it is the optimal reaction function or two alternatives scenarios.


**********************************************************************************************************
→ The second folder contains the m.file to use in Matlab. They are here to reproduce the graphs shown in the thesis.

The code assemble.m is used with tightfig to produce a figure displaying all graphs produced by Dynare after running a given model. It sets the figure size to have it corresponds to a A4 paper and it uses the tex names produced by the mod files. 

The code groupplot_assemble.m does the same as assemble.m, but to combine different outputs from different Dynare outputs. These Dynare outputs need to have the same variables names and in the same order. 

For tightfig.m, I retrieved the code from :
Richard Crozier (2022). tightfig(hfig) (https://www.mathworks.com/matlabcentral/fileexchange/34055-tightfig-hfig), MATLAB Central File Exchange. Retrieved June 30, 2022.
I only changed it very marginally to include the names of the variables and the x-axis. 

Comparisons_Millner_ADDICE.mlx is the live script used to explore the properties of the different adaptation function. 
