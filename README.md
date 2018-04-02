# asexual-lineage-adaptation
Matlab code for simulations and inference of DFE (distribution of fitness effects) described in "The fates of mutant lineages and the distribution of fitness effects of beneficial mutations in laboratory budding yeast populations" EM Frenkel, BH Good, MM Desai. Genetics, 2014 (http://www.genetics.org/content/196/4/1217.long)

To get started, add the directories 'custom Matlab functs' and 'simulation code' to the path and open the file 'script_showMeSomeTrajs.m' found in the directory 'scripts to test code'.  The comments there and in the files simTraj.m and returnGenMutFunct.m should explain how to run the simulation and plot its results. Note that the tester scripts are broken up into sections, whose boundaries are defined by the double %% (see https://www.mathworks.com/help/matlab/matlab_prog/run-sections-of-programs.html) and only individual sections are meant to be run at a time. 

The DFE inference was done by performing many simulations of the data in parallel on a compute cluster and iteratively analyzing the results to determine which parameter values require further simulation. The code for this aspect of the work is not in a very portable state, so is not included in the repository. Likewise the code for the adaptation rate modeling (Fig. 6) needs to be cleaned up for distribution. I'd be happy to do so if anyone would like to use them, so please feel free to contact me.

-Evgeni Frenkel, genya@wi.mit.edu
