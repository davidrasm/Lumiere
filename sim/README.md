Simulation code is provided to generate mock phylogenies to test the implementation of the the marginal fitness birth-death model in BEAST2.

The function main_fitBD_sim.m can be run to a single stochastic simulation and generate a phylogeny under a birth-death model with mutational fitness effects at multiple sites:
-This function expects a parameter file as an input argument. Such a parameter file can be genererated by modifying write_randomParams_file.m
-run_stochTreeSim_seqEvo.m is the actual workhorse function that runs the tree simulation.
-simmaptree_plot.m is provided to visualize the resulting simulated trees

The function batch_treeSim.m can be used to automate running multiple stochastic simulations and generate multiple trees at once.
-writeBeastXMLRandomParams.m will convert the simulation output into an XML file that can be run in BEAST.

The function test_margFitBirthDeath_likelihood.m can be used to test an independent implmentation of the marginal fitness birth-death model in Matlab
-The function compute_margFitBirthDeath_likelihood.m computes the likelihood of a phylogeny given a set of parameters, a phylo tree object and a fitness model. Examples of how to generate each of these are given in test_margFitBirthDeath_likelihood.m
-The function compute_birthDeathLikelihood.m implements the original multi-type birth-death model of Stadler & Bonhoeffer. It can be used to compare the MFBD approximation to the exact model tracking all possible genotypes in state space.

Phylo.m is a class for phylogenetic tree objects in Matlab. It contains a number of useful methods for working with phylogenetic trees.


