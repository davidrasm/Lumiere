# Lumiere

BEAST 2 package for phylodynamics with adaptive molecular evolution using the marginal fitness birth-death model.

## System Requirements

Lumiere is a BEAST2 package. The latest version of the package has been compiled against BEAST2 v2.5.2.

Lumiere requires:
* BEAST2 - available for download at the [BEAST2 site](https://www.beast2.org/). Required to be between versions v2.4.0 and v2.5.2, older versions can be found on the [GitHub page](https://github.com/CompEvol/beast2/releases).
* MASTER - available for download on the [GitHub page](https://github.com/tgvaughan/MASTER).
* MultiTypeTree - available for download on the [GitHub page](https://github.com/tgvaughan/MultiTypeTree).
* beast-classic - available for download on the [GitHub page](https://github.com/BEAST2-Dev/beast-classic).
* BEASTLabs - available for download on the [GitHub page](https://github.com/BEAST2-Dev/BEASTLabs).
* sampled-ancestors - download can be found on the [GitHub page](https://github.com/CompEvol/sampled-ancestors).
* Eclipse - software for package development and running, download can be found on [Eclipse website](https://www.eclipse.org/downloads/).


Users should check the relevant version.xml files for each package to ensure dependecies are fulfilled.

## Installation guide for Eclipse

To compile BEAST with the Lumiere source code, we recommend using the Eclipse IDE.

Create folder in any accessible directory, this will serve as your eclipse workspace. For our purposes we will assume this directory is 
```
\home\eclipse-workspace
```
Place unzipped folders for all of the BEAST2 related packages (BEAST2, MASTER, MultiTypeTree, beast-classic, BEASTLabs, and sampled-ancestors) here. 

Once you have installed Eclipse, open and choose \home\eclipse-workspace as working directory. Add a new java project for each BEAST2 related package via File > New > Java Project, subsequently de-selecting "Use default location", and browsing for the folder of interest.

Once all folders have been added, they should be present on the Project Explorer tab where we can explore their subdirectories. 

Subsequently, for each package in the table below we need to do the following:
* Go to package properties, add the listed packages (found in parentheses in the table below) to the project references and the Projects tab in the java build path
* Under the Libraries tab, add all the JAR files(found in the table below). The package they are found in preceeds them in parentheses, all will be found in their packages respective lib folder.

| Package              | Project Build Path / References| Libraries
| -------------        | -------------                  | -------------
| (1) beast-classic    | (2)-(7)                        | (1).mtj
| (2) beast2           | (1)                            | all jar files found in (2)/lib
| (3) BEAST LABS       | (2)                            | (3).junit
| (4) Lumiere          | (1)-(3),(5)-(7)                | (4).commons-lang3, (2).commons-math3, (1).mtj, (6).jblas
| (5) MASTER           | (2)                            | all jar files found in (5)/lib
| (6) MultiTypeTree    |                                | all jar files found in (6)/lib
| (7) SampledAncestors | (2)                            |

## Test Example

Test that you have set up Lumiere correctly by running any of the xml files found in 
```
Lumiere/examples
```

This can be done by navigating to Run > Run Configurations > Java Application and selecting the New Launch Configuration button at the top right of the pop-up window. Subsequently ensure that under the Main tab the Project references the Lumiere folder and that the Main class is beast.app.BeastMCMC. Under the Arguments tab, reference the selected xml (eg. examples/ebola_makona_multFitnessBD_9genos_exactP0s.xml). Finally, navigate to the Classpath tab and add all of the projects to the Bootstrap Entries and all of the jar files should automatically become listed under User Entries.

## Running simulation tests

To run your own simulations please follow the directions found in 
```
Lumiere/sim/README.md
```


## Running on clusters

Running MCMC algorithms on personal computers can be cumbersome. Here we have included additional steps for how one can run this algorithm using a cluster.

* File > Export > Runnable Jar File
* Browse > Choose where to save > Save As 
* SFTP shell file, xml, and .jar file to cluster
* Utilize cluster specific submission protocols (sbatch filename.sh)


## Authors

* **David Rasmussen** with **Marco Hamins-Puertolas** 

