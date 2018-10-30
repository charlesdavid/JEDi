## JEDi: Java Essential Dynamics Inspector
### JEDi 1.0 RELEASE NOTES: This is the new and updated version of JED 1.0

* JEDi is a powerful tool for examining the dynamics of proteins from trajectories derived from MD or Geometric simulations, or NMA.
* There are five types of PCA:
	* Residue / Hierarchical
	* Cartesian All-Atom
	* Cartesian Heavy-Atom
	* Cartesian Alpha-Carbon
	* Distance-Pair
* There are three statistical models: COV, CORR, and PCORR, for the Cartesian and Distance-Pair PCA.
* The idea behind the development of JEDi was to not only analyze a trajectory, but also compare it to other dynamical analyses.
* To do this, JEDi uses subspace analysis to assess the similarity of essential vector spaces defined by each dynamical method.
* Additionally, JEDi computes the Free Energy of a trajectory from two Order Parameters (weighted DVPs) so that a FE surface can be constructed.
* Finally, JEDi produces PDB frames and script files for viewing movies in PyMol(TM) of both individual PCA modes and a superposition of Essential Modes.


### Thank you for using JED and JEDi software:  

* Please report bugs by submitting a new issue tagged with BUG:  
	* Be sure to provide enough data to reproduce the error  
	* Include input files, errors, and Java stack traces 

* Please suggest new features by submitting a new issue tagged with NEW FEATURE or ENHANCEMENT

* Please contact me if you wish to become a collaborator on this software:  

Charles David: cdavid@carolina.rr.com  

##### Last Update: November, 2018
* Major Feature Enhancements:
	* Now offering 5 types of PCA, with 3 PCA Models (except Hierarchical which only uses COV Model).
		* The new Hierarchical approach allows for all-atom resolution analysis of very large proteins in reasonable computational time.
		* Cartesian PCA now does multiple levels of coarse-graining from all atom to alpha carbon only.
			* Each of the different types of analyses uses its own residue subsets, and its own number of modes.
	* Created new JEDI input file: Now a simple file of KEY-VALUE pairs.
	* Fixed PyMol(TM) script output (now works with educational license)
* New JavaDocs
* New Executable Jars
* New Tests

