## JEDi: Java Essential Dynamics Inspector
### JEDi 1.0 RELEASE NOTES: This is the new and updated version of JED 1.0

* JEDi is a powerful tool for examining the dynamics of proteins from trajectories derived from MD or Geometric simulations, or NMA.
* There are 4 Levels of Coarse-Graining:
	* All-Atom
	* Heavy-Atom
	* BackkBone
	* Alpha-Carbon
* There are 4 Types of PCA:
	* Hierarchical
	* Cartesian
	* Displacement
	* Distance-Pair
* There are 3 PCA Models:
	* Covariance
	* Correlation
	* Partial-Correlation
* Kernel PCA post-processing of the reduced data is available, using an assortment of kernels including:
     * Mutual Information
     * Neural-Net
     * Gaussian
     * Euclidean
     * Polynomial Degrees 2, 3, 4, Difference of Cubes, Difference of Squares
     * Sine-Cosine
     * Linear
* NEW: JEDi creates PNG images of the scatterplots of the top two PCs or DVPs for every analysis, making the review of results simple.
* NEW: JEDi uses subspace analysis to assess the similarity of essential vector spaces defined by each dynamical method:
     * Inter-Model within each type of PCA
     * Inter-Type when subsets yield subspaces with matching dimensions
* Additionally, JEDi computes the Free Energy of a trajectory from two PCs (normed projections) so that a FE surface can be constructed.
* Finally, JEDi produces PDB frames and script files for viewing the frames in PyMol(TM) as movies:
     * Individual modes
     * The superposition of the Top Essential Modes.


### Thank you for using JED and JEDi software:

###### https://github.com/charlesdavid/JED
###### https://github.com/charlesdavid/JEDi (new source code will be available upon submission of our paper)

* Please report bugs by submitting a new issue tagged with BUG:  
     * Be sure to provide enough data to reproduce the error  
     * Include JED/JEDi input files, JED/JEDi errors, and associated Java stack traces  

* Please suggest new features by submitting a new issue tagged with NEW FEATURE or ENHANCEMENT

* Please contact me if you wish to become a collaborator on this software:  

Charles David: ccdavid@protonmail.com  

##### Last Update: December, 2018
* New Features
* New JavaDocs
* New Jars
* Major BUG fixes
* Added KPCA, PLOT_XY, and FES to new JEDi Driver
* Created new input file: Now a list of KEY=VALUE pairs
* Fixed the visualization class to handle the 4 levels of coarse graining
