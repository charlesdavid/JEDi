## JEDi: Java Essential Dynamics Inspector
### JEDi RELEASE NOTES: JEDi is the updated version of JED
##### Beta Testing Completed
##### Now in release 1.0

* JEDi is a powerful tool for inspecting the dynamics of proteins from trajectories derived from MD or Geometric simulations, or NMA.
* All JEDi features are controlled via an input file containing a list of KEY=VALUE pairs
* There are 5 Levels of Coarse-Graining:
	* All-Atom by residue
	* Heavy-Atom by residue
	* Back-Bone by residue
	* Alpha-Carbon by residue
	* Atom List
* There are 4 Types of PCA:
	* Hierarchical Cartesian
	* Direct Cartesian
	* Distance-Pair
	* Residue Pair Interaction
* There are 3 PCA Models:
	* Covariance
	* Correlation
	* Partial-Correlation
* There is the option to perform Outlier Processing:
	* Outliers Removal is perfomed at a specified cutoff
	* Outliers Selection is perfomed at a specified cutoff
	* The results of outlier removal and selection are compared
* There is the option to perform Sparsification of the Correlation and Partial Correlation Matrices:
	* Separate cutoffs are provided for each matrix
	* When both Correlation and Partial Correlation models are selected, an adjacency matrix is provided to help identify activator and suppressor variables
	* The results of sparse and non-sparse analyses are compared
* All JEDi output is provided as high quality PNG images for immediate inspection:
	* PCs and DVPs
	* Square Modes
	* Scree Plots
	* Variable statistics
	* MSA scores
	* Eigenvector Collectivity
	* Iterated RMSIPs with comparison to random and with Z-scores for significance
* Kernel PCA post-processing of the reduced data is available, using an assortment of kernels including:
     * Mutual Information
     * Neural-Net (Tanh)
     * Gaussian
     * Log
     * Circular
     * Mahalanobis
     * Euclidean
     * Polynomial Degrees 2, 3, 4, and cross XY
     * Linear
* JEDi uses subspace analysis to assess the similarity of essential vector spaces defined by each dynamical method:
     * Inter-Model within each type of PCA
     * Inter-Type when subsets yield subspaces with matching dimensions
* JEDi computes the Free Energy of a trajectory from two PCs (normed projections) so that a FE surface can be constructed.
* JEDi produces PDB frames that capture the PCA modes and script files for viewing those frames in PyMol(TM) as movies:
     * Individual modes
     * The superposition of the Top Essential Modes
* For each JEDi run, a comprehensive log file is generated with time stamp recoding the parameters and key results of each analysis

### Thank you for using JED and JEDi software:

###### https://github.com/charlesdavid/JED
###### https://github.com/charlesdavid/JEDi

* Please report bugs by submitting a new issue tagged with BUG:  
     * Be sure to provide enough data to reproduce the error  
     * Include JED/JEDi input files, JED/JEDi errors, and associated Java stack traces  

* Please suggest new features by submitting a new issue tagged with NEW FEATURE or ENHANCEMENT

* Please contact me if you wish to become a collaborator on this software:

* Check out this link for a demonstration video for using JEDi:

  * https://www.youtube.com/watch?v=3I08oGQShwg

Charles David: charles.david@plantandfood.co.nz 

##### Last Software Update: December, 2020
* New Java Docs
* New executable JAR files
