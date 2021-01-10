package support;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;
import supportIO.Input_Parameters;

/**
 * Class PCA: Provides multiple methods to construct the COV(Q), CORR(R), and PCORR(P) matrices.
 * 
 * These methods include the use of shrinkage for generating optimal estimators, and numerical stability.
 * 
 * Methods to get R from Q, P from R-inv or Q-inv, and construct the N/3 Reduced matrices are also provided.
 * 
 * Provides a hook for the Jama eigenvalue decomposition.
 * 
 * Note: IN ALL CASES, the assumption is that ROWS are variables and COLS are observations (frames).
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class PCA
{
	final int ROWS, COLS;
	final static double FLOOR = Input_Parameters.FLOOR, Threshold_COV = Input_Parameters.THRESHOLD_COV;
	double lambda_cov, lambda_corr, percent_sparse;
	final Matrix data_matrix, centered_data_matrix, data_means, data_sigmas, data_variances;
	Matrix Q, R, K, Q_sample, R_sample, K_Sample, S, T_COV, T_CORR;
	final Descriptive_Stats ds;
	final Row_Center_Data rcd;

	/* ****************************** CONSTRUCTORS ********************************************************************************* */

	/**
	 * Constructor for the PCA support class
	 * 
	 * The input is a data matrix, normally the aligned coordinates
	 * 
	 * Note: The format of the data matrix is: ROWS are variables, and COLS are instances (frames)
	 * 
	 * @param data The data matrix (aligned coordinates)
	 */
	public PCA(Matrix data)
	{
		this.data_matrix = data; // The ROWS are the VARIABLES, the COLS are the FRAMES
		this.COLS = data.getColumnDimension(); // Number of FRAMES
		this.ROWS = data.getRowDimension(); // Number of VARIABLES

		this.rcd = new Row_Center_Data(data_matrix);
		this.data_means = rcd.get_variable_means(); // Matrix of row means
		this.data_sigmas = rcd.get_variable_standard_deviations(); // Matrix of row standard deviations
		this.data_variances = rcd.get_variable_variances(); // Matrix of row variances
		this.centered_data_matrix = rcd.get_row_centered_data(); // Mean-Centered data matrix
		this.ds = new Descriptive_Stats();

		this.Q = new Matrix(ROWS, ROWS); // Covariance Matrix
		this.K = new Matrix(COLS, COLS); // Kernel Matrix

		this.Q_sample = new Matrix(ROWS, ROWS); // Sample Covariance Matrix
		this.R_sample = new Matrix(ROWS, ROWS); // Sample Correlation Matrix

		this.T_COV = new Matrix(ROWS, ROWS); // Target for shrinking COV matrix: Diagonal, Un-Equal Variances.
		this.S = new Matrix(ROWS, ROWS); // Optimized Estimator of the Population Covariance.

		this.T_CORR = Matrix.identity(ROWS, ROWS); // Target for shrinking CORR matrix: Identity.
		this.R = new Matrix(ROWS, ROWS); // Optimized Estimator of the Population Correlation.
	}

	/* ****************************** METHODS ********************************************************************************* */

	/**
	 * Method to calculate a symmetric square matrix, i.e., a covariance matrix (A * A-transpose)/N --> ROW Space:
	 * 
	 * Note: Using MEAN CENTERED data
	 * 
	 * @return Returns the Covariance Matrix, Q
	 */
	public Matrix get_C_matrix()
	{
		Q = centered_data_matrix.times(centered_data_matrix.transpose());
		Q.timesEquals(Math.pow((COLS), -1));
		return Q;
	}

	/**
	 *
	 * Method to calculate a symmetric square matrix, i.e., a kernel matrix (A-transpose * A )/N --> COL Space:
	 * 
	 * Note: Using MEAN CENTERED data
	 * 
	 * @return Returns the Covariance Matrix, Q
	 */
	public Matrix get_K_matrix()
	{
		K = (centered_data_matrix.transpose()).times(centered_data_matrix);
		K.timesEquals(Math.pow((COLS), -1));
		return K;
	}

	/**
	 * Method to shrink a covariance matrix to the Target: "Diagonal Un-Equal Variances" Note: Using the non-optimized specified shrinkage intensity
	 * 
	 * @return Returns the Shrunk Covariance Matrix Estimator, S
	 */
	public Matrix Shrink_COV_Matrix(Matrix cov, double intensity)
	{
		int rows = cov.getRowDimension();
		Matrix Target = new Matrix(rows, rows); // This is the target for shrinking.
		Matrix S = new Matrix(rows, rows); // This is the new estimator of the population covariance.

		for (int i = 0; i < rows; i++)
			{
				double s_ii = cov.get(i, i);
				Target.set(i, i, s_ii);
			}
		S = (Target.times(intensity)).plus(cov.times(1 - intensity));
		return S;
	}

	/**
	 * Method to shrink a correlation or partial correlation matrix to the Target: "Identity matrix" Note: Using the non-optimized specified shrinkage intensity
	 * 
	 * @return Returns the Shrunk Correlation Matrix
	 */
	public Matrix Shrink_CORR_Matrix(Matrix corr, double intensity)
	{
		int rows = corr.getRowDimension();
		Matrix Target = Matrix.identity(rows, rows); // This is the target for shrinking.
		Matrix s = (Target.times(intensity)).plus(corr.times(1 - intensity));
		return s;
	}

	/**
	 * Method to shrink a kernel matrix to the Target: "Diagonal Un-Equal Variances" Note: Using the non-optimized specified shrinkage intensity
	 * 
	 * @return Returns the Shrunk Kernel Matrix Estimator, K
	 */
	public Matrix Shrink_K_Matrix(Matrix kernel, double intensity)
	{
		int rows = kernel.getRowDimension();
		Matrix Target = new Matrix(rows, rows); // This is the target for shrinking.
		Matrix S = new Matrix(rows, rows); // This is the new estimator of the population covariance.

		for (int i = 0; i < rows; i++)
			{
				double s_ii = kernel.get(i, i);
				Target.set(i, i, s_ii);
			}
		S = (Target.times(intensity)).plus(kernel.times(1 - intensity));
		return S;
	}

	/**
	 * Method to construct a covariance estimator by applying optimal shrinkage to the sample covariance matrix, using the target "Diagonal Unequal Variances"
	 * 
	 * Note: The optimal shrinkage intensity (lambda) is determined
	 * 
	 * Also computes the raw sample covariance and correlation matrices
	 * 
	 * @return Returns the Covariance Matrix Estimator, S
	 */
	public Matrix get_Covariance_Estimator_Using_Optimized_Shrinkage()
	{
		double sum_Sij_squared = 0;
		double sum_var_Sij = 0;

		for (int i = 0; i < ROWS; i++) // Construct the cov and corr matrices
			{
				double var = data_variances.get(i, 0);
				double sigmaX = data_sigmas.get(i, 0);
				Q_sample.set(i, i, var);
				R_sample.set(i, i, 1d);
				T_COV.set(i, i, var);
				Matrix var_X = centered_data_matrix.getMatrix(i, i, 0, COLS - 1);
				for (int j = i + 1; j < ROWS; j++) // Using property of symmetric matrix... Need i!=j for lambda calculation
					{
						Matrix var_Y = centered_data_matrix.getMatrix(j, j, 0, COLS - 1);
						double sigmaY = data_sigmas.get(j, 0);
						double[] var_XY = (var_X.arrayTimes(var_Y)).getRowPackedCopy();
						double mean_XY = ds.get_mean(var_XY);
						double cov = (mean_XY);
						double corr = cov / (sigmaX * sigmaY);
						sum_Sij_squared += (cov * cov);
						Q_sample.set(i, j, cov);
						Q_sample.set(j, i, cov);
						R_sample.set(i, j, corr);
						R_sample.set(j, i, corr);
					}
			}
		for (int k = 0; k < COLS; k++) // Get the variances of the individual elements of the covariance matrix
			{
				for (int i = 0; i < ROWS; i++)
					{
						double x = centered_data_matrix.get(i, k);
						for (int j = i + 1; j < ROWS; j++) // Need i!=j for lambda calculation
							{
								double y = centered_data_matrix.get(j, k);
								double Wkij = x * y;
								double Wij = Q_sample.get(i, j);
								double V = (Wkij - Wij) * (Wkij - Wij);
								sum_var_Sij += V;
							}
					}
			}
		double VAR_Sij = sum_var_Sij * COLS / Math.pow((COLS - 1), 3);
		lambda_cov = (VAR_Sij / sum_Sij_squared);  // Calculate the optimal shrinkage intensity
		lambda_cov = Math.max(0, Math.min(1, lambda_cov));
		S = (T_COV.times(lambda_cov)).plus(Q_sample.times(1 - lambda_cov)); // Calculating the optimized covariance estimator
		return S;
	}

	/**
	 * Method to construct a covariance estimator by applying optimal shrinkage to the sample covariance matrix, using the target "Diagonal Unequal Variances"
	 * 
	 * Note: The optimal shrinkage intensity (lambda) is determined
	 * 
	 * Also computes the raw sample covariance and correlation matrices
	 * 
	 * @return Returns the Covariance Matrix Estimator, S
	 */
	public Matrix get_Correlation_Estimator_Using_Optimized_Shrinkage()
	{
		double sum_Rij_squared = 0;
		double sum_var_Rij = 0;

		for (int i = 0; i < ROWS; i++)
			{
				double sigmaX = data_sigmas.get(i, 0);
				Matrix var_X = centered_data_matrix.getMatrix(i, i, 0, COLS - 1);
				for (int j = i + 1; j < ROWS; j++) // Using property of symmetric matrix...
					{
						Matrix var_Y = centered_data_matrix.getMatrix(j, j, 0, COLS - 1);
						double sigmaY = data_sigmas.get(j, 0);
						double[] var_XY = (var_X.arrayTimes(var_Y)).getRowPackedCopy();
						double mean_XY = ds.get_mean(var_XY);
						double cov = (mean_XY);
						double corr = cov / (sigmaX * sigmaY);
						sum_Rij_squared += (corr * corr);
					}
			}
		for (int k = 0; k < COLS; k++) // Get the variances of the individual elements of the correlation matrix
			{
				for (int i = 0; i < ROWS; i++)
					{
						double x = centered_data_matrix.get(i, k);
						for (int j = i + 1; j < ROWS; j++) // Need i!=j for lambda calculation
							{
								double y = centered_data_matrix.get(j, k);
								double Wkij = x * y;
								double Wij = R_sample.get(i, j);
								double V = (Wkij - Wij) * (Wkij - Wij);
								sum_var_Rij += V;
							}
					}
			}
		double VAR_Rij = sum_var_Rij * COLS / Math.pow((COLS - 1), 3);
		lambda_corr = (VAR_Rij / sum_Rij_squared);  // Calculate the optimal shrinkage intensity
		lambda_corr = Math.max(0, Math.min(1, lambda_corr));
		R = (T_CORR.times(lambda_corr)).plus(R_sample.times(1 - lambda_corr)); // Calculating the optimized covariance estimator
		return R;
	}

	/**
	 * Method to numerically stabilize a covariance matrix:
	 * 
	 * Values less than FLOOR are set to zero Values greater than FLOOR but less than 'Threshold_COV' are set to 'Threshold_COV'
	 * 
	 * @return Returns the modified Covariance Matrix
	 */
	public Matrix Threshold_COV_Matrix(Matrix cov_matrix)
	{
		int dim = cov_matrix.getColumnDimension();
		for (int i = 0; i < dim; i++)
			{
				for (int j = i; j < dim; j++)
					{
						double val = cov_matrix.get(i, j);
						double test = Math.abs(val);
						double sign = Math.signum(val);
						if (test < Threshold_COV)
							{
								if (test < FLOOR) val = 0d;
								else
									{
										val = (sign * Threshold_COV);
										cov_matrix.set(i, j, val);
										cov_matrix.set(j, i, val);
									}
							}
					}
			}
		return cov_matrix;
	}

	/**
	 * Method to calculate the sample covariance matrix (without matrix multiplication): 'Q_sample'
	 * 
	 * Note: Using Mean-Centered data.
	 * 
	 * The sample correlation 'R_sample' matrix is also constructed.
	 * 
	 * @return Returns the sample Covariance Matrix, 'Q_sample'
	 */
	public Matrix get_Covariance_Matrix_by_Definition()
	{
		for (int i = 0; i < ROWS; i++)
			{
				double var = data_variances.get(i, 0);
				double sigmaX = data_sigmas.get(i, 0);
				Q_sample.set(i, i, var);
				R_sample.set(i, i, 1d);
				Matrix var_X = centered_data_matrix.getMatrix(i, i, 0, COLS - 1);
				for (int j = i + 1; j < ROWS; j++)
					{
						Matrix var_Y = centered_data_matrix.getMatrix(j, j, 0, COLS - 1);
						double sigmaY = data_sigmas.get(j, 0);
						double[] var_XY = (var_X.arrayTimes(var_Y)).getRowPackedCopy();
						double mean_XY = ds.get_mean(var_XY);
						double cov = (mean_XY);
						double corr = cov / (sigmaX * sigmaY);
						Q_sample.set(i, j, cov);
						Q_sample.set(j, i, cov);
						R_sample.set(i, j, corr);
						R_sample.set(j, i, corr);
					}
			}
		return Q_sample;
	}

	/**
	 * Method to calculate a covariance matrix using QR Factorization of the centered data matrix
	 * 
	 * @return Returns a sample Covariance Matrix, Q
	 */
	public Matrix get_covariance_matrix_QRF()
	{
		QRDecomposition QR = new QRDecomposition((centered_data_matrix.transpose()));
		Matrix R = QR.getR();
		Matrix COV = R.transpose().times(R);
		return COV;
	}

	/**
	 * /** Method to calculate the reduced C-matrix (covariance or correlation).
	 * 
	 * Input is a 3Nx3N C-matrix.
	 * 
	 * Output is a NxN reduced C-matrix
	 * 
	 * @return Returns the Reduced Covariance Matrix, rCOV
	 */
	public Matrix get_Reduced_C_Matrix(Matrix Q)
	{
		int m = Q.getRowDimension() / 3;
		int nx = 0;
		int ny = m;
		int nz = 2 * m;

		Matrix rCOV = new Matrix(m, m);

		for (int i = 0; i < m; i++)
			{
				int a1 = nx + i;
				int a2 = ny + i;
				int a3 = nz + i;

				for (int j = 0; j < m; j++)
					{
						int b1 = nx + j;
						int b2 = ny + j;
						int b3 = nz + j;
						double cov = (Q.get(a1, b1) + Q.get(a2, b2) + Q.get(a3, b3));
						rCOV.set(i, j, cov);
					}
			}
		return rCOV;
	}

	/**
	 * Method to calculate a reduced Dynamical matrix (inverse C-matrix).
	 * 
	 * Input is a 3Nx3N Dynamical matrix, DYN
	 * 
	 * Output is a NxN Reduced Dynamical matrix, rDYN
	 * 
	 * @return Returns the Reduced Dynamical matrix, rDYN
	 */
	public Matrix get_Reduced_DYN_Matrix(Matrix DYN)
	{
		int m = DYN.getRowDimension() / 3;
		int nx = 0;
		int ny = m;
		int nz = 2 * m;

		Matrix rDYN = new Matrix(m, m);

		for (int i = 0; i < m; i++)
			{
				int a1 = nx + i;
				int a2 = ny + i;
				int a3 = nz + i;

				for (int j = 0; j < m; j++)
					{
						int b1 = nx + j;
						int b2 = ny + j;
						int b3 = nz + j;
						double dyn = -(DYN.get(a1, b1) + DYN.get(a2, b2) + DYN.get(a3, b3));
						rDYN.set(i, j, dyn);
					}
			}
		return rDYN;
	}

	/**
	 * General method to calculate a reduced C-matrix (covariance or correlation) with specified number of components, m.
	 * 
	 * Input is a mNxmN C-matrix, and number of components, m
	 * 
	 * Output is a NxN reduced C-matrix
	 * 
	 * @return Returns the Reduced Covariance Matrix, rCOV
	 */
	public Matrix get_Reduced_C_Matrix(Matrix Q, int components)
	{
		int i = 0, j = 0;
		int m = Q.getRowDimension() / components;
		double cov = 0;
		Matrix rCOV = new Matrix(m, m);

		for (int k = 0; k < components; k++)
			{
				int nk = k * components;

				for (i = 0; i < m; i++)
					{
						int a1 = nk + i;

						for (j = 0; j < m; j++)
							{
								int b1 = nk + j;

								cov = (Q.get(a1, b1));
							}
					}
				rCOV.set(i, j, cov);
			}
		return rCOV;
	}

	/**
	 * General method to calculate the reduced Dynamical matrix (inverse C-matrix) with specified number of components, m.
	 * 
	 * Input is a mNxmN Dynamical matrix, and number of components, m.
	 * 
	 * Output is a NxN reduced Dynamical matrix, rDYN
	 * 
	 * @return Returns the Reduced Dynamical Matrix, rDYN
	 */
	public Matrix get_Reduced_DYN_Matrix(Matrix DYN, int components)
	{
		int m = DYN.getRowDimension() / components;
		int nx = 0;
		int ny = m;
		int nz = 2 * m;

		Matrix rDYN = new Matrix(m, m);

		for (int i = 0; i < m; i++)
			{
				int a1 = nx + i;
				int a2 = ny + i;
				int a3 = nz + i;

				for (int j = 0; j < m; j++)
					{
						int b1 = nx + j;
						int b2 = ny + j;
						int b3 = nz + j;
						double dyn = -(DYN.get(a1, b1) + DYN.get(a2, b2) + DYN.get(a3, b3));
						rDYN.set(i, j, dyn);
					}
			}
		return rDYN;
	}

	/**
	 * 
	 * Method to calculate a covariance matrix: (A * A-transpose/N)
	 * 
	 * Note: The data are assumed to be mean centered
	 * 
	 * @return Returns the sample Covariance Matrix, Q
	 */
	public Matrix get_C_Matrix_NO_CENTERING(Matrix A)
	{
		Matrix Q = A.times(A.transpose());
		Q.timesEquals(Math.pow((A.getColumnDimension()), -1));
		return Q;
	}

	/**
	 * 
	 * Method to calculate a kernel matrix: (A-transpose * A/N)
	 * 
	 * Note: The data are assumed to be mean centered
	 * 
	 * @return Returns the Kernel Matrix, Q
	 */
	public Matrix get_K_Matrix_NO_CENTERING(Matrix A)
	{
		Matrix Q = (A.transpose()).times(A);
		Q.timesEquals(Math.pow((A.getColumnDimension()), -1));
		return Q;
	}

	/**
	 * Method to calculate the sample correlation matrix R from a sample covariance matrix Q
	 * 
	 * @return Returns the sample Correlation Matrix, R
	 */
	public Matrix get_Correlation_Matrix_from_Covariance_Matrix(Matrix cov)
	{
		int rows = cov.getRowDimension();
		Matrix R = new Matrix(rows, rows);
		for (int i = 0; i < rows; i++)
			{
				double sigmaX = Math.sqrt(cov.get(i, i));
				if (sigmaX == 0)
					{
						System.err.println("Warning: Diagonal element of covariance matrix is zero... Setting to Numerical Floor, default = 1.000E-16");
						sigmaX = FLOOR;
					}
				R.set(i, i, 1d);
				for (int j = i + 1; j < rows; j++)
					{
						double sigmaY = Math.sqrt(cov.get(j, j));
						if (sigmaY == 0)
							{
								System.err.println("Warning: Diagonal element of covariance matrix is zero... Setting to Numerical Floor, default = 1.000E-16");
								sigmaY = FLOOR;
							}
						double element = cov.get(i, j) / (sigmaX * sigmaY);
						R.set(j, i, element);
						R.set(i, j, element);
					}
			}
		return R;
	}

	/**
	 * Method to calculate the correlation matrix from the data.
	 * 
	 * Note: The data are MEAN CENTERED before the calculation.
	 * 
	 * @return Returns the Correlation Matrix, R
	 */
	public Matrix get_Correlation_Matrix()
	{
		get_Covariance_Matrix_by_Definition();
		return getR();
	}

	/**
	 * Method to calculate the Partial-Correlation Matrix from a Precision matrix (inverse covariance matrix) or Anti_Image matrix (inverse correlation matrix)
	 * 
	 * Input is inverse covariance or correlation matrix
	 * 
	 * @return Returns the Partial Correlation Matrix, PCORR
	 */
	public Matrix get_Partial_Correlation_Matrix(Matrix invC)
	{
		int rows = invC.getRowDimension();
		Matrix pcorr = new Matrix(rows, rows);
		for (int i = 0; i < rows; i++)
			{
				double x = invC.get(i, i);
				pcorr.set(i, i, 1d);
				for (int j = 0; j < i; j++)
					{
						double y = invC.get(j, j);
						double r = invC.get(i, j);
						double xy = (x * y);
						double p = (-r / (Math.sqrt(xy)));
						pcorr.set(j, i, p);
						pcorr.set(i, j, p);
					}
			}
		return pcorr;
	}

	/**
	 * Wrapper Method to calculate the eigenvalue decomposition of a matrix.
	 * 
	 * Note: The matrix argument should be a positive-definite matrix for real eigenvalues.
	 * 
	 * See the JAMA package for more details
	 * 
	 * @param PDM The matrix to factor: Should be positive-definite
	 * @return The eigenvalue decomposition holding the eigenvalues and eigenvectors.
	 */
	public EigenvalueDecomposition get_eigenvalue_decomposition(Matrix PDM)
	{
		EigenvalueDecomposition evd = new EigenvalueDecomposition(PDM);
		return evd;
	}

	public Matrix Sparsify_Matrix(Matrix NSM, double threshold)
	{
		int rows = NSM.getRowDimension();
		double total_count = rows * rows;

		int count = 0;
		Matrix sparse = new Matrix(rows, rows);
		for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < rows; j++)
					{
						double r = NSM.get(i, j);
						double val = Math.abs(r);
						if (val < threshold)
							{
								sparse.set(i, j, 0);
								count++;
							}
						if (val >= threshold) sparse.set(i, j, r);
					}
			}
		percent_sparse = (count / total_count) * 100.0;
		return sparse;
	}

	public Matrix Sparsify_COV_Matrix_Using_CORR(Matrix Q, Matrix R, double R_threshold)
	{
		int rows = Q.getRowDimension();
		int total_count = rows * rows;
		int count = 0;
		Matrix sparse = new Matrix(rows, rows);
		for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < rows; j++)
					{
						double q = Q.get(i, j);
						double r = R.get(i, j);
						double val = Math.abs(r);
						if (val < R_threshold)
							{
								sparse.set(i, j, 0);
								count++;
							}
						if (val >= R_threshold) sparse.set(i, j, q);
					}
			}
		percent_sparse = (count / total_count) * 100;
		return sparse;
	}

	// ************************** GETTERS ******************************************************

	/**
	 * @return Returns the optimal shrinkage intensity for COV matrix
	 */
	public double get_lambda_cov()
	{
		return lambda_cov;
	}

	/**
	 * @return Returns the optimal shrinkage intensity for CORR matrix
	 */
	public double get_lambda_corr()
	{
		return lambda_corr;
	}

	/**
	 * @return Returns the means of the centered variables
	 */
	public Matrix getData_means()
	{
		return data_means;
	}

	/**
	 * @return Returns the standard deviations of the centered variables
	 */
	public Matrix getData_sigmas()
	{
		return data_sigmas;
	}

	/**
	 * @return Returns the variances of the centered variables
	 */
	public Matrix getData_variances()
	{
		return data_variances;
	}

	/**
	 * @return Returns the mean centered variables
	 */
	public Matrix getCentered_data_Matrix()
	{
		return centered_data_matrix;
	}

	/**
	 * @return Returns the sample covariance matrix, Q_sample
	 */
	public Matrix getQ()
	{
		return Q_sample;
	}

	/**
	 * @return Returns the sample correlation matrix, R_sample
	 */
	public Matrix getR()
	{
		return R_sample;
	}

	/**
	 * @return Returns the Target matrix for shrinkage, T_COV
	 */
	public Matrix getT_COV()
	{
		return T_COV;
	}

	/**
	 * @return Returns the Target matrix for shrinkage, T_CORR
	 */
	public Matrix getT_CORR()
	{
		return T_CORR;
	}

	/**
	 * @return Returns the shrunk covariance matrix, S
	 */
	public Matrix getS()
	{
		return S;
	}

	public double get_Percent_Sparse()
	{
		return percent_sparse;
	}
}