package jedi;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * JED class PCA: Constructs the COV(Q), CORR(R), and PCORR(P) matrices. Provides a hook for the Jama eigenvalue decomposition. Note: The assumption is that ROWS are variables and
 * COLS are observations.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 */

public class PCA
	{
		int ROWS, COLS;
		final double ST_Threshold = 1.00E-9;
		Matrix data_matrix, centered_data_matrix, data_means, data_sigmas;

		/* ****************************** CONSTRUCTOR ********************************************************************************* */

		/**
		 * Constructor to perform PCA on data: Note: The data matrix takes ROWS as variables and COLS as instances.
		 * 
		 * @param data The data matrix A
		 */
		public PCA(Matrix data)
			{
				data_matrix = data;
				COLS = data.getColumnDimension();
				ROWS = data.getRowDimension();

			}

		/* ****************************** METHODS ********************************************************************************* */

		/**
		 * Method to calculate a covariance matrix: AA' Note: The data are ROW CENTERED before the calculation. No Normalization is done.
		 * 
		 * @return Returns the Covariance Matrix, Q
		 */
		Matrix get_covariance_matrix()
			{
				System.out.println("\n\tCentering Data: ");
				Row_Center_Data rcd = new Row_Center_Data(data_matrix);
				centered_data_matrix = rcd.get_row_centered_data();
				data_means = rcd.get_variable_means();
				data_sigmas = rcd.get_variable_sigmas();
				System.out.println("\tGetting COV Matrix: ");
				Matrix Q = centered_data_matrix.times(centered_data_matrix.transpose());
				// Q.timesEquals(Math.pow((COLS - 1), -1));
				return Q;
			}

		/**
		 * Method to threshold-stabilize a covariance matrix by setting to 0 any entry whose absolute value is less than ST_Threshold = 1.00E-9
		 * 
		 * @return Returns the Covariance Matrix, Q
		 */
		Matrix stabilize_Covariance_Matrix(Matrix Covariance_Matrix)
			{
				System.out.println("\tStabilizing Covariance Matrix: ");

				for (int i = 0; i < Covariance_Matrix.getColumnDimension(); i++)
					{
						for (int j = 0; j < Covariance_Matrix.getColumnDimension(); j++)
							{
								double test = Covariance_Matrix.get(i, j);
								if (Math.abs(test) < ST_Threshold)
									{
										test = 0;
										Covariance_Matrix.set(i, j, test);
									}
							}
					}

				return Covariance_Matrix;
			}

		/**
		 * Method to calculate the covariance matrix (without matrix multiplication): Note: The data are ROW CENTERED before the calculation.
		 * 
		 * @return Returns the Covariance Matrix, Q
		 */
		Matrix get_covariance_matrix_elegant()
			{
				System.out.println("\tGetting the Covariance Matrix without matrix multiplication: ");
				Row_Center_Data rcd = new Row_Center_Data(data_matrix);
				data_means = rcd.get_variable_means();
				data_sigmas = rcd.get_variable_sigmas();
				centered_data_matrix = rcd.get_row_centered_data();

				Matrix Q = new Matrix(ROWS, ROWS);

				for (int i = 0; i < ROWS; i++)
					{
						double s = data_sigmas.get(i, 0);
						double var = (s * s);
						Q.set(i, i, var);
						double mean_X = data_means.get(i, 0);
						Matrix var_X = data_matrix.getMatrix(i, i, 0, COLS - 1);
						for (int j = i + 1; j < ROWS; j++)
							{
								Matrix var_Y = data_matrix.getMatrix(j, j, 0, COLS - 1);
								double[] var_XY = (var_X.arrayTimes(var_Y)).getRowPackedCopy();
								double mean_Y = data_means.get(j, 0);
								double mean_XY = Descriptive_Stats.get_mean(var_XY);
								double cov = (mean_XY - (mean_X * mean_Y));
								Q.set(i, j, cov);
								Q.set(j, i, cov);
							}
					}
				data_matrix = null;
				System.gc();
				return Q;
			}

		/**
		 * Method to calculate the reduced C-matrix (covariance or correlation). Input is a 3Nx3N covariance matrix. Output is a NxN reduced C-matrix
		 * 
		 * @return Returns the Reduced Covariance Matrix, rCOV
		 */
		static Matrix get_reduced_C_matrix(Matrix Q)
			{
				System.out.println("\tGetting the Reduced Covariance Matrix from Covariance Matrix: ");
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
				Q = null;
				System.gc();
				return rCOV;
			}

		/**
		 * Method to calculate the reduced Dynamical matrix (inverse covariance or inverse correlation). Input is a 3Nx3N covariance matrix. Output is a NxN reduced Dynamical
		 * matrix
		 * 
		 * @return Returns the Reduced Covariance Matrix, rCOV
		 */
		static Matrix get_reduced_DYN_matrix(Matrix DYN)
			{
				System.out.println("\tGetting the Reduced Dynamical Matrix from Dynamical Matrix: ");
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
				DYN = null;
				System.gc();
				return rDYN;
			}

		/**
		 * /** General method to calculate a reduced C-matrix (covariance or correlation) with specified number of components, m. Input is a mNxmN covariance matrix. Output is a
		 * NxN reduced C-matrix
		 * 
		 * @return Returns the Reduced Covariance Matrix, rCOV
		 */
		static Matrix get_reduced_C_matrix(Matrix Q, int components)
			{
				System.out.println("\tGetting the Reduced Covariance Matrix from Covariance Matrix: ");
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
				Q = null;
				System.gc();
				return rCOV;
			}

		/**
		 * General method to calculate the reduced Dynamical matrix (inverse covariance or inverse correlation) with specified number of components, m. Input is a mNxmN covariance
		 * matrix, and number of components, m. Output is a NxN reduced Dynamical matrix
		 * 
		 * @return Returns the Reduced Covariance Matrix, rCOV
		 */
		static Matrix get_reduced_DYN_matrix(Matrix DYN, int components)
			{
				System.out.println("\tGetting the Reduced Dynamical Matrix from Dynamical Matrix: ");
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
				DYN = null;
				System.gc();
				return rDYN;
			}

		/**
		 * 
		 * Method to calculate the covariance matrix: AA(T) Note: The data are NOT CENTERED before the calculation: Assumed centered...
		 * 
		 * @return Returns the Covariance Matrix, Q
		 */

		Matrix get_covariance_matrix_NO_CENTERING()
			{
				System.out.println("\tGetting the Covariance Matrix from data, with NO ROW Centering: ");
				Matrix Q = data_matrix.times(data_matrix.transpose());
				// data_matrix = null;
				System.gc();
				Q.timesEquals(Math.pow((COLS - 1), -1));
				return Q;
			}

		/**
		 * Method to calculate the correlation matrix from a covariance matrix, Q.
		 * 
		 * @return Returns the Correlation Matrix, R
		 */
		Matrix get_R_from_Q(Matrix Cov_Matrix)
			{
				System.out.println("\tGetting the Correlation Matrix from specified Covariance Matrix: ");
				int num_Vars = Cov_Matrix.getColumnDimension();
				Matrix R = new Matrix(num_Vars, num_Vars);
				for (int i = 0; i < num_Vars; i++)
					{
						double sigma = Math.sqrt(Cov_Matrix.get(i, i));
						R.set(i, i, 1d);
						for (int j = 0; j < i; j++)
							{
								double element = Cov_Matrix.get(i, j) / (sigma * Math.sqrt(Cov_Matrix.get(j, j)));
								R.set(j, i, element);
								R.set(i, j, element);
							}
					}
				return R;
			}

		/**
		 * Method to calculate the correlation matrix from the data. Note: The data are ROW CENTERED before the calculation.
		 * 
		 * @return Returns the Correlation Matrix, R
		 */
		Matrix get_correlation_matrix()
			{
				System.out.println("\tGetting the Correlation Matrix from data, with ROW Centering: ");
				Matrix COV = get_covariance_matrix_elegant();
				Matrix R = get_R_from_Q(COV);
				return R;
			}

		/**
		 * Method to calculate the correlation matrix from the data. Note: The data are NOT CENTERED before the calculation.
		 * 
		 * @return Returns the Correlation Matrix, R
		 */
		Matrix get_correlation_matrix_NO_CENTERING()
			{
				System.out.println("\tGetting the Correlation Matrix from data, with NO ROW Centering: ");
				Matrix COV = get_covariance_matrix_NO_CENTERING();
				Matrix R = get_R_from_Q(COV);
				return R;
			}

		/**
		 * Method to calculate the partial-correlation matrix from a precision matrix, e.g., inverse covariance matrix.
		 * 
		 * @return Returns the Partial Correlation Matrix, P_CORR
		 */
		static Matrix get_partial_correlation_matrix(Matrix precision)
			{
				System.out.println("\tGetting the Partial Correlation Matrix from the specified Precision Matrix: ");
				int num_Vars = precision.getColumnDimension();
				Matrix P_CORR = new Matrix(num_Vars, num_Vars);
				for (int i = 0; i < num_Vars; i++)
					{
						double x = precision.get(i, i);
						P_CORR.set(i, i, -1d);
						for (int j = 0; j < i; j++)
							{
								double y = precision.get(j, j);
								double r = precision.get(i, j);
								double xy = (x * y);
								double p = 0;
								if (xy > 0)
									p = (-r / (Math.sqrt(xy)));
								if (xy < 0)
									p = (-r / (Math.sqrt(Math.abs(xy))));
								P_CORR.set(j, i, p);
								P_CORR.set(i, j, p);
							}
					}

				return P_CORR;
			}

		/**
		 * Method to calculate the eigenvalue decomposition of a matrix. Note: The matrix argument should be a positive-definite matrix for real eigenvalues. Based in the JAMA
		 * library
		 * 
		 * @param q The matrix to factor
		 * @return The eigenvalue decomposition holding the eigenvalues and eigenvectors.
		 */
		static EigenvalueDecomposition get_eigenvalue_decomposition(Matrix q)
			{
				System.out.println("\tGetting the Eigenvalue Decomposition: ");
				EigenvalueDecomposition evd = new EigenvalueDecomposition(q);
				return evd;
			}

		/**
		 * @return Returns the means of the centered variables.
		 */
		public Matrix getData_means()
			{
				return data_means;
			}

		/**
		 * @return Returns the standard deviations of the centered variables.
		 */
		public Matrix getData_sigmas()
			{
				return data_sigmas;
			}

		public Matrix getCentered_data_Matrix()
			{
				return centered_data_matrix;
			}
	}