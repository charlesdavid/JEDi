package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * JED class JED_Get_Generalized_Coordinate_PCA: Gets the PCA using COV for the Cartesian subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Get_Hierarchical_PCA
	{

		String directory, out_dir, description, file_name_head, path;
		int number_of_modes, number_of_residues, ROWS, COLS, number_of_modes_Residue;
		double trace_COV, cond_COV, rank_COV, det_COV;
		List<Double> eigenvalues_COV, top_eigenvalues_COV;
		double[] pca_mode_COV_min, pca_mode_COV_max;
		Matrix input_data, centered_input_data, cov, cov_ST, rcov, top_evectors_COV, square_pca_modes_COV, weighted_square_pca_modes_COV, weighted_pca_modes_COV, pca_modes_COV,
				data_Means, data_Sigmas;
		EigenvalueDecomposition evd;
		NumberFormat nf;
		RoundingMode rm;
		PCA pca;
		boolean success, exist;

		/**
		 * Constructor to perform PCA on generalized residue coordinates
		 *
		 * @param data      The subset coordinates
		 * @param num_modes The number of PCA modes to process
		 * @param dir       The working directory
		 * @param des       The job description
		 */
		JEDi_Get_Hierarchical_PCA(Matrix data, int num_modes, String dir, String des, int res_modes)
			{

				nf = NumberFormat.getInstance();
				rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(12);
				nf.setMinimumFractionDigits(12);

				this.input_data = data;
				this.number_of_modes = num_modes;
				this.directory = dir;
				this.description = des;
				this.number_of_modes_Residue = res_modes;

				ROWS = input_data.getRowDimension();
				COLS = number_of_modes;
				number_of_residues = (ROWS / number_of_modes_Residue);

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Hierarchical_PCA" + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();

				file_name_head = out_dir + "ss_" + number_of_residues;
			}

		/* **************************************************** PUBLIC METHOD ****************************************************** */

		/**
		 * Performs analysis using the COV PCA model.
		 */
		public void get_Hierarchical_PCA()
			{
				pca = new PCA(input_data);
				Do_Cov_Model();
			}

		/* **************************************************** PRIVATE METHODS ****************************************************** */

		private void Do_Cov_Model()
			{
				cov = pca.get_covariance_matrix();
				cov_ST = pca.stabilize_Covariance_Matrix(cov); // ********************************NEW
				cov = cov_ST;

				// cov = pca.get_covariance_matrix_elegant();
				centered_input_data = pca.getCentered_data_Matrix();
				// rcov = PCA.get_reduced_DYN_matrix(cov);
				rank_COV = cov.rank();
				data_Means = pca.getData_means();
				data_Sigmas = pca.getData_sigmas();


				path = file_name_head + "_means_of_variables.txt";
				Matrix_IO.write_Matrix(data_Means, path, 12, 6);
				path = file_name_head + "_std_devs_of_variables.txt";
				Matrix_IO.write_Matrix(data_Sigmas, path, 12, 6);


				Matrix_IO.write_Matrix(cov, file_name_head + "_covariance_matrix.txt", 12, 6);
				// Matrix_IO.write_Matrix(rcov, file_name_head + "_reduced_covariance_matrix.txt", 12, 6);

				evd = PCA.get_eigenvalue_decomposition(cov);
				get_eigenvalues_COV();
				write_top_evals_COV();
				get_top_evects_and_reverse_COV();
				construct_PCA_Modes_COV();

				evd = null;
				System.gc();
			}

		/* **************************************** COV METHODS ************************************************************* */

		private void get_eigenvalues_COV()
			{

				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_COV = new ArrayList<>();
				trace_COV = 0;
				det_COV = 1;
				for (double k : ss_evals)
					{
						if (Math.abs(k) < 0.000000000001)
							k = 0;
						eigenvalues_COV.add(k);
						trace_COV += k;
						det_COV *= k;
					}
				Collections.sort(eigenvalues_COV, Collections.reverseOrder());
				int size = eigenvalues_COV.size();
				double first = eigenvalues_COV.get(0);
				double last = eigenvalues_COV.get(size - 1);
				cond_COV = (first / last);
				List_IO.write_Double_List(eigenvalues_COV, file_name_head + "_eigenvalues_COV.txt", 6);
			}

		private void write_top_evals_COV()
			{
				try
					{
						File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt");
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_COV = new ArrayList<>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_COV.get(i);
								double normed_val = (val / trace_COV);
								cumulative_variance += normed_val;
								top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_COV.add(val);
							}
						top_ss_evals_writer.close();
					} catch (IOException io)
					{
						System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt");
						io.getMessage();
						io.getStackTrace();
					}
			}

		private void get_top_evects_and_reverse_COV()
			{

				Matrix ss_evectors = evd.getV();
				top_evectors_COV = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
				Matrix modes_reversed = new Matrix(ROWS, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_COV.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
					}
				top_evectors_COV = modes_reversed;

				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt";
				Matrix_IO.write_Matrix(top_evectors_COV, path, 12, 9);

				ss_evectors = null;
				System.gc();
			}

		private void construct_PCA_Modes_COV()
			{

				pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
				weighted_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
				square_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
				weighted_square_pca_modes_COV = new Matrix(number_of_residues, number_of_modes);
				pca_mode_COV_max = new double[number_of_modes];
				pca_mode_COV_min = new double[number_of_modes];
				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < number_of_residues; b++)
							{
								int offset = b * number_of_modes_Residue;
								// System.out.println("offset: " + offset);
								double sum_of_Squares = 0;
								for (int k = 0; k < number_of_modes_Residue; k++)
									{
										double q = top_evectors_COV.get(offset + k, a);
										double sq = q * q;
										sum_of_Squares += sq;
									}
								double sqrt_sq = Math.sqrt(sum_of_Squares);
								double value = eigenvalues_COV.get(a);
								double sqrt_val = Math.sqrt(value);
								double w_sq = sum_of_Squares * value;
								double w_m = sqrt_sq * sqrt_val;
								pca_modes_COV.set(b, a, sqrt_sq);
								weighted_pca_modes_COV.set(b, a, w_m);
								square_pca_modes_COV.set(b, a, sum_of_Squares);
								weighted_square_pca_modes_COV.set(b, a, w_sq);
								if (sum_of_Squares > max)
									max = sum_of_Squares;
								if (sum_of_Squares < min)
									min = sum_of_Squares;
							}
						pca_mode_COV_max[a] = max;
						pca_mode_COV_min[a] = min;
					}
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_COV.txt";
				Array_IO.write_Double_Array(pca_mode_COV_max, path, 3);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_COV.txt";
				Array_IO.write_Double_Array(pca_mode_COV_min, path, 3);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(pca_modes_COV, path, 12, 3);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(weighted_pca_modes_COV, path, 12, 3);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(square_pca_modes_COV, path, 12, 3);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(weighted_square_pca_modes_COV, path, 12, 3);

			}

		/* ******************************************* GETTERS ***************************************************** */

		public double get_cond_COV()
			{

				return cond_COV;
			}

		public double get_trace_COV()
			{

				return trace_COV;
			}

		public List<Double> getEigenvalues_COV()
			{

				return eigenvalues_COV;
			}

		public double get_det_COV()
			{

				return det_COV;
			}

		public double get_rank_COV()
			{

				return rank_COV;
			}

		public double[] get_pca_mode_COV_min()
			{

				return pca_mode_COV_min;
			}

		public double[] get_pca_mode_COV_max()
			{

				return pca_mode_COV_max;
			}

		public Matrix getTop_evectors_COV()
			{

				return top_evectors_COV;
			}

		public Matrix getSquare_pca_modes_COV()
			{

				return square_pca_modes_COV;
			}

		public Matrix getWeighted_square_pca_modes_COV()
			{

				return weighted_square_pca_modes_COV;
			}

		public Matrix getWeighted_pca_modes_COV()
			{

				return weighted_pca_modes_COV;
			}

		public Matrix getPca_modes_COV()
			{

				return pca_modes_COV;
			}

		public List<Double> getTop_eigenvalues_COV()
			{

				return top_eigenvalues_COV;
			}

		public Matrix getData_Means()
			{
				return data_Means;
			}

		public Matrix getData_Sigmas()
			{
				return data_Sigmas;
			}

		public Matrix getCentered_input_data()
			{
				return centered_input_data;
			}
	}
