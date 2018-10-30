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
 * JED class JED_Get_Distance_Pair_PCA: Gets the COV and CORR PCA for the Residue Distance Pairs subset. Copyright (C) 2012 Dr. Charles David
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Distance_Pair_PCA
	{

		String directory, out_dir, description, file_name_head, path;
		int ROWS_DP, COLS, number_of_modes, number_of_pairs;
		double trace_COV, trace_CORR, trace_PCORR, cond_COV, cond_CORR, cond_PCORR, det_COV, det_CORR, det_PCORR, rank_COV, rank_CORR, rank_PCORR;
		List<Double> eigenvalues_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR, eigenvalues_PCORR, top_eigenvalues_PCORR;
		double[] dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max, dist_pca_mode_PCORR_min, dist_pca_mode_PCORR_max;
		Matrix distances, centered_data, COV_dist, COV_dist_ST, CORR_dist, pcorr_dist, precision_cov, precision_corr, top_evectors_dist_COV, top_evectors_dist_CORR,
				top_evectors_dist_PCORR, residue_means, residue_std_devs;
		NumberFormat nf;
		RoundingMode rm;
		EigenvalueDecomposition evd;
		PCA pca;
		boolean success, exist;

		/* ********************************** CONSTRUCTOR * ************************************************************************* */

		/**
		 * Constructor to perform the Distance Pair PCA analysis
		 *
		 * @param dist  The Distance Pair subset
		 * @param dir   The working directory
		 * @param des   The job description
		 * @param modes The number of PCA modes to process
		 * @param pairs The number of distance pairs
		 */
		JEDi_Get_Distance_Pair_PCA(Matrix dist, String dir, String des, int modes, int pairs)
			{

				this.directory = dir;
				this.description = des;

				this.distances = dist;
				this.number_of_pairs = pairs;
				this.number_of_modes = modes;
				this.ROWS_DP = number_of_pairs;
				this.COLS = number_of_modes;

				nf = NumberFormat.getInstance();
				rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
			}

		/* *********************************** DRIVER METHODS ********************************************************************************** */

		/**
		 * Performs the Distance Pair PCA using both COV and CORR (Q and R) models.
		 */
		public void get_Distance_Pair_PCA()
			{

				pca = new PCA(distances);

				distances = null;
				System.gc();

				Do_Cov("COV");
				Do_Corr("CORR");
				Do_PCorr("PCORR");
			}

		/**
		 * Performs the COV analysis.
		 */
		private void Do_Cov(String model_pca)
			{

				COV_dist = pca.get_covariance_matrix();
				// COV_dist = pca.get_covariance_matrix_elegant();
				COV_dist_ST = pca.stabilize_Covariance_Matrix(COV_dist); // NEW FEATURE
				COV_dist = COV_dist_ST;

				rank_COV = COV_dist.rank();

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + model_pca + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();
				file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";

				residue_means = pca.getData_means();
				path = file_name_head + "_Means_of_Variables.txt";
				Matrix_IO.write_Matrix(residue_means, path, 12, 6);

				residue_std_devs = pca.getData_sigmas();
				path = file_name_head + "_Std_Devs_of_Variables.txt";
				Matrix_IO.write_Matrix(residue_std_devs, path, 12, 6);

				path = file_name_head + "_covariance_matrix.txt";
				Matrix_IO.write_Matrix(COV_dist, path, 12, 6);

				evd = PCA.get_eigenvalue_decomposition(COV_dist);

				get_eigenvalues_COV();
				write_top_evals_COV();
				get_top_evects_and_reverse_COV();
				construct_PCA_Modes_COV();

				System.gc();
			}

		/**
		 * Performs the CORR analysis.
		 */
		private void Do_Corr(String model_pca)
			{

				CORR_dist = pca.get_R_from_Q(COV_dist);
				COV_dist = null;
				System.gc();

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + model_pca + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();

				file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";

				path = file_name_head + "_correlation_matrix.txt";
				Matrix_IO.write_Matrix(CORR_dist, path, 12, 6);

				evd = PCA.get_eigenvalue_decomposition(CORR_dist);

				get_eigenvalues_CORR();
				write_top_evals_CORR();
				get_top_evects_and_reverse_CORR();
				construct_PCA_Modes_CORR();

				System.gc();
			}

		/**
		 * Performs the PCORR analysis.
		 */
		private void Do_PCorr(String model_pca)
			{

				pcorr_dist = PCA.get_partial_correlation_matrix(precision_cov);

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "dpPCA" + File.separatorChar + model_pca + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();

				file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";
				Matrix_IO.write_Matrix(pcorr_dist, file_name_head + "_partial_correlation_matrix.txt", 12, 6);

				evd = PCA.get_eigenvalue_decomposition(pcorr_dist);

				get_eigenvalues_PCORR();
				write_top_evals_PCORR();
				get_top_evects_and_reverse_PCORR();
				construct_PCA_Modes_PCORR();

				precision_cov = null;
				pcorr_dist = null;
				evd = null;
				System.gc();
			}

		/* ************************************* COV METHODS ************************************************************************* ************ */

		/**
		 * Gets the eigenvalues from the COV analysis.
		 */
		private void get_eigenvalues_COV()
			{
				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_COV = new ArrayList<>();
				trace_COV = 0;
				det_COV = 1;
				for (double k : ss_evals)
					{
						if (k < 1.0E-12) k = 0;
						eigenvalues_COV.add(k);
						trace_COV += k;
						det_COV *= k;
					}
				Collections.sort(eigenvalues_COV, Collections.reverseOrder());
				double first = eigenvalues_COV.get(0);
				double last = eigenvalues_COV.get(eigenvalues_COV.size() - 1);
				cond_COV = (first / last);

				path = file_name_head + "_eigenvalues_COV.txt";
				List_IO.write_Double_List(eigenvalues_COV, path, 6);
			}

		/**
		 * Gets the top eigenvalues from the COV analysis.
		 */
		private void write_top_evals_COV()
			{
				try
					{
						path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt";
						File top_ss_evals_cov = new File(path);
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

					}
				catch (IOException io)
					{
						System.err.println("IOException thrown. Could not write the file: " + path);
						io.printStackTrace();
					}

			}

		/**
		 * Gets the top eigenvectors from the COV analysis.
		 */
		private void get_top_evects_and_reverse_COV()
			{

				Matrix ss_evectors = evd.getV();
				Matrix D = evd.getD();
				precision_cov = ss_evectors.times(D.inverse()).times(ss_evectors.inverse());
				evd = null;

				top_evectors_dist_COV = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
				Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_dist_COV.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
					}
				top_evectors_dist_COV = modes_reversed;

				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt";
				Matrix_IO.write_Matrix(top_evectors_dist_COV, path, 12, 6);
				path = file_name_head + "_inverse_covariance_matrix.txt";
				Matrix_IO.write_Matrix(precision_cov, path, 12, 6);
				ss_evectors = null;
			}

		/**
		 * Computes the Distance-Pair PCA modes from the COV analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
		 */
		private void construct_PCA_Modes_COV()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_COV_max = new double[number_of_modes];
				dist_pca_mode_COV_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
							{
								double d = top_evectors_dist_COV.get(b, a);
								double d_abs = Math.abs(d);
								double sq = (d * d);
								double value = top_eigenvalues_COV.get(a);
								if (value < 0) value = 0;
								double sqrt_val = Math.sqrt(value);
								double wm = sqrt_val * d_abs;
								double w_sq = value * sq;
								SS_pca_modes.set(b, a, d_abs);
								SS_weighted_pca_modes.set(b, a, wm);
								SS_pca_square_modes.set(b, a, sq);
								SS_weighted_pca_square_modes.set(b, a, w_sq);
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}

						dist_pca_mode_COV_max[a] = max;
						dist_pca_mode_COV_min[a] = min;
					}
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_COV.txt";
				Array_IO.write_Double_Array(dist_pca_mode_COV_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_COV.txt";
				Array_IO.write_Double_Array(dist_pca_mode_COV_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}

		/* ************************************** CORR METHODS ************************************************************************* ********** */

		/**
		 * Gets the eigenvalues from the CORR analysis.
		 */
		private void get_eigenvalues_CORR()
			{
				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_CORR = new ArrayList<>();
				trace_CORR = 0;
				for (double k : ss_evals)
					{
						if (k < 1.0E-12) k = 0;
						eigenvalues_CORR.add(k);
						trace_CORR += k;
					}
				Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
				path = file_name_head + "_eigenvalues_CORR.txt";
				List_IO.write_Double_List(eigenvalues_CORR, path, 6);
			}

		/**
		 * Gets the top eigenvalues from the CORR analysis.
		 */
		private void write_top_evals_CORR()
			{
				try
					{
						path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR.txt";
						File top_ss_evals_cov = new File(path);
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_CORR = new ArrayList<>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_CORR.get(i);
								double normed_val = (val / trace_CORR);
								cumulative_variance += normed_val;
								top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_CORR.add(val);
							}
						top_ss_evals_writer.close();

					}
				catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the file: " + path);
						io.printStackTrace();
					}
			}

		/**
		 * Gets the top eigenvectors from the CORR analysis.
		 */
		private void get_top_evects_and_reverse_CORR()
			{

				Matrix ss_evectors = evd.getV();
				Matrix D = evd.getD();
				// ss_evectors.times(D.inverse()).times(ss_evectors.transpose());
				precision_corr = ss_evectors.times(D.inverse()).times(ss_evectors.inverse());
				evd = null;
				System.gc();

				top_evectors_dist_CORR = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
				Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
					}
				top_evectors_dist_CORR = modes_reversed;

				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR.txt";
				Matrix_IO.write_Matrix(top_evectors_dist_CORR, path, 12, 6);
				path = file_name_head + "_inverse_correlation_matrix.txt";
				Matrix_IO.write_Matrix(precision_corr, path, 12, 6);
				ss_evectors = null;
				System.gc();
			}

		/**
		 * Computes the Distance-Pair PCA modes from the CORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
		 */
		private void construct_PCA_Modes_CORR()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_CORR_max = new double[number_of_modes];
				dist_pca_mode_CORR_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
							{
								double d = top_evectors_dist_CORR.get(b, a);
								double d_abs = Math.abs(d);
								double sq = (d * d);
								double value = top_eigenvalues_CORR.get(a);
								if (value < 0) value = 0;
								double sqrt_val = Math.sqrt(value);
								double wm = sqrt_val * d_abs;
								double w_sq = value * sq;
								SS_pca_modes.set(b, a, d_abs);
								SS_weighted_pca_modes.set(b, a, wm);
								SS_pca_square_modes.set(b, a, sq);
								SS_weighted_pca_square_modes.set(b, a, w_sq);
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}

						dist_pca_mode_CORR_max[a] = max;
						dist_pca_mode_CORR_min[a] = min;
					}
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_CORR_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_CORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_CORR_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}

		/* *********************************** PCORR METHODS ******************************************************************** */

		/**
		 * Gets the eigenvalues from the PCORR analysis.
		 */
		private void get_eigenvalues_PCORR()
			{

				double[] ss_evals = evd.getRealEigenvalues();
				eigenvalues_PCORR = new ArrayList<>();
				for (double k : ss_evals)
					{
						if (Math.abs(k) < 1.0E-9) k = 0;
						eigenvalues_PCORR.add(k);
					}
				Collections.sort(eigenvalues_PCORR, Collections.reverseOrder());
				List_IO.write_Double_List(eigenvalues_PCORR, file_name_head + "_eigenvalues_PCORR.txt", 6);
			}

		/**
		 * Gets the top eigenvalues from the PCORR analysis.
		 */
		private void write_top_evals_PCORR()
			{
				try
					{
						File top_ss_evals_cov = new File(file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR.txt");
						BufferedWriter top_ss_evals_writer = new BufferedWriter(new FileWriter(top_ss_evals_cov));
						top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));
						top_eigenvalues_PCORR = new ArrayList<>();
						double cumulative_variance = 0;
						for (int i = 0; i < number_of_modes; i++)
							{
								double val = eigenvalues_PCORR.get(i);
								double normed_val = (val / trace_PCORR);
								cumulative_variance += normed_val;
								top_ss_evals_writer.write(String.format("%-16s%-16s%-16s%n", nf.format(val), nf.format(normed_val), nf.format(cumulative_variance)));
								top_eigenvalues_PCORR.add(val);
							}
						top_ss_evals_writer.close();
					}
				catch (IOException io)
					{
						System.err.println("Could not write to the file: " + file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR.txt");
						io.getMessage();
						io.getStackTrace();
					}
			}

		/**
		 * Gets the top eigenvectors from the PCORR analysis.
		 */
		private void get_top_evects_and_reverse_PCORR()
			{

				Matrix ss_evectors = evd.getV();

				evd = null;
				System.gc();

				top_evectors_dist_PCORR = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
				Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
				for (int r = 0; r < COLS; r++)
					{
						Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
						modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
					}
				top_evectors_dist_CORR = modes_reversed;

				path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR.txt";
				Matrix_IO.write_Matrix(top_evectors_dist_PCORR, path, 12, 6);
				ss_evectors = null;
				System.gc();
			}

		/**
		 * Computes the Distance-Pair PCA modes from the PCORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
		 */
		private void construct_PCA_Modes_PCORR()
			{

				Matrix SS_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_modes = new Matrix(ROWS_DP, number_of_modes);
				Matrix SS_weighted_pca_square_modes = new Matrix(ROWS_DP, number_of_modes);
				dist_pca_mode_PCORR_max = new double[number_of_modes];
				dist_pca_mode_PCORR_min = new double[number_of_modes];

				for (int a = 0; a < number_of_modes; a++)
					{
						double max = 0;
						double min = 1;

						for (int b = 0; b < ROWS_DP; b++)
							{
								double d = top_evectors_dist_PCORR.get(b, a);
								double d_abs = Math.abs(d);
								double sq = (d * d);
								double value = top_eigenvalues_PCORR.get(a);
								double value_abs = Math.abs(value);
								double sqrt_val = Math.sqrt(value_abs);
								double wm = sqrt_val * d_abs;
								double w_sq = value_abs * sq;
								SS_pca_modes.set(b, a, d_abs);
								SS_weighted_pca_modes.set(b, a, wm);
								SS_pca_square_modes.set(b, a, sq);
								SS_weighted_pca_square_modes.set(b, a, w_sq);
								if (d_abs >= max) max = d_abs;
								if (d_abs <= min) min = d_abs;
							}
						dist_pca_mode_PCORR_max[a] = max;
						dist_pca_mode_PCORR_min[a] = min;
					}
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MAXES_PCORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_PCORR_max, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_mode_MINS_PCORR.txt";
				Array_IO.write_Double_Array(dist_pca_mode_PCORR_min, path, 6);
				path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR.txt";
				Matrix_IO.write_Matrix(SS_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR.txt";
				Matrix_IO.write_Matrix(SS_pca_square_modes, path, 12, 6);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR.txt";
				Matrix_IO.write_Matrix(SS_weighted_pca_square_modes, path, 12, 6);
			}
		/* ************************************** GETTERS ************************************************************************* **************** */

		/**
		 * @return the directory
		 */
		public String getDirectory()
			{
				return directory;
			}

		/**
		 * @return the description
		 */
		public String getDescription()
			{
				return description;
			}

		/**
		 * @return the file_name_head
		 */
		public String getFile_name_head()
			{
				return file_name_head;
			}

		/**
		 * @return the path
		 */
		public String getPath()
			{
				return path;
			}

		/**
		 * @return the ROWS_DP
		 */
		public int getROWS_DP()
			{
				return ROWS_DP;
			}

		/**
		 * @return the COLS
		 */
		public int getCOLS()
			{
				return COLS;
			}

		/**
		 * @return the number_of_modes
		 */
		public int getNumber_of_modes()
			{
				return number_of_modes;
			}

		/**
		 * @return the number_of_pairs
		 */
		public int getNumber_of_pairs()
			{
				return number_of_pairs;
			}

		/**
		 * @return the trace_COV
		 */
		public double getTrace_COV()
			{
				return trace_COV;
			}

		/**
		 * @return the trace_CORR
		 */
		public double getTrace_CORR()
			{
				return trace_CORR;
			}

		/**
		 * @return the trace_PCORR
		 */
		public double getTrace_PCORR()
			{
				return trace_PCORR;
			}

		/**
		 * @return the cond_COV
		 */
		public double getCond_COV()
			{
				return cond_COV;
			}

		/**
		 * @return the cond_CORR
		 */
		public double getCond_CORR()
			{
				return cond_CORR;
			}

		/**
		 * @return the cond_PCORR
		 */
		public double getCond_PCORR()
			{
				return cond_PCORR;
			}

		/**
		 * @return the det_COV
		 */
		public double getDet_COV()
			{
				return det_COV;
			}

		/**
		 * @return the det_CORR
		 */
		public double getDet_CORR()
			{
				return det_CORR;
			}

		/**
		 * @return the det_PCORR
		 */
		public double getDet_PCORR()
			{
				return det_PCORR;
			}

		/**
		 * @return the rank_COV
		 */
		public double getRank_COV()
			{
				return rank_COV;
			}

		/**
		 * @return the rank_CORR
		 */
		public double getRank_CORR()
			{
				return rank_CORR;
			}

		/**
		 * @return the rank_PCORR
		 */
		public double getRank_PCORR()
			{
				return rank_PCORR;
			}

		/**
		 * @return the eigenvalues_COV
		 */
		public List<Double> getEigenvalues_COV()
			{
				return eigenvalues_COV;
			}

		/**
		 * @return the top_eigenvalues_COV
		 */
		public List<Double> getTop_eigenvalues_COV()
			{
				return top_eigenvalues_COV;
			}

		/**
		 * @return the eigenvalues_CORR
		 */
		public List<Double> getEigenvalues_CORR()
			{
				return eigenvalues_CORR;
			}

		/**
		 * @return the top_eigenvalues_CORR
		 */
		public List<Double> getTop_eigenvalues_CORR()
			{
				return top_eigenvalues_CORR;
			}

		/**
		 * @return the eigenvalues_PCORR
		 */
		public List<Double> getEigenvalues_PCORR()
			{
				return eigenvalues_PCORR;
			}

		/**
		 * @return the top_eigenvalues_PCORR
		 */
		public List<Double> getTop_eigenvalues_PCORR()
			{
				return top_eigenvalues_PCORR;
			}

		/**
		 * @return the dist_pca_mode_COV_min
		 */
		public double[] getDist_pca_mode_COV_min()
			{
				return dist_pca_mode_COV_min;
			}

		/**
		 * @return the dist_pca_mode_COV_max
		 */
		public double[] getDist_pca_mode_COV_max()
			{
				return dist_pca_mode_COV_max;
			}

		/**
		 * @return the dist_pca_mode_CORR_min
		 */
		public double[] getDist_pca_mode_CORR_min()
			{
				return dist_pca_mode_CORR_min;
			}

		/**
		 * @return the dist_pca_mode_CORR_max
		 */
		public double[] getDist_pca_mode_CORR_max()
			{
				return dist_pca_mode_CORR_max;
			}

		/**
		 * @return the dist_pca_mode_PCORR_min
		 */
		public double[] getDist_pca_mode_PCORR_min()
			{
				return dist_pca_mode_PCORR_min;
			}

		/**
		 * @return the dist_pca_mode_PCORR_max
		 */
		public double[] getDist_pca_mode_PCORR_max()
			{
				return dist_pca_mode_PCORR_max;
			}

		/**
		 * @return the distances
		 */
		public Matrix getDistances()
			{
				return distances;
			}

		/**
		 * @return the centered_data
		 */
		public Matrix getCentered_data()
			{
				return centered_data;
			}

		/**
		 * @return the cOV_dist
		 */
		public Matrix getCOV_dist()
			{
				return COV_dist;
			}

		/**
		 * @return the cORR_dist
		 */
		public Matrix getCORR_dist()
			{
				return CORR_dist;
			}

		/**
		 * @return the pcorr_dist
		 */
		public Matrix getPcorr_dist()
			{
				return pcorr_dist;
			}

		/**
		 * @return the precision_cov
		 */
		public Matrix getPrecision_cov()
			{
				return precision_cov;
			}

		/**
		 * @return the precision_corr
		 */
		public Matrix getPrecision_corr()
			{
				return precision_corr;
			}

		/**
		 * @return the top_evectors_dist_COV
		 */
		public Matrix getTop_evectors_dist_COV()
			{
				return top_evectors_dist_COV;
			}

		/**
		 * @return the top_evectors_dist_CORR
		 */
		public Matrix getTop_evectors_dist_CORR()
			{
				return top_evectors_dist_CORR;
			}

		/**
		 * @return the top_evectors_dist_PCORR
		 */
		public Matrix getTop_evectors_dist_PCORR()
			{
				return top_evectors_dist_PCORR;
			}

		/**
		 * @return the residue_means
		 */
		public Matrix getResidue_means()
			{
				return residue_means;
			}

		/**
		 * @return the residue_std_devs
		 */
		public Matrix getResidue_std_devs()
			{
				return residue_std_devs;
			}

		/**
		 * @return the nf
		 */
		public NumberFormat getNf()
			{
				return nf;
			}

		/**
		 * @return the rm
		 */
		public RoundingMode getRm()
			{
				return rm;
			}

		/**
		 * @return the evd
		 */
		public EigenvalueDecomposition getEvd()
			{
				return evd;
			}

		/**
		 * @return the pca
		 */
		public PCA getPca()
			{
				return pca;
			}

		/**
		 * @return the success
		 */
		public boolean isSuccess()
			{
				return success;
			}

		/**
		 * @return the exist
		 */
		public boolean isExist()
			{
				return exist;
			}
	}
