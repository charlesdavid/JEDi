package jedi;

import java.io.File;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import support.Collectivity;
import support.PCA;
import supportIO.Input_Parameters;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportIO.Write_Top_Eigenvalues;
import supportPlot.COLLECTIVITY_Plot;
import supportPlot.HeatMap;
import supportPlot.SCREE_Plot;

/**
 * JED class JED_Get_Cartesian_PCA: Gets the PCA using COV and/or CORR and PCORR for the Cartesian subset. Copyright (C) 2012 Dr. Charles David
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
public class JEDi_Get_Cartesian_PCA
{
	boolean success, exist;
	final boolean doCORR, doPCORR, doReduce, doSparsify, verbose;
	final String directory, description, type;
	String out, out_dir, file_name_head, path, model;
	final int number_of_modes, number_of_atoms, ROWS, COLS;
	int number_of_residues, rank;
	double shrinkage, trace_COV, trace_CORR, trace_PCORR, trace_CORR_SPARSE, trace_PCORR_SPARSE, cond, det, KMO, percent_sparse_COV, percent_sparse_CORR, percent_sparse_PCORR;
	final double st_threshold, FLOOR, NOISE_LEVEL;
	List<Double> eigenvalues_COV, eigenvalues_prime_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR, eigenvalues_PCORR, top_eigenvalues_PCORR;
	List<Double> eigenvalues_CORR_SPARSE, top_eigenvalues_CORR_SPARSE, eigenvalues_PCORR_SPARSE, top_eigenvalues_PCORR_SPARSE;
	double[] pca_mode_mins_COV, pca_mode_maxes_COV, pca_mode_mins_CORR, pca_mode_maxes_CORR, pca_mode_mins_PCORR, pca_mode_maxes_PCORR;
	double[] pca_mode_mins_CORR_SPARSE, pca_mode_maxes_CORR_SPARSE, pca_mode_mins_PCORR_SPARSE, pca_mode_maxes_PCORR_SPARSE;
	Matrix centered_input_data, cov, rcov, inv_cov, r_inv_cov, corr, rcorr, pcorr, rpcorr, corr_sparse, rcorr_sparse, pcorr_sparse, rpcorr_sparse, diff_RP, rdiff_RP;
	Matrix top_evectors_COV, square_pca_modes_COV, weighted_square_pca_modes_COV, weighted_pca_modes_COV, top_evectors_CORR, square_pca_modes_CORR, weighted_square_pca_modes_CORR,
			weighted_pca_modes_CORR, pca_modes_COV, pca_modes_CORR, top_evectors_PCORR, square_pca_modes_PCORR, weighted_square_pca_modes_PCORR, weighted_pca_modes_PCORR,
			pca_modes_PCORR;
	Matrix top_evectors_CORR_SPARSE, square_pca_modes_CORR_SPARSE, weighted_pca_modes_CORR_SPARSE, pca_modes_CORR_SPARSE, top_evectors_PCORR_SPARSE, square_pca_modes_PCORR_SPARSE,
			weighted_pca_modes_PCORR_SPARSE, pca_modes_PCORR_SPARSE, weighted_square_pca_modes_CORR_SPARSE, weighted_square_pca_modes_PCORR_SPARSE;
	Matrix collectivity_cov, collectivity_corr, collectivity_pcorr, collectivity_corr_SPARSE, collectivity_pcorr_SPARSE;
	final Matrix input_data;
	EigenvalueDecomposition evd_cov, evd_corr, evd_pcorr, evd_corr_sparse, evd_pcorr_sparse;
	final NumberFormat nf, nf6, nf12;
	final RoundingMode rm;
	final PCA pca;

	/**
	 * Constructor to perform the Cartesian PCA
	 *
	 * @param data      The Cartesian subset coordinates
	 * @param num_modes The number of PCA modes to process
	 * @param dir       The working directory
	 * @param des       The job description
	 * @param type_pca  The type of PCA analysis
	 */
	JEDi_Get_Cartesian_PCA(Matrix data, int num_modes, String type_pca)
	{
		this.nf = NumberFormat.getInstance();
		this.rm = RoundingMode.HALF_UP;
		this.nf.setRoundingMode(rm);
		this.nf.setMaximumFractionDigits(0);
		this.nf.setMinimumFractionDigits(0);

		this.nf6 = NumberFormat.getInstance();
		this.nf6.setRoundingMode(rm);
		this.nf6.setMaximumFractionDigits(6);
		this.nf6.setMinimumFractionDigits(6);

		this.nf12 = NumberFormat.getInstance();
		this.nf12.setRoundingMode(rm);
		this.nf12.setMaximumFractionDigits(12);
		this.nf12.setMinimumFractionDigits(12);

		this.directory = Input_Parameters.DIRECTORY;
		this.description = Input_Parameters.DESCRIPTION;
		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doReduce = Input_Parameters.doREDUCE;
		this.doSparsify = Input_Parameters.doSPARSIFY;
		this.FLOOR = Input_Parameters.FLOOR;
		this.NOISE_LEVEL = Input_Parameters.NOISE_LEVEL;
		this.st_threshold = Input_Parameters.THRESHOLD_COV;
		this.verbose = Input_Parameters.verbose;

		this.number_of_modes = num_modes;
		this.input_data = data;
		this.type = type_pca;
		this.pca = new PCA(input_data);

		this.ROWS = input_data.getRowDimension();
		this.COLS = number_of_modes;
		this.number_of_atoms = (ROWS / 3);
	}

	/* *************************************** PRIMARY METHODS ****************************************************** */

	public void Do_Cov_PCA()
	{
		model = "COV";
		cov = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		centered_input_data = pca.getCentered_data_Matrix();
		shrinkage = pca.get_lambda_cov();

		if (st_threshold > 0)
			{
				Matrix cov_st = pca.Threshold_COV_Matrix(cov);
				cov = cov_st;
			}

		out = out_dir + model + File.separatorChar;
		create_Directory(out);

		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(cov, file_name_head + "_covariance_matrix.txt.bz2");
		HeatMap hm_cov = new HeatMap(cov, out, "Covariance_Matrix", "Variable 1", "Variable 2");
		hm_cov.createMatrixHeatmap();

		if (doReduce)
			{
				rcov = pca.get_Reduced_C_Matrix(cov);
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(rcov, file_name_head + "_reduced_covariance_matrix.txt.bz2");
				HeatMap hm_rcov = new HeatMap(rcov, out, "Reduced_Covariance_Matrix", "Atom 1", "Atom 2");
				hm_rcov.createMatrixHeatmap();
			}
		evd_cov = pca.get_eigenvalue_decomposition(cov);
		get_eigenvalues_COV();
		write_top_evals_COV();
		get_top_evects_and_reverse_COV();
		construct_PCA_Modes_COV();

		corr = pca.get_Correlation_Matrix_from_Covariance_Matrix(cov);
		pcorr = pca.get_Partial_Correlation_Matrix(inv_cov);
	}

	public void Do_Corr_PCA()
	{
		model = "CORR";
		out = out_dir + model + File.separatorChar;
		create_Directory(out);
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(corr, file_name_head + "_correlation_matrix.txt.bz2");
		HeatMap hm_corr = new HeatMap(corr, out, "Correlation_Matrix", "Variable 1", "Variable 2");
		hm_corr.createMatrixHeatmap();

		if (doReduce)
			{
				rcorr = pca.get_Reduced_C_Matrix(corr);
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(rcorr, file_name_head + "_reduced_correlation_matrix.txt.bz2");
				HeatMap hm_rcorr = new HeatMap(rcorr, out, "Reduced_Correlation_Matrix", "Atom 1", "Atom 2");
				hm_rcorr.createMatrixHeatmap();
			}

		evd_corr = pca.get_eigenvalue_decomposition(corr);

		get_eigenvalues_CORR();
		write_top_evals_CORR();
		get_top_evects_and_reverse_CORR();
		construct_PCA_Modes_CORR();
	}

	public void Do_CORR_SPARSE_PCA()
	{
		out = out_dir + model + File.separatorChar + "sparse" + File.separatorChar;
		create_Directory(out);

		corr_sparse = pca.Sparsify_Matrix(corr, Input_Parameters.THRESHOLD_CORR);
		percent_sparse_CORR = pca.get_Percent_Sparse();

		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(corr_sparse, file_name_head + "_sparsified_correlation_matrix.txt.bz2");
		String name = "Sparsified_Correlation_Matrix_" + nf.format(percent_sparse_CORR) + "_Percent";
		HeatMap hm_scorr = new HeatMap(corr_sparse, out, name, "Variable 1", "Variable 2");
		hm_scorr.createMatrixHeatmap();

		if (doReduce)
			{
				rcorr_sparse = pca.get_Reduced_C_Matrix(corr_sparse);
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(rcorr_sparse, file_name_head + "_reduced_sparsified_correlation_matrix.txt.bz2");
				HeatMap hm_rcorr = new HeatMap(rcorr_sparse, out, "Reduced_Sparsified_Correlation_Matrix", "Atom 1", "Atom 2");
				hm_rcorr.createMatrixHeatmap();
			}

		evd_corr_sparse = pca.get_eigenvalue_decomposition(corr_sparse);

		get_eigenvalues_CORR_SPARSE();
		write_top_evals_CORR_SPARSE();
		get_top_evects_and_reverse_CORR_SPARSE();
		construct_PCA_Modes_CORR_SPARSE();
	}

	public void Do_PCorr_PCA()
	{
		model = "PCORR";
		out = out_dir + model + File.separatorChar;
		create_Directory(out);

		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pcorr, file_name_head + "_partial_correlation_matrix.txt.bz2");
		HeatMap hm_pcorr = new HeatMap(pcorr, out, "Partial_Correlation_Matrix", "Variable 1", "Variable 2");
		hm_pcorr.createMatrixHeatmap();

		if (doReduce)
			{
				rpcorr = pca.get_Reduced_C_Matrix(pcorr);
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(rpcorr, file_name_head + "_Reduced_Partial_Correlation_Matrix.txt.bz2");
				HeatMap hm_rpcorr = new HeatMap(rpcorr, out, "Reduced_partial_correlation_matrix", "Atom 1", "Atom 2");
				hm_rpcorr.createMatrixHeatmap();
			}

		evd_pcorr = pca.get_eigenvalue_decomposition(pcorr);

		get_eigenvalues_PCORR();
		write_top_evals_PCORR();
		get_top_evects_and_reverse_PCORR();
		construct_PCA_Modes_PCORR();
	}

	public void Do_PCORR_SPARSE_PCA()
	{
		out = out_dir + model + File.separatorChar + "sparse" + File.separatorChar;
		create_Directory(out);
		exist = new File(out).exists();

		pcorr_sparse = pca.Sparsify_Matrix(pcorr, Input_Parameters.THRESHOLD_PCORR);
		percent_sparse_PCORR = pca.get_Percent_Sparse();
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pcorr_sparse, file_name_head + "_sparsified_partial_correlation_matrix.txt.bz2");
		String name = "Sparsified_Partial_Correlation_Matrix_" + nf.format(percent_sparse_PCORR) + "_Percent";
		HeatMap hm_pcorrS = new HeatMap(pcorr_sparse, out, name, "Variable 1", "Variable 2");
		hm_pcorrS.createMatrixHeatmap();

		if (doReduce)
			{
				rpcorr_sparse = pca.get_Reduced_C_Matrix(pcorr_sparse);
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(rpcorr_sparse, file_name_head + "_Reduced_Sparsified_Partial_Correlation_Matrix.txt.bz2");
				HeatMap hm_rpcorrS = new HeatMap(rpcorr_sparse, out, "Reduced_Sparsified_Partial_Correlation_Matrix", "Atom 1", "Atom 2");
				hm_rpcorrS.createMatrixHeatmap();

				Matrix rdiff = rcorr_sparse.minus(rpcorr_sparse);
				Matrix rdiff_abs = new Matrix(rcorr_sparse.getRowDimension(), rcorr_sparse.getRowDimension());
				Matrix radjacency = new Matrix(rcorr_sparse.getRowDimension(), rcorr_sparse.getRowDimension());

				for (int i = 0; i < rcorr_sparse.getRowDimension(); i++)
					{
						for (int j = 0; j < i; j++)
							{
								double r = rcorr_sparse.get(i, j);
								double r_abs = Math.abs(r);
								double p = rpcorr_sparse.get(i, j);
								double p_abs = Math.abs(p);

								rdiff_abs.set(i, j, (r_abs - p_abs));
								rdiff_abs.set(j, i, (r_abs - p_abs));
							}
					}

				HeatMap hm = new HeatMap(rdiff, out, "Difference of Reduced Correlation and  Partial Correlation Matrices", "Variable 1", "Variable 2");
				hm.createMatrixHeatmap();

				for (int i = 0; i < rcorr_sparse.getRowDimension(); i++)
					{
						for (int j = 0; j < i; j++)
							{
								double difference = rdiff.get(i, j);
								double val = rdiff_abs.get(i, j);
								if (val > Input_Parameters.THRESHOLD_PCORR & difference > 0) // This means that corr > pcorr --> activation
									{
										radjacency.set(i, j, 1d);
									}
								if (val > Input_Parameters.THRESHOLD_PCORR & difference < 0) // This means that corr < pcorr --> suppression
									{
										radjacency.set(i, j, -1d);
									}
							}
					}

				hm = new HeatMap(radjacency, out, "Adjacency Matrix derived from difference of Reduced Correlation and  Partial Correlation Matrices", "Atom 1", "Atom 2");
				hm.createAdjacencyHeatmap();
			}

		Matrix diff = corr_sparse.minus(pcorr_sparse);
		Matrix diff_abs = new Matrix(corr_sparse.getRowDimension(), corr_sparse.getRowDimension());
		Matrix adjacency = new Matrix(corr_sparse.getRowDimension(), corr_sparse.getRowDimension());

		for (int i = 0; i < corr_sparse.getRowDimension(); i++)
			{
				for (int j = 0; j < i; j++)
					{
						double r = corr_sparse.get(i, j);
						double r_abs = Math.abs(r);
						double p = pcorr_sparse.get(i, j);
						double p_abs = Math.abs(p);

						diff_abs.set(i, j, (r_abs - p_abs));
						diff_abs.set(j, i, (r_abs - p_abs));
					}
			}

		HeatMap hm = new HeatMap(diff, out, "Difference of Correlation and  Partial Correlation Matrices", "Variable 1", "Variable 2");
		hm.createMatrixHeatmap();

		for (int i = 0; i < corr_sparse.getRowDimension(); i++)
			{
				for (int j = 0; j < i; j++)
					{
						double difference = diff.get(i, j);
						double val = diff_abs.get(i, j);
						if (val > Input_Parameters.THRESHOLD_PCORR & difference > 0) // This means that corr > pcorr --> activation
							{
								adjacency.set(i, j, 1d);
							}
						if (val > Input_Parameters.THRESHOLD_PCORR & difference < 0) // This means that corr < pcorr --> suppression
							{
								adjacency.set(i, j, -1d);
							}
					}
			}

		hm = new HeatMap(adjacency, out, "Adjacency Matrix derived from difference of Correlation and Partial Correlation Matrices", "Variable 1", "Variable 2");
		hm.createAdjacencyHeatmap();

		evd_pcorr_sparse = pca.get_eigenvalue_decomposition(pcorr_sparse);

		get_eigenvalues_PCORR_SPARSE();
		write_top_evals_PCORR_SPARSE();
		get_top_evects_and_reverse_PCORR_SPARSE();
		construct_PCA_Modes_PCORR_SPARSE();
	}

	/* **************************************** COV METHODS ************************************************************* */

	private void get_eigenvalues_COV()
	{
		double[] ss_evals = evd_cov.getRealEigenvalues();
		eigenvalues_COV = new ArrayList<Double>();
		eigenvalues_prime_COV = new ArrayList<Double>();
		trace_COV = 0;
		det = 1;
		int rankCount = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				if (k > 0) rankCount++;
				eigenvalues_COV.add(k);
				trace_COV += k;
				det *= k;
				if (k < NOISE_LEVEL) eigenvalues_prime_COV.add(NOISE_LEVEL);
				else
					eigenvalues_prime_COV.add(k);
			}
		rank = rankCount;
		Collections.sort(eigenvalues_COV, Collections.reverseOrder());
		Collections.sort(eigenvalues_prime_COV, Collections.reverseOrder());

		int size = eigenvalues_COV.size();
		double first = eigenvalues_COV.get(0);
		double last = eigenvalues_COV.get(size - 1);
		cond = (first / last);

		if (verbose) List_IO.write_Double_List(eigenvalues_COV, file_name_head + "_eigenvalues.txt", 12);
		if (verbose) List_IO.write_Double_List(eigenvalues_prime_COV, file_name_head + "_eigenvalues_prime.txt", 12);
	}

	private void write_top_evals_COV()
	{
		Matrix screePlot = new Matrix(number_of_modes, 2);
		top_eigenvalues_COV = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes; i++)
			{
				double val = eigenvalues_COV.get(i);
				double normed_val = (val / trace_COV);
				cumulative_variance += normed_val;
				top_eigenvalues_COV.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_COV.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues_COV, trace_COV);
		SCREE_Plot.createChart(out, "Scree_Plot_COV", screePlot);
	}

	private void get_top_evects_and_reverse_COV()
	{
		Matrix ss_evectors = evd_cov.getV();
		Matrix D = new Matrix(ROWS, ROWS);
		for (int i = 0; i < ROWS; i++)
			{
				double eval = eigenvalues_prime_COV.get(ROWS - i - 1);
				D.set(i, i, eval);
			}
		inv_cov = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());

		top_evectors_COV = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_COV.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_COV = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix(top_evectors_COV, path, 15, 12);

		path = file_name_head + "_inverse_covariance_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(inv_cov, path);

		if (doReduce)
			{
				r_inv_cov = pca.get_Reduced_DYN_Matrix(inv_cov);
				path = file_name_head + "_reduced_inverse_covariance_matrix.txt.bz2";
				if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(r_inv_cov, path);
			}
	}

	private void construct_PCA_Modes_COV()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors_COV, eigenvalues_COV);

		pca_modes_COV = gPCA.get_PCA_modes();
		weighted_pca_modes_COV = gPCA.get_Weighted_PCA_modes();
		square_pca_modes_COV = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes_COV = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes_COV = gPCA.get_PCA_mode_maxs();
		pca_mode_mins_COV = gPCA.get_PCA_mode_mins();

		collectivity_cov = Collectivity.get_Collectivity(square_pca_modes_COV);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_cov, path);
		COLLECTIVITY_Plot.createPlot(out, "Eigenvector_Collectivity_COV", collectivity_cov);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes_COV, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes_COV, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes_COV, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes_COV, path);
	}

	/* ***************************************** CORR METHODS ************************************************************** */

	private void get_eigenvalues_CORR()
	{
		double[] ss_evals = evd_corr.getRealEigenvalues();
		eigenvalues_CORR = new ArrayList<>();
		trace_CORR = 0;
		for (double k : ss_evals)
			{
				if (k < 0) k = 0;
				eigenvalues_CORR.add(k);
				trace_CORR += k;
			}
		Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
		if (verbose) List_IO.write_Double_List(eigenvalues_CORR, file_name_head + "_eigenvalues_CORR.txt", 12);
	}

	private void write_top_evals_CORR()
	{
		Matrix screePlot = new Matrix(number_of_modes, 2);
		top_eigenvalues_CORR = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes; i++)
			{
				double val = eigenvalues_CORR.get(i);
				double normed_val = (val / trace_CORR);
				cumulative_variance += normed_val;
				top_eigenvalues_CORR.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues_CORR, trace_CORR);
		SCREE_Plot.createChart(out, "Scree_Plot_CORR", screePlot);
	}

	private void get_top_evects_and_reverse_CORR()
	{
		Matrix ss_evectors = evd_corr.getV();
		top_evectors_CORR = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_CORR.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_CORR = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix(top_evectors_CORR, path, 15, 12);
	}

	private void construct_PCA_Modes_CORR()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors_CORR, eigenvalues_CORR);

		pca_modes_CORR = gPCA.get_PCA_modes();
		weighted_pca_modes_CORR = gPCA.get_Weighted_PCA_modes();
		square_pca_modes_CORR = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes_CORR = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes_CORR = gPCA.get_PCA_mode_maxs();
		pca_mode_mins_CORR = gPCA.get_PCA_mode_mins();

		collectivity_corr = Collectivity.get_Collectivity(square_pca_modes_CORR);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_corr, path);
		COLLECTIVITY_Plot.createPlot(out, "Eigenvector_Collectivity_CORR", collectivity_corr);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes_CORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes_CORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes_CORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes_CORR, path);
	}

	/* ***************************************** SPARSE CORR METHODS ************************************************************** */

	private void get_eigenvalues_CORR_SPARSE()
	{
		double[] ss_evals = evd_corr_sparse.getRealEigenvalues();
		eigenvalues_CORR_SPARSE = new ArrayList<>();
		trace_CORR_SPARSE = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues_CORR_SPARSE.add(k);
				trace_CORR_SPARSE += k;
			}
		Collections.sort(eigenvalues_CORR_SPARSE, Collections.reverseOrder());
		if (verbose) List_IO.write_Double_List(eigenvalues_CORR_SPARSE, file_name_head + "_eigenvalues_CORR_SPARSE.txt", 12);
	}

	private void write_top_evals_CORR_SPARSE()
	{
		Matrix screePlot = new Matrix(number_of_modes, 2);
		top_eigenvalues_CORR_SPARSE = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes; i++)
			{
				double val = eigenvalues_CORR_SPARSE.get(i);
				double normed_val = (val / trace_CORR_SPARSE);
				cumulative_variance += normed_val;
				top_eigenvalues_CORR_SPARSE.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_CORR_SPARSE.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues_CORR_SPARSE, trace_CORR_SPARSE);
		SCREE_Plot.createChart(file_name_head + "_", "Scree_Plot_CORR_SPARSE", screePlot);
	}

	private void get_top_evects_and_reverse_CORR_SPARSE()
	{
		Matrix ss_evectors = evd_corr_sparse.getV();
		top_evectors_CORR_SPARSE = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_CORR_SPARSE.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_CORR_SPARSE = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_CORR_SPARSE, path, 15, 12);
	}

	private void construct_PCA_Modes_CORR_SPARSE()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors_CORR_SPARSE, eigenvalues_CORR_SPARSE);

		pca_modes_CORR_SPARSE = gPCA.get_PCA_modes();
		weighted_pca_modes_CORR_SPARSE = gPCA.get_Weighted_PCA_modes();
		square_pca_modes_CORR_SPARSE = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes_CORR_SPARSE = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes_CORR_SPARSE = gPCA.get_PCA_mode_maxs();
		pca_mode_mins_CORR_SPARSE = gPCA.get_PCA_mode_mins();

		collectivity_corr_SPARSE = Collectivity.get_Collectivity(square_pca_modes_CORR_SPARSE);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_corr, path);
		COLLECTIVITY_Plot.createPlot(file_name_head + "_", "Eigenvector_Collectivity_CORR_SPARSE", collectivity_corr_SPARSE);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes_CORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes_CORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes_CORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes_CORR_SPARSE, path);
	}

	/* ********************************************* PCORR METHODS ************************************************************** */

	private void get_eigenvalues_PCORR()
	{
		double[] ss_evals = evd_pcorr.getRealEigenvalues();
		eigenvalues_PCORR = new ArrayList<>();
		trace_PCORR = 0;
		for (double k : ss_evals)
			{
				if (k < 0) k = 0;
				eigenvalues_PCORR.add(k);
				trace_PCORR += k;
			}
		Collections.sort(eigenvalues_PCORR, Collections.reverseOrder());
		if (verbose) List_IO.write_Double_List(eigenvalues_PCORR, file_name_head + "_eigenvalues_PCORR.txt", 12);
	}

	private void write_top_evals_PCORR()
	{
		Matrix screePlot = new Matrix(number_of_modes, 2);
		top_eigenvalues_PCORR = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes; i++)
			{
				double val = eigenvalues_PCORR.get(i);
				double normed_val = (val / trace_PCORR);
				cumulative_variance += normed_val;
				top_eigenvalues_PCORR.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues_PCORR, trace_PCORR);
		SCREE_Plot.createChart(out, "Scree_Plot_PCORR", screePlot);
	}

	private void get_top_evects_and_reverse_PCORR()
	{
		Matrix ss_evectors = evd_pcorr.getV();
		top_evectors_PCORR = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);

		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_PCORR.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_PCORR = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix(top_evectors_PCORR, path, 15, 12);
	}

	private void construct_PCA_Modes_PCORR()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors_PCORR, eigenvalues_PCORR);

		pca_modes_PCORR = gPCA.get_PCA_modes();
		weighted_pca_modes_PCORR = gPCA.get_Weighted_PCA_modes();
		square_pca_modes_PCORR = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes_PCORR = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes_PCORR = gPCA.get_PCA_mode_maxs();
		pca_mode_mins_PCORR = gPCA.get_PCA_mode_mins();

		collectivity_pcorr = Collectivity.get_Collectivity(square_pca_modes_PCORR);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_pcorr, path);
		COLLECTIVITY_Plot.createPlot(out, "Eigenvector_Collectivity_PCORR", collectivity_pcorr);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes_PCORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes_PCORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes_PCORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes_PCORR, path);
	}

	/* ********************************************* SPARSE PCORR METHODS ************************************************************** */

	private void get_eigenvalues_PCORR_SPARSE()
	{
		double[] ss_evals = evd_pcorr_sparse.getRealEigenvalues();

		eigenvalues_PCORR_SPARSE = new ArrayList<>();

		trace_PCORR_SPARSE = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues_PCORR_SPARSE.add(k);
				trace_PCORR_SPARSE += k;
			}
		Collections.sort(eigenvalues_PCORR_SPARSE, Collections.reverseOrder());
		if (verbose) List_IO.write_Double_List(eigenvalues_PCORR_SPARSE, file_name_head + "_eigenvalues_PCORR_SPARSE.txt", 12);
	}

	private void write_top_evals_PCORR_SPARSE()
	{
		Matrix screePlot = new Matrix(number_of_modes, 2);
		top_eigenvalues_PCORR_SPARSE = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes; i++)
			{
				double val = eigenvalues_PCORR_SPARSE.get(i);
				double normed_val = (val / trace_PCORR_SPARSE);
				cumulative_variance += normed_val;
				top_eigenvalues_PCORR_SPARSE.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes + "_eigenvalues_PCORR_SPARSE.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues_PCORR_SPARSE, trace_PCORR_SPARSE);
		SCREE_Plot.createChart(file_name_head + "_", "Scree_Plot_PCORR_SPARSE", screePlot);
	}

	private void get_top_evects_and_reverse_PCORR_SPARSE()
	{
		Matrix ss_evectors = evd_pcorr_sparse.getV();
		top_evectors_PCORR_SPARSE = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);

		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_PCORR_SPARSE.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_PCORR_SPARSE = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_PCORR_SPARSE, path, 15, 12);
	}

	private void construct_PCA_Modes_PCORR_SPARSE()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors_PCORR_SPARSE, eigenvalues_PCORR_SPARSE);

		pca_modes_PCORR_SPARSE = gPCA.get_PCA_modes();
		weighted_pca_modes_PCORR_SPARSE = gPCA.get_Weighted_PCA_modes();
		square_pca_modes_PCORR_SPARSE = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes_PCORR_SPARSE = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes_PCORR_SPARSE = gPCA.get_PCA_mode_maxs();
		pca_mode_mins_PCORR_SPARSE = gPCA.get_PCA_mode_mins();

		collectivity_pcorr_SPARSE = Collectivity.get_Collectivity(square_pca_modes_PCORR_SPARSE);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_pcorr_SPARSE, path);
		COLLECTIVITY_Plot.createPlot(file_name_head + "_", "Eigenvector_Collectivity_PCORR_SPARSE", collectivity_pcorr_SPARSE);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes_PCORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes_PCORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes_PCORR_SPARSE, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes_PCORR_SPARSE, path);
	}

	/* ******************************************* SETTERS ***************************************************** */

	private void create_Directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
		file_name_head = out + "ss_" + number_of_residues + "_" + number_of_atoms;
	}

	public void set_Out_dir(String dir)
	{
		this.out_dir = dir;
		create_Directory(out_dir);
		file_name_head = out + "ss_" + number_of_residues + "_" + number_of_atoms;
	}

	public void setNumber_of_residues(int number_of_residues)
	{
		this.number_of_residues = number_of_residues;
	}

	/* ******************************************* GETTERS ***************************************************** */


	/**
	 * Returns the rank of the covariance matrix
	 *
	 * @return
	 */
	public int get_rank_COV()
	{
		return rank;
	}

	// ----------------------------------------------------------------------------------------------------------- //

	/**
	 * Returns the determinant of the covariance matrix
	 *
	 * @return
	 */
	public double get_det_COV()
	{
		return det;
	}

	/**
	 * Returns the condition number of the covariance matrix
	 *
	 * @return
	 */
	public double get_cond_COV()
	{
		return cond;
	}

	/**
	 * Returns the shrinkage intensity
	 *
	 * @return
	 */
	public double getShrinkage()
	{
		return shrinkage;
	}

	/**
	 * Returns the Trace of the covariance matrix
	 *
	 * @return
	 */
	public double get_trace_COV()
	{
		return trace_COV;
	}

	/**
	 * Returns the Trace of the correlation matrix
	 *
	 * @return
	 */
	public double get_trace_CORR()
	{
		return trace_CORR;
	}

	/**
	 * Returns the Trace of the partial correlation matrix
	 *
	 * @return
	 */
	public double get_trace_PCORR()
	{
		return trace_PCORR;
	}

	// ----------------------------------------------------------------------------------------------------------- //

	/**
	 * Returns the Eigenvalues of the covariance matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_COV()
	{
		return eigenvalues_COV;
	}

	/**
	 * Returns the Eigenvalues of the correlation matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_CORR()
	{
		return eigenvalues_CORR;
	}

	/**
	 * Returns the Eigenvalues of the partial correlation matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_PCORR()
	{
		return eigenvalues_PCORR;
	}

	/**
	 * Returns the Top eigenvalues from the COV analysis
	 *
	 * @return
	 */
	public List<Double> getTop_eigenvalues_COV()
	{

		return top_eigenvalues_COV;
	}

	/**
	 * Returns the Top eigenvalues from the CORR analysis
	 *
	 * @return
	 */
	public List<Double> getTop_eigenvalues_CORR()
	{

		return top_eigenvalues_CORR;
	}

	/**
	 * Returns the Top eigenvalues from the PCORR analysis
	 *
	 * @return
	 */
	public List<Double> getTop_eigenvalues_PCORR()
	{

		return top_eigenvalues_PCORR;
	}

	// ----------------------------------------------------------------------------------------------------------- //

	/**
	 * Returns the array of the PCA mode minimums from the COV analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_mins_COV()
	{
		return pca_mode_mins_COV;
	}

	/**
	 * Returns the array of the PCA mode maximums from the COV analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_maxes_COV()
	{
		return pca_mode_maxes_COV;
	}

	/**
	 * Returns the array of the PCA mode minimums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_mins_CORR()
	{
		return pca_mode_mins_CORR;
	}

	/**
	 * Returns the array of the PCA mode maximums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_maxes_CORR()
	{

		return pca_mode_maxes_CORR;
	}

	/**
	 *
	 * /** Returns the array of the PCA mode minimums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_mins_PCORR()
	{

		return pca_mode_mins_PCORR;
	}

	/**
	 * Returns the array of the PCA mode maximums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_maxes_PCORR()
	{

		return pca_mode_maxes_PCORR;
	}


	/**
	 * Returns the array of the PCA mode minimums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_mins_CORR_SPARSE()
	{
		return pca_mode_mins_CORR_SPARSE;
	}

	/**
	 * Returns the array of the PCA mode maximums from the CORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_maxes_CORR_SPARSE()
	{

		return pca_mode_maxes_CORR_SPARSE;
	}

	/**
	 *
	 * /** Returns the array of the PCA mode minimums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_mins_PCORR_SPARSE()
	{

		return pca_mode_mins_PCORR_SPARSE;
	}

	/**
	 * Returns the array of the PCA mode maximums from the PCORR analysis
	 *
	 * @return
	 */
	public double[] get_pca_mode_maxes_PCORR_SPARSE()
	{

		return pca_mode_maxes_PCORR_SPARSE;
	}

	// ----------------------------------------------------------------------------------------------------------- //

	public Matrix getCentered_input_data()
	{
		return centered_input_data;
	}

	/**
	 * Returns the top eigenvectors from the COV analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_COV()
	{

		return top_evectors_COV;
	}

	/**
	 * Returns the Top eigenvectors from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_CORR()
	{

		return top_evectors_CORR;
	}

	/**
	 * Returns the Top eigenvectors from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getTop_evectors_PCORR()
	{
		return top_evectors_PCORR;
	}

	/**
	 * Returns the Square PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_COV()
	{

		return square_pca_modes_COV;
	}

	/**
	 * Returns the Weighted Square PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_COV()
	{

		return weighted_square_pca_modes_COV;
	}

	/**
	 * Returns the Weighted PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_COV()
	{

		return weighted_pca_modes_COV;
	}

	/**
	 * Returns the Square PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_CORR()
	{

		return square_pca_modes_CORR;
	}

	/**
	 * Returns the Weighted Square PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_CORR()
	{

		return weighted_square_pca_modes_CORR;
	}

	/**
	 * Returns the Weighted PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_CORR()
	{

		return weighted_pca_modes_CORR;
	}

	/**
	 * Returns the PCA modes from the COV analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_COV()
	{

		return pca_modes_COV;
	}

	/**
	 * Returns the PCA modes from the CORR analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_CORR()
	{

		return pca_modes_CORR;
	}

	/**
	 * Returns the Square PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getSquare_pca_modes_PCORR()
	{

		return square_pca_modes_PCORR;
	}

	/**
	 * Returns the Weighted Square PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_square_pca_modes_PCORR()
	{

		return weighted_square_pca_modes_PCORR;
	}

	/**
	 * Returns the Weighted PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getWeighted_pca_modes_PCORR()
	{

		return weighted_pca_modes_PCORR;
	}

	/**
	 * Returns the PCA modes from the PCORR analysis
	 *
	 * @return
	 */
	public Matrix getPca_modes_PCORR()
	{

		return pca_modes_PCORR;
	}

	public double getKMO()
	{
		return KMO;
	}

	public List<Double> getEigenvalues_CORR_SPARSE()
	{
		return eigenvalues_CORR_SPARSE;
	}

	public List<Double> getTop_eigenvalues_CORR_SPARSE()
	{
		return top_eigenvalues_CORR_SPARSE;
	}

	public List<Double> getEigenvalues_PCORR_SPARSE()
	{
		return eigenvalues_PCORR_SPARSE;
	}

	public List<Double> getTop_eigenvalues_PCORR_SPARSE()
	{
		return top_eigenvalues_PCORR_SPARSE;
	}

	public double[] getPca_mode_mins_COV()
	{
		return pca_mode_mins_COV;
	}

	public double[] getPca_mode_maxes_COV()
	{
		return pca_mode_maxes_COV;
	}

	public double[] getPca_mode_mins_CORR()
	{
		return pca_mode_mins_CORR;
	}

	public double[] getPca_mode_maxes_CORR()
	{
		return pca_mode_maxes_CORR;
	}

	public double[] getPca_mode_mins_PCORR()
	{
		return pca_mode_mins_PCORR;
	}

	public double[] getPca_mode_maxes_PCORR()
	{
		return pca_mode_maxes_PCORR;
	}

	public double[] getPca_mode_mins_CORR_SPARSE()
	{
		return pca_mode_mins_CORR_SPARSE;
	}

	public double[] getPca_mode_maxes_CORR_SPARSE()
	{
		return pca_mode_maxes_CORR_SPARSE;
	}

	public double[] getPca_mode_mins_PCORR_SPARSE()
	{
		return pca_mode_mins_PCORR_SPARSE;
	}

	public double[] getPca_mode_maxes_PCORR_SPARSE()
	{
		return pca_mode_maxes_PCORR_SPARSE;
	}

	public Matrix getTop_evectors_CORR_SPARSE()
	{
		return top_evectors_CORR_SPARSE;
	}

	public Matrix getSquare_pca_modes_CORR_SPARSE()
	{
		return square_pca_modes_CORR_SPARSE;
	}

	public Matrix getWeighted_pca_modes_CORR_SPARSE()
	{
		return weighted_pca_modes_CORR_SPARSE;
	}

	public Matrix getPca_modes_CORR_SPARSE()
	{
		return pca_modes_CORR_SPARSE;
	}

	public Matrix getTop_evectors_PCORR_SPARSE()
	{
		return top_evectors_PCORR_SPARSE;
	}

	public Matrix getSquare_pca_modes_PCORR_SPARSE()
	{
		return square_pca_modes_PCORR_SPARSE;
	}

	public Matrix getWeighted_pca_modes_PCORR_SPARSE()
	{
		return weighted_pca_modes_PCORR_SPARSE;
	}

	public Matrix getPca_modes_PCORR_SPARSE()
	{
		return pca_modes_PCORR_SPARSE;
	}

	public Matrix getWeighted_square_pca_modes_CORR_SPARSE()
	{
		return weighted_square_pca_modes_CORR_SPARSE;
	}

	public Matrix getWeighted_square_pca_modes_PCORR_SPARSE()
	{
		return weighted_square_pca_modes_PCORR_SPARSE;
	}
}
