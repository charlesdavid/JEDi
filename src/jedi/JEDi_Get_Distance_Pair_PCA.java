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
 * JED class JED_Get_Distance_Pair_PCA: Gets the COV and CORR PCA for the Residue Distance Pairs subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Distance_Pair_PCA
{
	boolean success, exist;
	final boolean doCORR, doPCORR, doSPARSIFY, verbose;
	String out_dir, out_dir_COV, out_dir_CORR, out_dir_CORR_SPARSE, out_dir_PCORR, out_dir_PCORR_SPARSE, file_name_head, path, model_PCA;
	final String directory, description;
	final int ROWS_DP, COLS, number_of_modes, number_of_pairs;
	int rank;
	double shrinkage, trace_COV, trace_CORR, trace_PCORR, trace_CORR_SPARSE, trace_PCORR_SPARSE, cond, det;
	final double st_threshold, corr_threshold, pcorr_threshold, FLOOR, NOISE_LEVEL;
	List<Double> eigenvalues_COV, eigenvalues_prime_COV, top_eigenvalues_COV, eigenvalues_CORR, top_eigenvalues_CORR, eigenvalues_PCORR, top_eigenvalues_PCORR,
			eigenvalues_CORR_SPARSE, top_eigenvalues_CORR_SPARSE, eigenvalues_PCORR_SPARSE, top_eigenvalues_PCORR_SPARSE;
	Matrix cov_dist, cov_dist_ST, corr_dist, pcorr_dist, corr_dist_sparse, pcorr_dist_sparse, inv_cov;
	Matrix top_evectors_dist_COV, top_evectors_dist_CORR, top_evectors_dist_PCORR, top_evectors_dist_CORR_SPARSE, top_evectors_dist_PCORR_SPARSE;
	Matrix square_PCA_modes_COV, square_PCA_modes_CORR, square_PCA_modes_PCORR, square_PCA_modes_CORR_SPARSE, square_PCA_modes_PCORR_SPARSE, weighted_Square_PCA_modes_COV,
			weighted_Square_PCA_modes_CORR, weighted_Square_PCA_modes_CORR_SPARSE, weighted_Square_PCA_modes_PCORR, weighted_Square_PCA_modes_PCORR_SPARSE;
	Matrix collectivity_cov, collectivity_corr, collectivity_pcorr, collectivity_corr_SPARSE, collectivity_pcorr_SPARSE;
	final Matrix distances;
	final NumberFormat nf, nf6, nf12;
	final RoundingMode rm;
	EigenvalueDecomposition evd_cov, evd_corr, evd_pcorr, evd_corr_SPARSE, evd_pcorr_SPARSE;
	final PCA pca;

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
	JEDi_Get_Distance_Pair_PCA(Matrix dist, int pairs)
	{
		this.distances = dist;
		this.number_of_pairs = pairs;

		this.directory = Input_Parameters.DIRECTORY;
		this.description = Input_Parameters.DESCRIPTION;
		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doSPARSIFY = Input_Parameters.doSPARSIFY;
		this.st_threshold = Input_Parameters.THRESHOLD_COV;
		this.corr_threshold = Input_Parameters.THRESHOLD_CORR;
		this.pcorr_threshold = Input_Parameters.THRESHOLD_PCORR;
		this.FLOOR = Input_Parameters.FLOOR;
		this.NOISE_LEVEL = Input_Parameters.NOISE_LEVEL;
		this.number_of_modes = Input_Parameters.MODES_DISTANCE_PAIRS;
		this.verbose = Input_Parameters.verbose;

		this.ROWS_DP = number_of_pairs;
		this.COLS = number_of_modes;

		this.pca = new PCA(distances);

		this.nf = NumberFormat.getInstance();
		this.rm = RoundingMode.HALF_UP;
		this.nf.setRoundingMode(rm);
		this.nf.setMaximumFractionDigits(3);
		this.nf.setMinimumFractionDigits(3);

		this.nf6 = NumberFormat.getInstance();
		this.nf6.setRoundingMode(rm);
		this.nf6.setMaximumFractionDigits(6);
		this.nf6.setMinimumFractionDigits(6);

		this.nf12 = NumberFormat.getInstance();
		this.nf12.setRoundingMode(rm);
		this.nf12.setMaximumFractionDigits(12);
		this.nf12.setMinimumFractionDigits(12);
	}

	/* *********************************** DRIVER METHODS ********************************************************************************** */

	/**
	 * Performs the Distance Pair PCA using COV, CORR, and PCORR (Q, R, and P) models.
	 */
	public void get_Distance_Pair_PCA()
	{
		Do_Cov();
		if (doCORR)
			{
				Do_Corr();
				if (doSPARSIFY)
					{
						Do_CORR_SPARSE_PCA();
					}
			}
		if (doPCORR)
			{
				Do_PCorr();
				if (doSPARSIFY)
					{
						Do_PCORR_SPARSE_PCA();
					}
			}
	}

	/**
	 * Performs the COV analysis.
	 */
	private void Do_Cov()
	{
		model_PCA = "COV";
		cov_dist = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		shrinkage = pca.get_lambda_cov();

		if (st_threshold > 0)
			{
				cov_dist_ST = pca.Threshold_COV_Matrix(cov_dist);
				cov_dist = cov_dist_ST;
			}

		out_dir_COV = out_dir + model_PCA + File.separatorChar;
		create_directory(out_dir_COV);
		path = file_name_head + "_covariance_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(cov_dist, path);

		HeatMap hm_cov_dist = new HeatMap(cov_dist, out_dir_COV, "Distance_Pair_Covariance_Matrix", "Atom 1", "Atom 2");
		hm_cov_dist.createMatrixHeatmap();

		evd_cov = pca.get_eigenvalue_decomposition(cov_dist);

		get_eigenvalues_COV();
		write_top_evals_COV();
		get_top_evects_and_reverse_COV();
		construct_PCA_Modes_COV();

		corr_dist = pca.get_Correlation_Matrix_from_Covariance_Matrix(cov_dist);
		pcorr_dist = pca.get_Partial_Correlation_Matrix(inv_cov);
	}

	/**
	 * Performs the CORR analysis.
	 */
	private void Do_Corr()
	{
		model_PCA = "CORR";
		out_dir_CORR = out_dir + model_PCA + File.separatorChar;
		create_directory(out_dir_CORR);
		path = file_name_head + "_correlation_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(corr_dist, path);

		HeatMap hm_corr_dist = new HeatMap(corr_dist, out_dir_CORR, "Distance_Pair_Correlation_Matrix", "Atom 1", "Atom 2");
		hm_corr_dist.createMatrixHeatmap();

		evd_corr = pca.get_eigenvalue_decomposition(corr_dist);

		get_eigenvalues_CORR();
		write_top_evals_CORR();
		get_top_evects_and_reverse_CORR();
		construct_PCA_Modes_CORR();
	}

	/**
	 * Performs the SPARSE CORR analysis.
	 */
	private void Do_CORR_SPARSE_PCA()
	{
		out_dir_CORR_SPARSE = out_dir + model_PCA + File.separatorChar + "sparse" + File.separatorChar;
		create_directory(out_dir_CORR_SPARSE);

		corr_dist_sparse = pca.Sparsify_Matrix(corr_dist, Input_Parameters.THRESHOLD_CORR);
		path = file_name_head + "_sparsified_correlation_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(corr_dist, path);

		HeatMap hm_corrS_dist = new HeatMap(corr_dist, out_dir_CORR_SPARSE, "Sparsified_Distance_Pair_Correlation_Matrix_" + Input_Parameters.THRESHOLD_CORR, "Atom Pair 1",
				"Atom Pair 2");
		hm_corrS_dist.createMatrixHeatmap();

		evd_corr_SPARSE = pca.get_eigenvalue_decomposition(corr_dist_sparse);

		get_eigenvalues_CORR_SPARSE();
		write_top_evals_CORR_SPARSE();
		get_top_evects_and_reverse_CORR_SPARSE();
		construct_PCA_Modes_CORR_SPARSE();
	}

	/**
	 * Performs the PCORR analysis.
	 */
	private void Do_PCorr()
	{
		model_PCA = "PCORR";
		out_dir_PCORR = out_dir + model_PCA + File.separatorChar;
		create_directory(out_dir_PCORR);

		path = file_name_head + "_partial_correlation_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pcorr_dist, path);

		HeatMap hm_pcorr_dist = new HeatMap(pcorr_dist, out_dir_PCORR, "Distance_Pair_Partial_Correlation_Matrix", "Atom 1", "Atom 2");
		hm_pcorr_dist.createMatrixHeatmap();

		evd_pcorr = pca.get_eigenvalue_decomposition(pcorr_dist);

		get_eigenvalues_PCORR();
		write_top_evals_PCORR();
		get_top_evects_and_reverse_PCORR();
		construct_PCA_Modes_PCORR();
	}

	/**
	 * Performs the SPARSE PCORR analysis.
	 */
	private void Do_PCORR_SPARSE_PCA()
	{
		out_dir_PCORR_SPARSE = out_dir + model_PCA + File.separatorChar + "sparse" + File.separatorChar;
		create_directory(out_dir_PCORR_SPARSE);

		pcorr_dist_sparse = pca.Sparsify_Matrix(pcorr_dist, Input_Parameters.THRESHOLD_PCORR);
		path = file_name_head + "_sparsified_partial_correlation_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(corr_dist_sparse, path);

		HeatMap hm_corrS_dist = new HeatMap(pcorr_dist_sparse, out_dir_PCORR_SPARSE, "Sparsified_Partial_Correlation_Matrix_" + Input_Parameters.THRESHOLD_PCORR, "Atom Pair 1",
				"Atom Pair 2");
		hm_corrS_dist.createMatrixHeatmap();


		Matrix diff = corr_dist_sparse.minus(pcorr_dist_sparse);
		Matrix diff_abs = new Matrix(corr_dist_sparse.getRowDimension(), corr_dist_sparse.getRowDimension());
		Matrix adjacency = new Matrix(corr_dist_sparse.getRowDimension(), corr_dist_sparse.getRowDimension());

		for (int i = 0; i < corr_dist_sparse.getRowDimension(); i++)
			{
				for (int j = 0; j < i; j++)
					{
						double r = corr_dist_sparse.get(i, j);
						double r_abs = Math.abs(r);
						double p = pcorr_dist_sparse.get(i, j);
						double p_abs = Math.abs(p);

						diff_abs.set(i, j, (r_abs - p_abs));
						diff_abs.set(j, i, (r_abs - p_abs));
					}
			}

		HeatMap hm = new HeatMap(diff, out_dir_PCORR_SPARSE, "Difference of Correlation and  Partial Correlation Matrices", "Variable 1", "Variable 2");
		// hm.createMatrixHeatmap();

		for (int i = 0; i < corr_dist_sparse.getRowDimension(); i++)
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

		hm = new HeatMap(adjacency, out_dir_PCORR_SPARSE, "Adjacency Matrix derived from difference of Correlation and Partial Correlation Matrices", "Variable 1", "Variable 2");
		hm.createAdjacencyHeatmap();


		evd_pcorr_SPARSE = pca.get_eigenvalue_decomposition(pcorr_dist_sparse);

		get_eigenvalues_PCORR_SPARSE();
		write_top_evals_PCORR_SPARSE();
		get_top_evects_and_reverse_PCORR_SPARSE();
		construct_PCA_Modes_PCORR_SPARSE();
	}

	/* ************************************* COV METHODS ************************************************************************* ************ */

	/**
	 * Gets the eigenvalues from the COV analysis.
	 */
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
				if (k > 0) rankCount++;
				if (k < FLOOR) k = FLOOR;
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

	/**
	 * Gets the top eigenvalues from the COV analysis.
	 */
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
		SCREE_Plot.createChart(out_dir_COV, "Scree_Plot_COV", screePlot);
	}

	/**
	 * Gets the top eigenvectors from the COV analysis.
	 */
	private void get_top_evects_and_reverse_COV()
	{
		Matrix ss_evectors = evd_cov.getV();
		Matrix D = new Matrix(ROWS_DP, ROWS_DP);
		for (int i = 0; i < ROWS_DP; i++)
			{
				double eval = eigenvalues_prime_COV.get(ROWS_DP - i - 1);
				D.set(i, i, eval);
			}
		inv_cov = ss_evectors.times(D.inverse()).times(ss_evectors.transpose());

		top_evectors_dist_COV = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);
		Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_dist_COV.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
			}
		top_evectors_dist_COV = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_COV.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_dist_COV, path, 15, 12);

		path = file_name_head + "_inverse_covariance_matrix.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(inv_cov, path);
	}

	/**
	 * Computes the Distance-Pair PCA modes from the COV analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
	 */
	private void construct_PCA_Modes_COV()
	{
		JEDi_Get_PCA_Modes gpca = new JEDi_Get_PCA_Modes(top_evectors_dist_COV, eigenvalues_COV);

		Matrix SS_pca_modes = gpca.get_PCA_modes();
		square_PCA_modes_COV = gpca.get_Square_PCA_modes();
		Matrix SS_weighted_pca_modes = gpca.get_Weighted_PCA_modes();
		weighted_Square_PCA_modes_COV = gpca.get_Weighted_Square_PCA_modes();

		collectivity_cov = Collectivity.get_Collectivity(square_PCA_modes_COV);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_cov, path);
		COLLECTIVITY_Plot.createPlot(out_dir_COV, "Eigenvector_Collectivity_COV", collectivity_cov);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_PCA_modes_COV, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_COV.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_Square_PCA_modes_COV, path);
	}

	/* ************************************** CORR METHODS ************************************************************************* ********** */

	/**
	 * Gets the eigenvalues from the CORR analysis.
	 */
	private void get_eigenvalues_CORR()
	{
		double[] ss_evals = evd_corr.getRealEigenvalues();
		eigenvalues_CORR = new ArrayList<>();
		trace_CORR = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues_CORR.add(k);
				trace_CORR += k;
			}
		Collections.sort(eigenvalues_CORR, Collections.reverseOrder());
		path = file_name_head + "_eigenvalues_CORR.txt";
		List_IO.write_Double_List(eigenvalues_CORR, path, 12);
	}

	/**
	 * Gets the top eigenvalues from the CORR analysis.
	 */
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
		SCREE_Plot.createChart(out_dir_CORR, "Scree_Plot_CORR", screePlot);
	}

	/**
	 * Gets the top eigenvectors from the CORR analysis.
	 */
	private void get_top_evects_and_reverse_CORR()
	{
		Matrix ss_evectors = evd_corr.getV();
		top_evectors_dist_CORR = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);

		Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
			}
		top_evectors_dist_CORR = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_dist_CORR, path, 15, 12);
	}

	/**
	 * Computes the Distance-Pair PCA modes from the CORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
	 */
	private void construct_PCA_Modes_CORR()
	{
		JEDi_Get_PCA_Modes gpca = new JEDi_Get_PCA_Modes(top_evectors_dist_CORR, eigenvalues_CORR);

		Matrix SS_pca_modes = gpca.get_PCA_modes();
		square_PCA_modes_CORR = gpca.get_Square_PCA_modes();
		Matrix SS_weighted_pca_modes = gpca.get_Weighted_PCA_modes();
		weighted_Square_PCA_modes_CORR = gpca.get_Weighted_Square_PCA_modes();

		collectivity_corr = Collectivity.get_Collectivity(square_PCA_modes_CORR);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_CORR.txt.bz2";
		Matrix_IO.write_Matrix_adaptive_spacing(collectivity_corr, path);
		COLLECTIVITY_Plot.createPlot(out_dir_CORR, "Eigenvector_Collectivity_CORR", collectivity_corr);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_PCA_modes_CORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_Square_PCA_modes_CORR, path);
	}

	/* *********************************** PCORR METHODS ******************************************************************** */

	/**
	 * Gets the eigenvalues from the PCORR analysis.
	 */
	private void get_eigenvalues_PCORR()
	{
		double[] ss_evals = evd_pcorr.getRealEigenvalues();
		eigenvalues_PCORR = new ArrayList<>();
		trace_PCORR = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues_PCORR.add(k);
				trace_PCORR += k;
			}
		Collections.sort(eigenvalues_PCORR, Collections.reverseOrder());
		if (verbose) List_IO.write_Double_List(eigenvalues_PCORR, file_name_head + "_eigenvalues_PCORR.txt", 12);
	}

	/**
	 * Gets the top eigenvalues from the PCORR analysis.
	 */
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
		SCREE_Plot.createChart(out_dir_PCORR, "Scree_Plot_PCORR", screePlot);
	}

	/**
	 * Gets the top eigenvectors from the PCORR analysis.
	 */
	private void get_top_evects_and_reverse_PCORR()
	{
		Matrix ss_evectors = evd_pcorr.getV();
		top_evectors_dist_PCORR = ss_evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);

		Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_dist_PCORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
			}
		top_evectors_dist_PCORR = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_dist_PCORR, path, 15, 12);
	}

	/**
	 * Computes the Distance-Pair PCA modes from the PCORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
	 */
	private void construct_PCA_Modes_PCORR()
	{
		JEDi_Get_PCA_Modes gpca = new JEDi_Get_PCA_Modes(top_evectors_dist_PCORR, eigenvalues_PCORR);

		Matrix SS_pca_modes = gpca.get_PCA_modes();
		square_PCA_modes_PCORR = gpca.get_Square_PCA_modes();
		Matrix SS_weighted_pca_modes = gpca.get_Weighted_PCA_modes();
		weighted_Square_PCA_modes_PCORR = gpca.get_Weighted_Square_PCA_modes();

		collectivity_pcorr = Collectivity.get_Collectivity(square_PCA_modes_PCORR);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_pcorr, path);
		COLLECTIVITY_Plot.createPlot(out_dir_PCORR, "Eigenvector_Collectivity_PCORR", collectivity_pcorr);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(SS_weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_PCA_modes_PCORR, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_Square_PCA_modes_PCORR, path);
	}

	/* ************************************** SPARSE CORR METHODS ************************************************************************* ********** */

	/**
	 * Gets the eigenvalues from the SPARSE CORR analysis.
	 */
	private void get_eigenvalues_CORR_SPARSE()
	{
		double[] ss_evals = evd_corr_SPARSE.getRealEigenvalues();
		eigenvalues_CORR_SPARSE = new ArrayList<>();
		trace_CORR_SPARSE = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues_CORR_SPARSE.add(k);
				trace_CORR_SPARSE += k;
			}
		Collections.sort(eigenvalues_CORR_SPARSE, Collections.reverseOrder());
		path = file_name_head + "_eigenvalues_CORR_SPARSE.txt";
		if (verbose) List_IO.write_Double_List(eigenvalues_CORR_SPARSE, path, 12);
	}

	/**
	 * Gets the top eigenvalues from the SPARSE CORR analysis.
	 */
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
		SCREE_Plot.createChart(out_dir_CORR_SPARSE, "Scree_Plot_CORR_SPARSE", screePlot);
	}

	/**
	 * Gets the top eigenvectors from the SPARSE CORR analysis.
	 */
	private void get_top_evects_and_reverse_CORR_SPARSE()
	{
		Matrix evectors = evd_corr_SPARSE.getV();
		top_evectors_dist_CORR_SPARSE = evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);

		Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_dist_CORR.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
			}
		top_evectors_dist_CORR_SPARSE = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_dist_CORR_SPARSE, path, 15, 12);
	}

	/**
	 * Computes the Distance-Pair PCA modes from the SPARSE CORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
	 */
	private void construct_PCA_Modes_CORR_SPARSE()
	{
		JEDi_Get_PCA_Modes gpca = new JEDi_Get_PCA_Modes(top_evectors_dist_CORR_SPARSE, eigenvalues_CORR_SPARSE);

		Matrix pca_modes = gpca.get_PCA_modes();
		square_PCA_modes_CORR_SPARSE = gpca.get_Square_PCA_modes();
		Matrix weighted_pca_modes = gpca.get_Weighted_PCA_modes();
		weighted_Square_PCA_modes_CORR_SPARSE = gpca.get_Weighted_Square_PCA_modes();

		collectivity_corr_SPARSE = Collectivity.get_Collectivity(square_PCA_modes_CORR);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_corr, path);
		COLLECTIVITY_Plot.createPlot(out_dir_CORR_SPARSE, "Eigenvector_Collectivity_CORR_SPARSE", collectivity_corr_SPARSE);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_CORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_PCA_modes_CORR_SPARSE, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_CORR.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_Square_PCA_modes_CORR_SPARSE, path);
	}

	/* *********************************** SPARSE PCORR METHODS ******************************************************************** */

	/**
	 * Gets the eigenvalues from the SPARSE PCORR analysis.
	 */
	private void get_eigenvalues_PCORR_SPARSE()
	{
		double[] ss_evals = evd_pcorr_SPARSE.getRealEigenvalues();
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

	/**
	 * Gets the top eigenvalues from the PCORR analysis.
	 */
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
		SCREE_Plot.createChart(out_dir_PCORR_SPARSE, "Scree_Plot_PCORR_SPARSE", screePlot);
	}

	/**
	 * Gets the top eigenvectors from the SPARSE PCORR analysis.
	 */
	private void get_top_evects_and_reverse_PCORR_SPARSE()
	{
		Matrix evectors = evd_pcorr_SPARSE.getV();
		top_evectors_dist_PCORR_SPARSE = evectors.getMatrix(0, ROWS_DP - 1, ROWS_DP - number_of_modes, ROWS_DP - 1);

		Matrix modes_reversed = new Matrix(ROWS_DP, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_dist_PCORR_SPARSE.getMatrix(0, ROWS_DP - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS_DP - 1, r, r, col);
			}
		top_evectors_dist_PCORR_SPARSE = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes + "_eigenvectors_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_evectors_dist_PCORR_SPARSE, path, 15, 12);
	}

	/**
	 * Computes the Distance-Pair PCA modes from the SPARSE PCORR analysis and writes them to files: PCA, weighted-PCA, squared, weighted-squared
	 */
	private void construct_PCA_Modes_PCORR_SPARSE()
	{
		JEDi_Get_PCA_Modes gpca = new JEDi_Get_PCA_Modes(top_evectors_dist_PCORR_SPARSE, eigenvalues_PCORR_SPARSE);

		Matrix pca_modes = gpca.get_PCA_modes();
		square_PCA_modes_PCORR_SPARSE = gpca.get_Square_PCA_modes();
		Matrix weighted_pca_modes = gpca.get_Weighted_PCA_modes();
		weighted_Square_PCA_modes_PCORR_SPARSE = gpca.get_Weighted_Square_PCA_modes();

		collectivity_pcorr_SPARSE = Collectivity.get_Collectivity(square_PCA_modes_PCORR_SPARSE);
		path = file_name_head + "_top_" + number_of_modes + "_Eigenvector_Collectivity_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity_pcorr_SPARSE, path);
		COLLECTIVITY_Plot.createPlot(out_dir_PCORR_SPARSE, "Eigenvector_Collectivity_PCORR_SPARSE", collectivity_pcorr);

		path = file_name_head + "_top_" + number_of_modes + "_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes + "_square_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_PCA_modes_PCORR_SPARSE, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_square_pca_modes_PCORR_SPARSE.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_Square_PCA_modes_PCORR_SPARSE, path);
	}

	/* ************************************** SETTERS ***************************************************************************************** */

	private void create_directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
		file_name_head = dir + "ss_" + number_of_pairs + "_Atom_Pairs";
	}


	public void set_Out_dir(String out)
	{
		this.out_dir = out;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
		file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";
	}

	/* ************************************** GETTERS ************************************************************************* **************** */

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
	 * @return the shrinkage
	 */
	public double getShrinkage()
	{
		return shrinkage;
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
		return cond;
	}

	/**
	 * @return the det_COV
	 */
	public double getDet_COV()
	{
		return det;
	}

	/**
	 * @return the rank_COV
	 */
	public int getRank_COV()
	{
		return rank;
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
	 * @return the distances
	 */
	public Matrix getDistances()
	{
		return distances;
	}

	/**
	 * @return the cOV_dist
	 */
	public Matrix getCOV_dist()
	{
		return cov_dist;
	}

	/**
	 * @return the cORR_dist
	 */
	public Matrix getCORR_dist()
	{
		return corr_dist;
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
		return inv_cov;
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

	public Matrix getWeightedSquarePCAmodesCOV()
	{
		return weighted_Square_PCA_modes_COV;
	}

	public Matrix getWeightedSquarePCAmodesCORR()
	{
		return weighted_Square_PCA_modes_CORR;
	}

	public Matrix getWeightedSquarePCAmodesPCORR()
	{
		return weighted_Square_PCA_modes_PCORR;
	}

	public double getTrace_CORR_SPARSE()
	{
		return trace_CORR_SPARSE;
	}

	public double getTrace_PCORR_SPARSE()
	{
		return trace_PCORR_SPARSE;
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

	public Matrix getTop_evectors_dist_CORR_SPARSE()
	{
		return top_evectors_dist_CORR_SPARSE;
	}

	public Matrix getTop_evectors_dist_PCORR_SPARSE()
	{
		return top_evectors_dist_PCORR_SPARSE;
	}

	public Matrix getSquare_PCA_modes_CORR_SPARSE()
	{
		return square_PCA_modes_CORR_SPARSE;
	}

	public Matrix getSquare_PCA_modes_PCORR_SPARSE()
	{
		return square_PCA_modes_PCORR_SPARSE;
	}

	public Matrix getWeighted_Square_PCA_modes_CORR_SPARSE()
	{
		return weighted_Square_PCA_modes_CORR_SPARSE;
	}

	public Matrix getWeighted_Square_PCA_modes_PCORR_SPARSE()
	{
		return weighted_Square_PCA_modes_PCORR_SPARSE;
	}
}
