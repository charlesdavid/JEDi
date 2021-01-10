package jedi;

import java.io.File;
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
 * JED class JED_Get_Generalized_Coordinate_PCA: Gets the PCA using COV for the Cartesian subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Get_Hierarchical_PCA
{
	boolean success, exist, verbose;
	double trace, cond, det, shrinkage;
	final double FLOOR, NOISE_LEVEL;
	final int number_of_modes_residue, number_of_modes_hierarchical, number_of_residues, ROWS, COLS;
	int number_of_atoms, rank;
	String out_dir, file_name_head, path;
	List<Double> eigenvalues, top_eigenvalues;
	double[] pca_mode_mins, pca_mode_maxes;
	Matrix centered_input_data, cov, top_evectors, square_pca_modes, weighted_square_pca_modes, weighted_pca_modes, pca_modes, data_Means, data_Sigmas;
	Matrix collectivity;
	final Matrix input_data;
	EigenvalueDecomposition evd;
	final PCA pca;

	/**
	 * Constructor to perform PCA on generalized residue coordinates
	 *
	 * @param data      The subset coordinates
	 * @param num_modes The number of PCA modes to process
	 * @param dir       The working directory
	 * @param des       The job description
	 */
	JEDi_Get_Hierarchical_PCA(Matrix data, int num_modes_res, int num_modes_hier)
	{
		this.input_data = data;
		this.number_of_modes_residue = num_modes_res;
		this.number_of_modes_hierarchical = num_modes_hier;
		this.ROWS = input_data.getRowDimension();
		this.COLS = number_of_modes_hierarchical;
		this.number_of_residues = (ROWS / number_of_modes_residue);

		this.FLOOR = Input_Parameters.FLOOR;
		this.NOISE_LEVEL = Input_Parameters.NOISE_LEVEL;
		this.verbose = Input_Parameters.verbose;

		pca = new PCA(input_data);
	}

	/* **************************************************** PUBLIC METHOD ****************************************************** */

	/**
	 * Performs analysis using the COV PCA model.
	 */
	public void get_Hierarchical_PCA()
	{
		Do_Cov_Model();
	}

	/* **************************************************** PRIVATE METHODS ****************************************************** */

	private void Do_Cov_Model()
	{
		cov = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		if (Input_Parameters.THRESHOLD_COV > 0)
			{
				Matrix cov_st = pca.Threshold_COV_Matrix(cov);
				cov = cov_st;
			}

		HeatMap hm_cov = new HeatMap(cov, out_dir, "Covariance_Matrix_of_Residue_PCs", "Variable 1", "Variable 2");
		hm_cov.createMatrixHeatmap();

		centered_input_data = pca.getCentered_data_Matrix();
		data_Means = pca.getData_means();
		data_Sigmas = pca.getData_sigmas();
		shrinkage = pca.get_lambda_cov();

		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(cov, file_name_head + "_covariance_matrix_of_residue_PCs.txt.bz2");

		evd = pca.get_eigenvalue_decomposition(cov);
		get_eigenvalues_COV();
		write_top_evals_COV();
		get_top_evects_and_reverse_COV();
		construct_PCA_Modes_COV();
	}

	/* **************************************** COV METHODS ************************************************************* */

	private void get_eigenvalues_COV()
	{
		double[] ss_evals = evd.getRealEigenvalues();
		eigenvalues = new ArrayList<Double>();
		trace = 0;
		det = 1;
		int rankCount = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				if (k > 0) rankCount++;
				eigenvalues.add(k);
				trace += k;
				det *= k;
			}
		rank = rankCount;
		Collections.sort(eigenvalues, Collections.reverseOrder());

		int size = eigenvalues.size();
		double first = eigenvalues.get(0);
		double last = eigenvalues.get(size - 1);
		cond = (first / last);

		if (verbose) List_IO.write_Double_List(eigenvalues, file_name_head + "_eigenvalues.txt", 12);
	}

	private void write_top_evals_COV()
	{
		Matrix screePlot = new Matrix(number_of_modes_hierarchical, 2);
		top_eigenvalues = new ArrayList<>();
		double cumulative_variance = 0;
		for (int i = 0; i < number_of_modes_hierarchical; i++)
			{
				double val = eigenvalues.get(i);
				double normed_val = (val / trace);
				cumulative_variance += normed_val;
				top_eigenvalues.add(val);
				screePlot.set(i, 0, val);
				screePlot.set(i, 1, cumulative_variance);
			}
		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_eigenvalues.txt";
		if (verbose) Write_Top_Eigenvalues.write_evals(path, top_eigenvalues, trace);
		SCREE_Plot.createChart(out_dir, "Scree_Plot_Hierarchical", screePlot);
	}

	private void get_top_evects_and_reverse_COV()
	{
		Matrix ss_evectors = evd.getV();
		top_evectors = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes_hierarchical, ROWS - 1);
		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors = modes_reversed;

		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_eigenvectors_of_residue_PCs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix(top_evectors, path, 15, 12);
	}

	private void construct_PCA_Modes_COV()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(top_evectors, eigenvalues);

		pca_modes = gPCA.get_PCA_modes();
		weighted_pca_modes = gPCA.get_Weighted_PCA_modes();
		square_pca_modes = gPCA.get_Square_PCA_modes();
		weighted_square_pca_modes = gPCA.get_Weighted_Square_PCA_modes();
		pca_mode_maxes = gPCA.get_PCA_mode_maxs();
		pca_mode_mins = gPCA.get_PCA_mode_mins();

		collectivity = Collectivity.get_Collectivity(square_pca_modes);
		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_Eigenvector_Collectivity.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity, path);
		COLLECTIVITY_Plot.createPlot(out_dir, "Eigenvector_Collectivity", collectivity);

		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_weighted_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_square_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(square_pca_modes, path);
		path = file_name_head + "_top_" + number_of_modes_hierarchical + "_weighted_square_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_square_pca_modes, path);
	}

	/* ******************************************* SETTERS ***************************************************** */

	public void set_Out_Dir(String dir) // be sure to set the number of atoms first!
	{
		this.out_dir = dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();

		this.file_name_head = out_dir + "ss_" + number_of_residues + "_" + number_of_atoms;
	}

	public void setNumber_of_atoms(int numberOfAtoms)
	{
		this.number_of_atoms = numberOfAtoms;
	}

	/* ******************************************* GETTERS ***************************************************** */

	public double get_shrinkage()
	{
		return shrinkage;
	}

	public double get_cond()
	{
		return cond;
	}

	public double get_trace()
	{
		return trace;
	}

	public List<Double> getEigenvalues()
	{

		return eigenvalues;
	}

	public double get_det()
	{
		return det;
	}

	public int get_rank()
	{
		return rank;
	}

	public double[] get_pca_mode_mins()
	{
		return pca_mode_mins;
	}

	public double[] get_pca_mode_maxes()
	{

		return pca_mode_maxes;
	}

	public Matrix getTop_evectors()
	{

		return top_evectors;
	}

	public Matrix getSquare_pca_modes()
	{
		return square_pca_modes;
	}

	public Matrix getWeighted_square_pca_modes()
	{
		return weighted_square_pca_modes;
	}

	public Matrix getWeighted_pca_modes()
	{
		return weighted_pca_modes;
	}

	public Matrix getPca_modes()
	{
		return pca_modes;
	}

	public List<Double> getTop_eigenvalues()
	{
		return top_eigenvalues;
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
