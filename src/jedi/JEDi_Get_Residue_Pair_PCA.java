package jedi;

import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import support.PCA;
import supportIO.Input_Parameters;

/**
 * JED class JED_Get_Cartesian_PCA: Gets the PCA using COV and CORR for the Cartesian subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Get_Residue_Pair_PCA
{
	boolean success, exist;
	double trace_COV, cond_COV, det_COV;
	final double FLOOR = Input_Parameters.FLOOR;
	int rank_COV;
	final int number_of_atoms, number_of_modes, ROWS, COLS;
	List<Double> eigenvalues_COV;

	Matrix centered_input_data, centroids, sigmas;
	Matrix cov;
	Matrix top_evectors_COV;

	final Matrix input_data;
	EigenvalueDecomposition evd;
	final NumberFormat nf, nf6, nf12;
	final RoundingMode rm;
	final PCA pca;

	/**
	 * Constructor to perform the Residue Pair PCA
	 *
	 * @param data           The Cartesian subset coordinates
	 * @param num_modes      The number of PCA modes to process
	 * @param dir            The working directory
	 * @param des            The job description
	 * @param res_pair_index The residue pair indices
	 */
	JEDi_Get_Residue_Pair_PCA(Matrix data, int num_modes)
	{
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

		this.input_data = data;
		this.number_of_modes = num_modes;

		this.ROWS = input_data.getRowDimension();
		this.COLS = number_of_modes;
		this.number_of_atoms = (ROWS / 3);

		this.pca = new PCA(input_data);
	}

	/* *************************************** PRIMARY METHODS ****************************************************** */

	/**
	 * Performs the COV, CORR, and PCORR PCA methods
	 */
	public void get_Residue_Pair_PCA()
	{
		Do_Cov_PCA();
	}

	private void Do_Cov_PCA()
	{
		cov = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();

		centroids = pca.getData_means();
		sigmas = pca.getData_sigmas();
		centered_input_data = pca.getCentered_data_Matrix();

		evd = pca.get_eigenvalue_decomposition(cov);
		get_eigenvalues_COV();
		get_top_evects_and_reverse_COV();
	}
	/* **************************************** COV METHODS ************************************************************* */

	private void get_eigenvalues_COV()
	{
		double[] ss_evals = evd.getRealEigenvalues();
		eigenvalues_COV = new ArrayList<Double>();
		trace_COV = 0;
		det_COV = 1;
		int rankCount = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				if (k > 0) rankCount++;
				eigenvalues_COV.add(k);
				trace_COV += k;
				det_COV *= k;
			}
		rank_COV = rankCount;
		Collections.sort(eigenvalues_COV, Collections.reverseOrder());

		int size = eigenvalues_COV.size();
		double first = eigenvalues_COV.get(0);
		double last = eigenvalues_COV.get(size - 1);
		cond_COV = (first / last);
	}

	private void get_top_evects_and_reverse_COV()
	{
		Matrix ss_evectors = evd.getV();
		Matrix D = new Matrix(ROWS, ROWS);
		for (int i = 0; i < ROWS; i++)
			{
				double eval = eigenvalues_COV.get(ROWS - i - 1);
				D.set(i, i, eval);
			}
		top_evectors_COV = ss_evectors.getMatrix(0, ROWS - 1, ROWS - number_of_modes, ROWS - 1);
		Matrix modes_reversed = new Matrix(ROWS, COLS);
		for (int r = 0; r < COLS; r++)
			{
				Matrix col = top_evectors_COV.getMatrix(0, ROWS - 1, COLS - 1 - r, COLS - 1 - r);
				modes_reversed.setMatrix(0, ROWS - 1, r, r, col);
			}
		top_evectors_COV = modes_reversed;
	}

	/* ******************************************* GETTERS ***************************************************** */

	/**
	 * Returns the condition number of the covariance matrix
	 *
	 * @return
	 */
	public double get_cond_COV()
	{
		return cond_COV;
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
	 * Returns the Eigenvalues of the covariance matrix
	 *
	 * @return
	 */
	public List<Double> getEigenvalues_COV()
	{
		return eigenvalues_COV;
	}

	/**
	 * Returns the determinant of the covariance matrix
	 *
	 * @return
	 */
	public double get_det_COV()
	{
		return det_COV;
	}

	/**
	 * Returns the rank of the covariance matrix
	 *
	 * @return
	 */
	public double get_rank_COV()
	{
		return rank_COV;
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

	public Matrix getCentered_input_data()
	{
		return centered_input_data;
	}
}
