package jedi;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import support.PCA;
import supportIO.Input_Parameters;

/**
 * JEDi class JEDi_Get_Residue_Pair_Coupling: Gets the interaction between pairs of residues using hierarchical pca on each residue pair. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Get_Residue_Pair_Coupling
{
	double trace, coupling_score;
	Matrix U_evectors;
	final double FLOOR, NOISE_LEVEL;
	final int number_of_eigenresidues, number_of_residues, ROWS, COLS;
	final String Residue_Pair_Index;
	final List<Double> eigenvalues;
	final Matrix cov;
	final Matrix residuePCs;
	final EigenvalueDecomposition evd;
	final PCA pca;
	NumberFormat nf0, nf3, nf6, df;
	RoundingMode rm;

	/**
	 * Constructor to perform PCA on generalized residue coordinates
	 *
	 * @param data      The subset coordinates
	 * @param num_modes The number of PCA modes to process
	 * @param dir       The working directory
	 * @param desc      The job description
	 */
	JEDi_Get_Residue_Pair_Coupling(Matrix res_Pair_PCs, String res_pair_index)
	{
		this.residuePCs = res_Pair_PCs;
		this.Residue_Pair_Index = res_pair_index;

		this.FLOOR = Input_Parameters.FLOOR;
		this.NOISE_LEVEL = Input_Parameters.NOISE_LEVEL;
		this.number_of_eigenresidues = Input_Parameters.MODES_EIGEN_RESIDUE_PAIRS;

		this.ROWS = residuePCs.getRowDimension();
		this.COLS = residuePCs.getColumnDimension();
		this.number_of_residues = (ROWS / number_of_eigenresidues);

		this.pca = new PCA(residuePCs);
		this.cov = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		this.evd = pca.get_eigenvalue_decomposition(cov);
		this.eigenvalues = new ArrayList<Double>();
		rm = RoundingMode.HALF_UP;
		df = new DecimalFormat("0.###E0");
		df.setRoundingMode(rm);

		nf0 = NumberFormat.getInstance();
		nf0.setRoundingMode(rm);
		nf0.setMaximumFractionDigits(0);
		nf0.setMinimumFractionDigits(0);

		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);
	}

	/* **************************************************** PUBLIC METHOD ****************************************************** */

	/**
	 * Performs Pairwise Residue analysis.
	 */
	public void get_Pairwise_Couping_Analysis()
	{
		get_eigenvalues();
		get_U_evects_and_reverse();
		get_Score();
	}
	/* **************************************************** PRIVATE METHODS ****************************************************** */

	private void get_eigenvalues()
	{
		double[] ss_evals = evd.getRealEigenvalues();
		trace = 0;
		for (double k : ss_evals)
			{
				if (k < FLOOR) k = FLOOR;
				eigenvalues.add(k);
				trace += k;
			}
		Collections.sort(eigenvalues, Collections.reverseOrder());
	}

	private void get_U_evects_and_reverse()
	{
		Matrix ss_evectors = evd.getV();
		U_evectors = ss_evectors;
		Matrix modes_reversed = new Matrix(ROWS, ROWS);
		for (int i = 0; i < ROWS; i++)
			{
				Matrix col = U_evectors.getMatrix(0, ROWS - 1, ROWS - 1 - i, ROWS - 1 - i);
				modes_reversed.setMatrix(0, ROWS - 1, i, i, col);
			}
		U_evectors = modes_reversed;
	}

	private double get_Score()
	{
		coupling_score = 0;
		for (int i = 0; i < 2 * number_of_eigenresidues; i++) // Loop iterating over the eigenvectors.
			{
				double lambda = eigenvalues.get(i);
				Matrix col1 = U_evectors.getMatrix(0, number_of_eigenresidues - 1, i, i);
				Matrix col2 = U_evectors.getMatrix(number_of_eigenresidues, 2 * number_of_eigenresidues - 1, i, i);
				double w1 = Math.pow(col1.normF(), 2);
				double w2 = Math.pow(col2.normF(), 2);
				double scale = lambda / trace;
				// double interaction = Math.abs(w1 - w2);
				// double Cij = (scale) * (1 - interaction); // old scoring function
				double score = (scale) * Math.exp(-8 * (w2 - w1) * (w2 - w1)); // Gaussian scoring function
				coupling_score += score;
				// coupling_score += Cij;
			}
		coupling_score *= 100;
		// System.out.println(nf3.format(coupling_score));
		return coupling_score;
	}

	/* ******************************************* GETTERS ***************************************************** */

	public double getCoupling_score()
	{
		return coupling_score;
	}

	public double get_trace()
	{
		return trace;
	}

	public List<Double> getEigenvalues()
	{
		return eigenvalues;
	}

	public Matrix getTop_evectors()
	{
		return U_evectors;
	}
}
