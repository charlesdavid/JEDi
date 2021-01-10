package jedi;

import java.util.List;

import Jama.Matrix;

public class JEDi_Get_PCA_Modes
{
	final int ROWS, COLS, number_atoms;
	final List<Double> eigenvalues;
	double[] pca_mode_mins, pca_mode_maxs;
	final Matrix evectors;
	Matrix pca_modes, square_pca_modes, weighted_square_pca_modes, weighted_pca_modes;

	public JEDi_Get_PCA_Modes(Matrix evects, List<Double> evals)
	{
		super();
		this.evectors = evects;
		this.eigenvalues = evals;
		this.ROWS = evectors.getRowDimension();
		this.COLS = evectors.getColumnDimension();
		this.number_atoms = (ROWS / 3);
		construct_PCA_Modes();
	}

	private void construct_PCA_Modes()
	{
		pca_modes = new Matrix(number_atoms, COLS);
		weighted_pca_modes = new Matrix(number_atoms, COLS);
		square_pca_modes = new Matrix(number_atoms, COLS);
		weighted_square_pca_modes = new Matrix(number_atoms, COLS);
		pca_mode_maxs = new double[COLS];
		pca_mode_mins = new double[COLS];
		for (int i = 0; i < COLS; i++)
			{
				double max = 0;
				double min = 1;
				for (int j = 0; j < number_atoms; j++)
					{
						double x = evectors.get(j, i);
						double y = evectors.get(j + number_atoms, i);
						double z = evectors.get((j + 2 * number_atoms), i);
						double sq = (x * x + y * y + z * z);
						double sqrt_sq = Math.sqrt(sq);
						double value = eigenvalues.get(i);
						double sqrt_val = Math.sqrt(value);
						double w_sq = sq * value;
						double w_m = sqrt_sq * sqrt_val;
						pca_modes.set(j, i, sqrt_sq);
						weighted_pca_modes.set(j, i, w_m);
						square_pca_modes.set(j, i, sq);
						weighted_square_pca_modes.set(j, i, w_sq);
						if (sq > max) max = sq;
						if (sq < min) min = sq;
					}
				pca_mode_maxs[i] = max;
				pca_mode_mins[i] = min;
			}
	}

	/* ******************************************* GETTERS ***************************************************** */

	public double[] get_PCA_mode_mins()
	{
		return pca_mode_mins;
	}

	public double[] get_PCA_mode_maxs()
	{
		return pca_mode_maxs;
	}

	public Matrix get_PCA_modes()
	{
		return pca_modes;
	}

	public Matrix get_Square_PCA_modes()
	{
		return square_pca_modes;
	}

	public Matrix get_Weighted_Square_PCA_modes()
	{
		return weighted_square_pca_modes;
	}

	public Matrix get_Weighted_PCA_modes()
	{
		return weighted_pca_modes;
	}
}
