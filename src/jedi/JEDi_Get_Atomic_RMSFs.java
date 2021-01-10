package jedi;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import support.Descriptive_Stats;

/**
 * JED class Residue_RMSD: This class computes the Residue RMSDs (RMSFs) from an entire trajectory. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Atomic_RMSFs
{
	final int COLS, ROWS, number_of_atoms;
	final Matrix coordinates, means, sum_of_sq_dev, sigmas;
	final Descriptive_Stats ds;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Constructor that uses the PDB coordinates to calculate the residue RMSFs.
	 */
	JEDi_Get_Atomic_RMSFs(Matrix input)
	{
		this.coordinates = input;
		this.COLS = coordinates.getColumnDimension();
		this.ROWS = coordinates.getRowDimension();
		this.number_of_atoms = (ROWS / 3);
		this.means = new Matrix(ROWS, 1);
		this.sum_of_sq_dev = new Matrix(ROWS, 1);
		this.sigmas = new Matrix(ROWS, 1);
		this.ds = new Descriptive_Stats();
	}

	/* ************************************** METHODS ******************************************************************************** */

	/**
	 * Method that determines the residue statistics.
	 */
	private void get_residue_stats()
	{
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double mean = ds.get_mean(row);
				double ssdevs = ds.get_sum_of_squared_deviations(row, mean);
				double sigma = ds.get_standard_deviation(row, mean, ssdevs);
				means.set(i, 0, mean);
				sum_of_sq_dev.set(i, 0, ssdevs);
				sigmas.set(i, 0, sigma);
			}
	}

	/**
	 * Method that returns a matrix of RMSFs.
	 */
	private Matrix get_residue_rmsfs_matrix()
	{
		get_residue_stats();

		Matrix res_rmsfs = new Matrix(number_of_atoms, 1);
		for (int i = 0; i < number_of_atoms; i++)
			{
				double x = sum_of_sq_dev.get(i, 0);
				double y = sum_of_sq_dev.get(i + number_of_atoms, 0);
				double z = sum_of_sq_dev.get(i + 2 * number_of_atoms, 0);
				double sum = ((x + y + z) / COLS);
				double rmsd = Math.sqrt(sum);
				res_rmsfs.set(i, 0, rmsd);
			}
		return res_rmsfs;
	}

	/**
	 * Method that returns a List of RMSFs
	 */
	public List<Double> get_residue_rmsfs()
	{
		Matrix rmsfs = get_residue_rmsfs_matrix();
		ArrayList<Double> residue_rmsfs = new ArrayList<>();
		double[] rrmsfs = rmsfs.getColumnPackedCopy();
		for (double d : rrmsfs)
			{
				residue_rmsfs.add(d);
			}
		return residue_rmsfs;
	}
}
