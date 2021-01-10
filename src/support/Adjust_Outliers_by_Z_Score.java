package support;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import supportIO.Input_Parameters;

/**
 * Adjust_Outliers_by_Z_Score: Adjusts variable outliers using a Z score cutoff.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Adjust_Outliers_by_Z_Score
{
	final int COLS, ROWS;
	final double z_threshold;
	final List<Double> residue_rmsd_list;
	final Matrix coordinates, coordinates_adjusted, z_scores, counts;
	final double[] means, std_deviations;
	final Descriptive_Stats ds;

	/* *************************** CONSTRUCTOR *************************************************************************************** */

	/**
	 * Constructor takes the original matrix of coordinates
	 * 
	 * @param input_data The matrix of coordinates
	 */
	public Adjust_Outliers_by_Z_Score(Matrix input_data)
	{
		this.coordinates = input_data;
		this.COLS = coordinates.getColumnDimension();
		this.ROWS = coordinates.getRowDimension();
		this.z_scores = new Matrix(ROWS, COLS);
		this.coordinates_adjusted = new Matrix(ROWS, COLS);
		this.residue_rmsd_list = new ArrayList<Double>();
		this.counts = new Matrix(ROWS, 1);
		this.means = new double[ROWS];
		this.std_deviations = new double[ROWS];
		this.z_threshold = Input_Parameters.Z_SCORE_CUTOFF;
		this.ds = new Descriptive_Stats();
	}

	/* *************************** METHODS *************************************************************************************** */

	/**
	 * Sets all variable values beyond |z-cut| to their mean value
	 */
	public void adjust_row_data()
	{
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double mean = ds.get_mean(row);
				double ssdevs = ds.get_sum_of_squared_deviations(row, mean);
				double sigma = ds.get_standard_deviation(row, mean, ssdevs);
				double[] z_s = ds.get_Z_scores(row, mean, sigma);

				int index = 0;
				int count = 0;
				for (double d : z_s)
					{
						double score = Math.abs(d);
						if (score <= z_threshold)
							{
								coordinates_adjusted.set(i, index, mean);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted.set(i, index, coord);
							}
						index++;
					}
				Matrix z = new Matrix(z_s, 1);
				z_scores.setMatrix(i, i, 0, COLS - 1, z);
				means[i] = mean;
				std_deviations[i] = sigma;
				counts.set(i, 0, count);
			}
	}

	/* *************************** GETTERS *************************************************************************************** */

	/**
	 * Returns the matrix of coordinates with outliers adjusted
	 * 
	 * @return
	 */
	public Matrix get_coorinates_adjusted()
	{

		return coordinates_adjusted;
	}

	/**
	 * Gets the matrix of Z scores, which has the same dimensions as the coordinates matrix.
	 * 
	 * @return z_scores
	 */
	public Matrix get_z_scores()
	{

		return z_scores;
	}

	/**
	 * Returns the number of variable adjustments (per variable)
	 * 
	 * @return The matrix of the number of adjustments per variable
	 */
	public Matrix get_counts()
	{

		return counts;
	}

	/**
	 * Returns the means
	 * 
	 * @return means
	 */
	public double[] get_means()
	{

		return means;
	}

	/**
	 * Returns the standard deviations
	 * 
	 * @return std_deviations
	 */
	public double[] get_std_deviations()
	{

		return std_deviations;
	}
}
