package support;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import supportIO.Input_Parameters;

/**
 * Adjust Outliers by Median Absolute Deviation: Adjusts variable outliers using a MAD cutoff.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Adjust_Outliers_by_MAD
{
	final int COLS, ROWS;
	final double mad_threshold;
	final List<Double> residue_rmsd_list;
	final Matrix coordinates, coordinates_adjusted, mad_scores, counts;
	double[] medians, mads;
	final Descriptive_Stats ds;

	/* *************************** CONSTRUCTOR *************************************************************************************** */

	/**
	 * Constructor takes the original matrix of coordinates
	 * 
	 * @param input_data The matrix of coordinates
	 */
	public Adjust_Outliers_by_MAD(Matrix input_data)
	{
		this.coordinates = input_data;
		this.COLS = coordinates.getColumnDimension();
		this.ROWS = coordinates.getRowDimension();
		this.mad_scores = new Matrix(ROWS, COLS);
		this.coordinates_adjusted = new Matrix(ROWS, COLS);
		this.residue_rmsd_list = new ArrayList<Double>();
		this.counts = new Matrix(ROWS, 1);
		this.medians = new double[ROWS];
		this.mads = new double[ROWS];
		this.mad_threshold = Input_Parameters.MAD_SCORE_CUTOFF;
		this.ds = new Descriptive_Stats();
	}

	/* *************************** METHODS *************************************************************************************** */

	/**
	 * Sets all variable values beyond |MAD-cut| to their median value
	 */
	public void adjust_row_data()
	{
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double median = ds.get_median(row);
				double mad = ds.get_mad(row);
				double[] madScores = ds.get_mad_scores(row, median, mad);
				int index = 0;
				int count = 0;
				for (double d : madScores)
					{
						double score = Math.abs(d);
						if (score <= mad_threshold)
							{
								coordinates_adjusted.set(i, index, median);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted.set(i, index, coord);
							}
						index++;
					}
				Matrix z = new Matrix(madScores, 1);
				mad_scores.setMatrix(i, i, 0, COLS - 1, z);
				medians[i] = median;
				mads[i] = mad;
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
	 * Gets the matrix of MAD scores, which has the same dimensions as the coordinates matrix.
	 * 
	 * @return mad_scores
	 */
	public Matrix get_mad_scores()
	{

		return mad_scores;
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
	 * Returns the medians
	 * 
	 * @return medians
	 */
	public double[] get_medians()
	{

		return medians;
	}

	/**
	 * Returns the median absolute deviations
	 * 
	 * @return mads
	 */
	public double[] get_mads()
	{

		return mads;
	}
}
