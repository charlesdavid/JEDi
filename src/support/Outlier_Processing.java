package support;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import supportIO.Input_Parameters;

/**
 * Class for Processing Outliers:
 * 
 * Outliers can be REMOVED or SELECTED.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Outlier_Processing
{
	final int COLS, ROWS;
	double z_threshold, mad_threshold;
	final List<Double> residue_rmsd_list;
	final Matrix coordinates, counts_remove_MAD, counts_remove_Z, counts_select_MAD, counts_select_Z, MAD_scores, Z_scores;
	final double[] means, medians, mads, std_deviations, sum_square_deviations, variances;
	final Descriptive_Stats ds;

	/* *************************** CONSTRUCTOR *************************************************************************************** */

	/**
	 * Constructor takes the original matrix of coordinates
	 * 
	 * @param input_data The matrix of coordinates
	 */
	public Outlier_Processing(Matrix input_data)
	{
		this.coordinates = input_data;
		this.COLS = coordinates.getColumnDimension();
		this.ROWS = coordinates.getRowDimension();
		this.ds = new Descriptive_Stats();

		this.MAD_scores = new Matrix(ROWS, COLS);
		this.Z_scores = new Matrix(ROWS, COLS);
		this.residue_rmsd_list = new ArrayList<Double>();
		this.counts_remove_MAD = new Matrix(ROWS, 1);
		this.counts_remove_Z = new Matrix(ROWS, 1);
		this.counts_select_MAD = new Matrix(ROWS, 1);
		this.counts_select_Z = new Matrix(ROWS, 1);
		this.means = new double[ROWS];
		this.medians = new double[ROWS];
		this.sum_square_deviations = new double[ROWS];
		this.std_deviations = new double[ROWS];
		this.variances = new double[ROWS];
		this.mads = new double[ROWS];

		this.mad_threshold = Input_Parameters.MAD_SCORE_CUTOFF;
		this.z_threshold = Input_Parameters.Z_SCORE_CUTOFF;
	}

	/* *************************** SETTERS *************************************************************************************** */

	public void set_Z_threshold(double z)
	{
		this.z_threshold = z;
	}

	public void set_Mad_threshold(double mad)
	{
		this.mad_threshold = mad;
	}

	/* *************************** METHODS *************************************************************************************** */


	/**
	 * Get the statistics for each variable in the subset.
	 */
	public Matrix getCoordinateStats()
	{
		Matrix stats = new Matrix(ROWS, 4);
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double mean = ds.get_mean(row);
				double median = ds.get_median(row);
				double mad = ds.get_mad(row);
				double skew = ds.get_skew(row);
				double kurtosis = ds.get_kurtosis(row);
				double ssdevs = ds.get_sum_of_squared_deviations(row, mean);
				double std_dev = ds.get_standard_deviation(row, mean, ssdevs);
				double var = std_dev * std_dev;

				means[i] = mean;
				medians[i] = median;
				sum_square_deviations[i] = ssdevs;
				std_deviations[i] = std_dev;
				variances[i] = var;
				mads[i] = mad;

				stats.set(i, 0, mean);
				stats.set(i, 1, var);
				stats.set(i, 2, skew);
				stats.set(i, 3, kurtosis);

				// System.out.println("Mean: " + mean + " Variance: " + var + " Skew: " + skew + " Kurtosis: " + kurtosis);
			}
		return stats;
	}


	/**
	 * Sets all variable values with MAD Score GREATER THAN |MAD-cut| to their median value. This effectively removes the outliers without discarding any data.
	 */
	public Matrix remove_outliers_MAD()
	{
		Matrix coordinates_adjusted_mad = new Matrix(ROWS, COLS);
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double median = medians[i];
				double mad = mads[i];

				double[] madScores = ds.get_mad_scores(row, median, mad);
				int index = 0;
				int count = 0;
				for (double d : madScores)
					{
						double score = Math.abs(d);
						if (score >= mad_threshold)
							{
								coordinates_adjusted_mad.set(i, index, median);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted_mad.set(i, index, coord);
							}
						index++;
					}
				Matrix z = new Matrix(madScores, 1);
				MAD_scores.setMatrix(i, i, 0, COLS - 1, z);
				counts_remove_MAD.set(i, 0, count);
			}
		return coordinates_adjusted_mad;
	}

	/**
	 * Sets all variable values with Z-Score GREATER THAN |z-cut| to their mean value. This effectively removes the outliers without discarding any data.
	 */
	public Matrix remove_outliers_Z()
	{
		Matrix coordinates_adjusted_Z = new Matrix(ROWS, COLS);
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double mean = means[i];
				double sigma = std_deviations[i];
				double[] z_s = ds.get_Z_scores(row, mean, sigma);
				int index = 0;
				int count = 0;
				for (double d : z_s)
					{
						double score = Math.abs(d);

						if (score >= z_threshold)
							{
								coordinates_adjusted_Z.set(i, index, mean);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted_Z.set(i, index, coord);
							}
						index++;
					}
				Matrix z = new Matrix(z_s, 1);
				Z_scores.setMatrix(i, i, 0, COLS - 1, z);
				counts_remove_Z.set(i, 0, count);
			}
		return coordinates_adjusted_Z;
	}


	/**
	 * Sets all variable values LESS THAN |MAD-cut| to their median value. This effectively emphasizes the affect of the outliers without discarding any data.
	 */
	public Matrix select_outliers_MAD()
	{
		Matrix coordinates_adjusted_mad = new Matrix(ROWS, COLS);
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double median = medians[i];
				double mad = mads[i];
				double[] madScores = ds.get_mad_scores(row, median, mad);
				int index = 0;
				int count = 0;
				for (double d : madScores)
					{
						double score = Math.abs(d);
						if (score < mad_threshold)
							{
								coordinates_adjusted_mad.set(i, index, median);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted_mad.set(i, index, coord);
							}
						index++;
					}
				counts_select_MAD.set(i, 0, count);
			}
		return coordinates_adjusted_mad;
	}

	/**
	 * Sets all variable values LESS THAN |z-cut| to their mean value. This effectively emphasizes the affect of the outliers without discarding any data.
	 */
	public Matrix select_outliers_Z()
	{
		Matrix coordinates_adjusted_Z = new Matrix(ROWS, COLS);
		for (int i = 0; i < ROWS; i++)
			{
				double[] row = coordinates.getMatrix(i, i, 0, COLS - 1).getRowPackedCopy();
				double mean = means[i];
				double sigma = std_deviations[i];
				double[] z_s = ds.get_Z_scores(row, mean, sigma);
				int index = 0;
				int count = 0;
				for (double d : z_s)
					{
						double score = Math.abs(d);
						if (score < z_threshold)
							{
								coordinates_adjusted_Z.set(i, index, mean);
								count++;
							}
						else
							{
								double coord = row[index];
								coordinates_adjusted_Z.set(i, index, coord);
							}
						index++;
					}
				counts_select_Z.set(i, 0, count);
			}
		return coordinates_adjusted_Z;
	}

	/* *************************** GETTERS *************************************************************************************** */


	/**
	 * Gets the matrix of MAD scores, which has the same dimensions as the coordinates matrix.
	 * 
	 * @return mad_scores
	 */
	public Matrix get_mad_scores()
	{
		return MAD_scores;
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

	/**
	 * Gets the matrix of Z scores, which has the same dimensions as the coordinates matrix.
	 * 
	 * @return z_scores
	 */
	public Matrix get_z_scores()
	{
		return Z_scores;
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

	public Matrix getCounts_remove_MAD()
	{
		return counts_remove_MAD;
	}

	public Matrix getCounts_remove_Z()
	{
		return counts_remove_Z;
	}

	public Matrix getCounts_select_MAD()
	{
		return counts_select_MAD;
	}

	public Matrix getCounts_select_Z()
	{
		return counts_select_Z;
	}
}
