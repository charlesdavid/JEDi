package support;

import Jama.Matrix;

/**
 * Row_Center_Data: This class centers the data in a matrix prior to PCA.
 * 
 * Note: The ROWS are the variables, the COLS are the frames.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */
public class Row_Center_Data
{
	final int ROWS, COLS;
	final Matrix data, row_centered_data, row_variances, row_standard_deviations, row_means;
	final Descriptive_Stats ds;

	/**
	 * Constructor to center the data prior to PCA based on rows being the variables.
	 * 
	 * @param data_matrix The un-centered data
	 */
	public Row_Center_Data(Matrix data_matrix)
	{
		data = data_matrix;
		ROWS = data.getRowDimension();
		COLS = data.getColumnDimension();
		row_centered_data = new Matrix((ROWS), (COLS));
		row_means = new Matrix(ROWS, (1));
		row_variances = new Matrix(ROWS, 1);
		row_standard_deviations = new Matrix(ROWS, 1);
		ds = new Descriptive_Stats();
		center_data();
	}

	/**
	 * Method that mean centers the data in each row of the input matrix.
	 */
	private void center_data()
	{
		for (int j = 0; j < (ROWS); j++)
			{
				Matrix row = data.getMatrix(j, j, 0, COLS - 1);
				double[] row_data = row.getRowPackedCopy();
				double mean = ds.get_mean(row_data);
				row_means.set(j, 0, mean);
				double variance = ds.get_variance(row_data, mean);
				row_variances.set(j, 0, variance);
				double stdev = ds.get_standard_deviation(row_data, mean);
				row_standard_deviations.set(j, 0, stdev);
				for (int i = 0; i < (COLS); i++)
					{
						double centered_row_element = (data.get(j, i) - mean);
						row_centered_data.set(j, i, centered_row_element);
					}
				if (j % 10 == 0) System.gc();
			}
	}

	/**
	 * @return The centered data
	 */
	public Matrix get_row_centered_data()
	{
		return row_centered_data;
	}

	/**
	 * @return The means of the variables (ROWS)
	 */
	public Matrix get_variable_means()
	{
		return row_means;
	}

	/**
	 * @return The Variances of the variables (ROWS)
	 */
	public Matrix get_variable_variances()
	{
		return row_variances;
	}

	/**
	 * @return The Standard Deviations of the variables (ROWS)
	 */
	public Matrix get_variable_standard_deviations()
	{
		return row_standard_deviations;
	}
}
