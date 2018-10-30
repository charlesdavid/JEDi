package jedi;

import Jama.Matrix;

/**
 * JED class Row_Center_Data: This class centers the data in a matrix prior to PCA. Note: The ROWS are the variables, the COLS are the values. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

public class Row_Center_Data
	{

		int ROWS, COLS;
		Matrix data, row_centered_data, row_standard_deviations, row_means;

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
				row_standard_deviations = new Matrix(ROWS, 1);
				center_data();
			}

		/**
		 * Method that centers the data along a row, subtracting off the mean.
		 */
		private void center_data()
			{

				for (int j = 0; j < (ROWS); j++)
					{
						Matrix row = data.getMatrix(j, j, 0, COLS - 1);
						double[] row_data = row.getRowPackedCopy();
						double mean = Descriptive_Stats.get_mean(row_data);
						row_means.set(j, 0, mean);
						double stdev = Descriptive_Stats.get_standard_deviation(row_data, mean);
						row_standard_deviations.set(j, 0, stdev);
						for (int i = 0; i < (COLS); i++)
							{
								double centered_row_element = (data.get(j, i) - mean);
								row_centered_data.set(j, i, centered_row_element);
							}
						if (j % 10 == 0)
							System.gc();
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
		 * @return The Standard Deviations of the variables (ROWS)
		 */
		public Matrix get_variable_sigmas()
			{

				return row_standard_deviations;
			}
	}
