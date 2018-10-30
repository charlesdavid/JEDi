package jedi;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class Adjust_Outliers_by_Z_Score: Adjusts variable outliers using a Z score cutoff. Copyright (C) 2012 Dr. Charles David
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
 */

public class Adjust_Outliers_by_Z_Score
	{

		int COLS, ROWS;
		static double z_threshold;
		List<Double> residue_rmsd_list;
		Matrix coordinates, coordinates_adjusted, z_scores, counts;
		double[] means, std_deviations;

		/* *************************** CONSTRUCTORS *************************************************************************************** */

		/**
		 * Constructor takes the original matrix of coordinates
		 * 
		 * @param input_data The matrix of coordinates
		 */
		public Adjust_Outliers_by_Z_Score(Matrix input_data)
			{

				coordinates = input_data;
				COLS = coordinates.getColumnDimension();
				ROWS = coordinates.getRowDimension();
				z_scores = new Matrix(ROWS, COLS);
				coordinates_adjusted = new Matrix(ROWS, COLS);
				residue_rmsd_list = new ArrayList<Double>();
				counts = new Matrix(ROWS, 1);
				means = new double[ROWS];
				std_deviations = new double[ROWS];
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
						double mean = Descriptive_Stats.get_mean(row);
						double ssdevs = Descriptive_Stats.get_sum_of_squared_deviations(row, mean);
						double sigma = Descriptive_Stats.get_standard_deviation(row, mean, ssdevs);
						double[] z_s = Descriptive_Stats.get_Z_scores(row, mean, sigma);
						int index = 0;
						int count = 0;
						for (double d : z_s)
							{
								double score = Math.abs(d);
								if (score >= z_threshold)
									{
										coordinates_adjusted.set(i, index, mean);
										count++;
									} else
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

		/* *************************** SETTERS *************************************************************************************** */

		/**
		 * Sets the Z-Cutoff
		 * 
		 * @param z_cutoff
		 */
		public void set_Z_threshold(double z_cutoff)
			{

				Adjust_Outliers_by_Z_Score.z_threshold = z_cutoff;
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

		/**
		 * Returns the Z-score cutoff
		 * 
		 * @return z_threshold
		 */
		public double get_Z_threshold()
			{

				return z_threshold;
			}

	}
