package jedi;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class Residue_RMSD: This class computes the Residue RMSDs (RMSFs) from an entire trajectory. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Residue_RMSFs
	{

		int COLS, ROWS, number_of_atoms;
		Matrix coordinates, means, sum_of_sq_dev, sigmas, z_scores;

		/* ************************************** CONSTRUCTORS ******************************************************************************** */

		/**
		 * Constructor that uses the PDB coordinates to calculate the residue RMSFs.
		 */
		JEDi_Get_Residue_RMSFs(Matrix input)
			{

				coordinates = input;
				COLS = coordinates.getColumnDimension();
				ROWS = coordinates.getRowDimension();
				number_of_atoms = (ROWS / 3);
				means = new Matrix(ROWS, 1);
				sum_of_sq_dev = new Matrix(ROWS, 1);
				sigmas = new Matrix(ROWS, 1);
				z_scores = new Matrix(ROWS, COLS);
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
						double mean = Descriptive_Stats.get_mean(row);
						double ssdevs = Descriptive_Stats.get_sum_of_squared_deviations(row, mean);
						double sigma = Descriptive_Stats.get_standard_deviation(row, mean, ssdevs);
						means.set(i, 0, mean);
						sum_of_sq_dev.set(i, 0, ssdevs);
						sigmas.set(i, 0, sigma);
						double[] z_s = Descriptive_Stats.get_Z_scores(row, mean, sigma);
						Matrix z = new Matrix(z_s, 1);
						z_scores.setMatrix(i, i, 0, COLS - 1, z);
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
		 * @return A Java List of RMSFs based on the entire trajectory.
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

		/**
		 * @return The variable Z scores.
		 */
		public Matrix get_z_scores()
			{
				return z_scores;
			}
	}
