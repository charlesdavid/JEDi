package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Cartesian_PCs: Constructs the Principle Components for the Cartesian PCA on the residues. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Residue_PCs
	{

		String directory, out_dir, description, model, file_name_head, path;
		int number_of_modes, number_of_conformations, ROWS, COLS, number_of_atoms, Residue_Index;
		List<Double> eigenvalues;
		Matrix centered_coords, top_evectors, projections, normed_projections, weighted_normed_projections, weighted_projections;
		boolean exist, success;

		/* ************************************** CONSTRUCTOR ******************************************************************************** */

		/**
		 * Constructor to create the delta vector series and the delta vector projections.
		 *
		 * @param centered_data The mean centered data
		 * @param evects        The top Cartesian eigenvectors
		 * @param evals         The top Cartesian eigenvalues
		 * @param dir           The working directory
		 * @param des           The job description
		 * @param pca_model     The type of PCA: COV, CORR, PCORR (Q, R, P)
		 */
		JEDi_Get_Residue_PCs(Matrix centered_data, Matrix evects, List<Double> evals, String dir, String des, String pca_model)
			{

				centered_coords = centered_data;
				top_evectors = evects;
				eigenvalues = evals;
				directory = dir;
				description = des;
				model = pca_model;

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Residue_cPCA" + File.separatorChar + model + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = new File(out_dir).mkdirs();

				number_of_modes = top_evectors.getColumnDimension();
				number_of_conformations = centered_coords.getColumnDimension();
				ROWS = top_evectors.getRowDimension();
				number_of_atoms = (ROWS / 3);
			}

		/* ************************************** Public Methods ******************************************************************************** */

		/**
		 * Computes the Principle Components and writes them to files: Un-normed, normed, weighted, weighted-normed
		 */
		public void get_PCs()
			{
				projections = new Matrix(number_of_conformations, number_of_modes);
				normed_projections = new Matrix(number_of_conformations, number_of_modes);
				weighted_projections = new Matrix(number_of_conformations, number_of_modes);
				weighted_normed_projections = new Matrix(number_of_conformations, number_of_modes);

				for (int outer = 0; outer < number_of_modes; outer++)
					{
						for (int inner = 0; inner < number_of_conformations; inner++)
							{
								int row_index_1 = 0;
								int row_index_2 = ROWS - 1;
								Matrix data1 = top_evectors.getMatrix(row_index_1, row_index_2, outer, outer);
								Matrix data2 = centered_coords.getMatrix(row_index_1, row_index_2, inner, inner);
								double weight = eigenvalues.get(outer);
								if (weight < 0)
									weight = 0; // Fixes the -0.00000000 problem.
								double sqrt_weight = Math.sqrt(weight);
								Matrix vector1 = Projector.get_Normed_array(data1);
								Matrix vector2 = Projector.get_Normed_array(data2);
								double dp = Projector.get_InnerProduct(data1, data2);
								double normed_dp = Projector.get_InnerProduct(vector1, vector2);
								if (dp == Double.NaN)
									dp = 0.000;
								if (normed_dp == Double.NaN)
									normed_dp = 0.000;
								double w_dp = sqrt_weight * dp;
								double weighted_normed_dp = sqrt_weight * normed_dp;
								projections.set(inner, outer, dp);
								normed_projections.set(inner, outer, normed_dp);
								weighted_projections.set(inner, outer, w_dp);
								weighted_normed_projections.set(inner, outer, weighted_normed_dp);
							}
					}
				String res_index = String.format("%03d", (Residue_Index + 1));
				file_name_head = out_dir + "Residue_" + res_index;

				path = file_name_head + "_top_" + number_of_modes + "_PCs_" + model + ".txt";
				Matrix_IO.write_Matrix(projections, path, 9, 6);

				path = file_name_head + "_top_" + number_of_modes + "_normed_PCs_" + model + ".txt";
				Matrix_IO.write_Matrix(normed_projections, path, 9, 6);

				path = file_name_head + "_top_" + number_of_modes + "_weighted_PCs_" + model + ".txt";
				Matrix_IO.write_Matrix(weighted_projections, path, 9, 6);

				path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_PCs_" + model + ".txt";
				Matrix_IO.write_Matrix(weighted_normed_projections, path, 9, 6);

				System.gc();
			}

		/* ************************************** SETTERS ******************************************************************************** */

		public void set_Residue_Index(int residue_Index)
			{
				Residue_Index = residue_Index;
			}

		/* ************************************** GETTERS ******************************************************************************** */

		public Matrix getProjections()
			{
				return projections;
			}

		public Matrix getNormed_projections()
			{
				return normed_projections;
			}

		public Matrix getWeighted_normed_projections()
			{
				return weighted_normed_projections;
			}

		public Matrix getWeighted_projections()
			{
				return weighted_projections;
			}

		public int get_Residue_Index()
			{
				return Residue_Index;
			}
	}
