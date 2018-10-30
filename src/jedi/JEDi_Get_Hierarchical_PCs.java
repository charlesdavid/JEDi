package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Generalized_Residue_Coordinate_PCs: Constructs the PCs for the RGC subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Hierarchical_PCs
	{

		String directory, out_dir, description, file_name_head, path;
		int number_of_modes, number_of_modes_residue, number_of_conformations;
		int ROWS, COLS, number_of_residues;
		List<Double> eigenvalues;
		Matrix centered_coords, top_evectors;
		Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;
		boolean exist, success;

		/**
		 * Constructor to create the PCs.
		 *
		 */
		JEDi_Get_Hierarchical_PCs(Matrix centered_data, Matrix evects, List<Double> evals, String dir, String des, int modes_resi)
			{

				this.centered_coords = centered_data;
				this.top_evectors = evects;
				this.eigenvalues = evals;
				this.directory = dir;
				this.description = des;
				this.number_of_modes_residue = modes_resi;
				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Hierarchical_PCA" + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();

				number_of_modes = top_evectors.getColumnDimension();
				number_of_conformations = centered_coords.getColumnDimension();
				ROWS = top_evectors.getRowDimension();
				number_of_residues = (ROWS / number_of_modes_residue);

				file_name_head = out_dir + "ss_" + number_of_residues;
			}

		/**
		 * Computes the Principle Components and writes them to files: Un-normed, normed, weighted, weighted-normed
		 */
		void get_PCs()
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
								Matrix vector1 = Projector.get_Normed_array(data1);
								Matrix vector2 = Projector.get_Normed_array(data2);
								double dp = Projector.get_InnerProduct(data1, data2);
								double normed_dp = Projector.get_InnerProduct(vector1, vector2);
								if (dp == Double.NaN)
									dp = 0.000;
								if (normed_dp == Double.NaN)
									normed_dp = 0.000;
								double w_dp = weight * dp;
								double weighted_normed_dp = weight * normed_dp;
								projections.set(inner, outer, dp);
								normed_projections.set(inner, outer, normed_dp);
								weighted_projections.set(inner, outer, w_dp);
								weighted_normed_projections.set(inner, outer, weighted_normed_dp);
							}
					}
				path = file_name_head + "_top_" + number_of_modes + "_PCs.txt";
				Matrix_IO.write_Matrix(projections, path, 12, 9);
				path = file_name_head + "_top_" + number_of_modes + "_normed_PCs.txt";
				Matrix_IO.write_Matrix(normed_projections, path, 12, 9);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_PCs.txt";
				Matrix_IO.write_Matrix(weighted_projections, path, 12, 9);
				path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_PCs.txt";
				Matrix_IO.write_Matrix(weighted_normed_projections, path, 12, 9);
				System.gc();
			}

		/* ******************************************* GETTERS ***************************************************** */

		public Matrix getPCs()
			{
				return projections;
			}

		public Matrix getNormed_PCs()
			{
				return normed_projections;
			}

		public Matrix getWeighted_normed_PCs()
			{
				return weighted_normed_projections;
			}

		public Matrix getWeighted_PCs()
			{
				return weighted_projections;
			}
	}
