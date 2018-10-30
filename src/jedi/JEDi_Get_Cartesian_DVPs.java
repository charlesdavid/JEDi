package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Cartesian_DVPs: Constructs the DVPs for the Cartesian subset. Note: DVPs are like Principle Components, except that the reference structure is not the mean
 * structure. Copyright (C) 2012 Dr. Charles David
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Cartesian_DVPs
	{

		String directory, out_dir, description, type, model, file_name_head, path;
		final String AA = "All_Atom_cPCA", AC = "Alpha_Carbon_cPCA", Res = "Residue_cPCA", Q = "COV", R = "CORR", P = "PCORR";
		int number_of_modes, number_of_conformations, ROWS, COLS, number_of_atoms;
		List<Double> eigenvalues;
		Matrix reference_PDB_coordinates, input_coords, delta_vector_series, top_evectors;
		Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;
		long startTime, endTime, totalTime;
		boolean exist;
		boolean success;

		/**
		 * Constructor to create the delta vector series and the delta vector projections.
		 *
		 * @param data       The transformed Coordinates Matrix
		 * @param ref_coords The Reference coordinates from the reference PDB file
		 * @param evects     The top Cartesian eigenvectors
		 * @param evals      The top Cartesian eigenvalues
		 * @param dir        The working directory
		 * @param des        The job description
		 * @param type_pca   The type of PCA: All-Atom, Alpha-Carbon, Residue (AA, AC, Res)
		 * @param model_pca  The PCA model: COV, CORR, PCORR (Q, R, P)
		 */
		JEDi_Get_Cartesian_DVPs(Matrix data, Matrix ref_coords, Matrix evects, List<Double> evals, String dir, String des, String type_pca, String model_pca)
			{

				this.input_coords = data;
				this.reference_PDB_coordinates = ref_coords;
				this.top_evectors = evects;
				this.eigenvalues = evals;
				this.directory = dir;
				this.description = des;
				this.type = type_pca;
				this.model = model_pca;

				this.number_of_modes = top_evectors.getColumnDimension();
				this.number_of_conformations = input_coords.getColumnDimension();
				this.ROWS = top_evectors.getRowDimension();
				this.number_of_atoms = (ROWS / 3);
				this.delta_vector_series = new Matrix(ROWS, number_of_conformations);

				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + type + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();
			}

		/**
		 * Computes the Cartesian delta vector series and writes it to file.
		 */
		public void get_Cartesian_DV_Series()
			{
				Matrix ref_col = reference_PDB_coordinates;
				Matrix col = new Matrix(ROWS, 1);
				for (int b = 0; b < number_of_conformations; b++)
					{
						col = input_coords.getMatrix(0, ROWS - 1, b, b);
						Matrix delta = col.minus(ref_col);
						delta_vector_series.setMatrix(0, ROWS - 1, b, b, delta);
					}
				// input_coords = null;
				System.gc();

				path = out_dir + "ss_" + number_of_atoms + "_delta_vectors.txt";
				Matrix_IO.write_Matrix(delta_vector_series, path, 9, 3);

				get_DVPs(delta_vector_series);
			}

		/**
		 * Computes the Cartesian delta vector projections and writes them to files: Un-normed, normed, weighted, weighted-normed
		 */
		private void get_DVPs(Matrix dv_series)
			{
				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + type + File.separatorChar + model + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();

				file_name_head = out_dir + "ss_" + number_of_atoms;

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
								Matrix data2 = dv_series.getMatrix(row_index_1, row_index_2, inner, inner);
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
				path = file_name_head + "_top_" + number_of_modes + "_DVPs_" + type + ".txt";
				Matrix_IO.write_Matrix(projections, path, 9, 3);
				// projections = null;
				path = file_name_head + "_top_" + number_of_modes + "_normed_DVPs_" + type + ".txt";
				Matrix_IO.write_Matrix(normed_projections, path, 9, 3);
				// normed_projections = null;
				path = file_name_head + "_top_" + number_of_modes + "_weighted_DVPs_" + type + ".txt";
				Matrix_IO.write_Matrix(weighted_projections, path, 9, 3);
				// weighted_projections = null;
				path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_DVPs_" + type + ".txt";
				Matrix_IO.write_Matrix(weighted_normed_projections, path, 9, 3);
				// weighted_normed_projections = null;

				System.gc();
			}

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

		// ********************************** SETTERS ***********************************************************

		public void set_Out_dir(String out_dir)
			{
				this.out_dir = out_dir;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();
			}
	}
