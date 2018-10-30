package jedi;

import java.io.File;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Transformed_Coordinates: Optimally aligns all frames of a trajectory to a reference frame using Quaternion operations. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Transformed_Residue_Coordinates
	{

		String directory, out_dir, description, res_index;
		double z_cut_off;
		int ROWS, COLS, number_of_atoms;
		List<Double> original_conformation_rmsds, transformed_conformation_rmsds, transformed_residue_rmsd_list, Z_Scores;
		List<Integer> conformation_outliers;
		Matrix subset_REF_PDB_coordinates, transformed_subset_REF_PDB_coordinates, subset_PDB_coordinates, transformed_subset_PDB_coordinates, adjusted_PDB_coordinates_rows,
				var_Z_scores, conf_Z_scores;
		NumberFormat nf;
		RoundingMode rm;
		JEDi_Transform_Coords tf_coords;

		/* ************************************** CONSTRUCTORS ******************************************************************************** */

		/**
		 * Constructor for transforming the original coordinates using quaternion operations and calculating the conformation and residue RMSDs. If specified, the coordinates
		 * outliers will be handled This class is used for Cartesian analysis, but not the distance analyses which use internal coordinates (distances)
		 * 
		 * @param data       The original PDB coordinates
		 * @param ref_coords The reference coordinates to use for the transformation
		 * @param dir        The working directory
		 * @param des        The job description
		 */
		JEDi_Get_Transformed_Residue_Coordinates(Matrix data, Matrix ref_coords, String dir, String des)
			{

				nf = NumberFormat.getInstance();
				rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);

				subset_PDB_coordinates = data;
				subset_REF_PDB_coordinates = ref_coords;
				directory = dir;
				description = des;

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Residue_cPCA" + File.separatorChar;
				transformed_residue_rmsd_list = new ArrayList<>();
				original_conformation_rmsds = new ArrayList<>();
				transformed_conformation_rmsds = new ArrayList<>();

				ROWS = subset_PDB_coordinates.getRowDimension();
				COLS = subset_PDB_coordinates.getColumnDimension();
				number_of_atoms = (ROWS / 3);
				transformed_subset_PDB_coordinates = new Matrix(ROWS, COLS);

				tf_coords = new JEDi_Transform_Coords(subset_REF_PDB_coordinates, subset_REF_PDB_coordinates);
				transformed_subset_REF_PDB_coordinates = tf_coords.get_transformed_coords();
			}

		/* ************************************** SETTERS ******************************************************************************** */

		/**
		 * Sets the Z cutoff for identifying outliers
		 * 
		 * @param z The Z-cutoff
		 */
		public void set_z_cutoff(double z)
			{

				z_cut_off = z;
			}

		/**
		 * Sets the output directory
		 * 
		 * @param dir The directory to use for output
		 */
		public void set_Output_Directory(String dir)
			{

				out_dir = dir;
			}

		/**
		 * Sets the residue index for base 1 numbering
		 * 
		 * @param res_index The res_index to use for writing output
		 */
		public void setRes_index(String res_index)
			{
				this.res_index = res_index;
			}

		/* ************************************** METHODS ******************************************************************************** */

		/**
		 * @return The Transformed coordinates for the specified subset
		 */
		public Matrix get_SS_Transformed_coords()
			{
				for (int z = 0; z < COLS; z++)
					{
						Matrix fc = subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);
						tf_coords = new JEDi_Transform_Coords(transformed_subset_REF_PDB_coordinates, fc);
						Matrix transformed_col_vector = tf_coords.get_transformed_coords();
						transformed_subset_PDB_coordinates.setMatrix(0, ROWS - 1, z, z, transformed_col_vector);
						tf_coords = null;
						if (z % 10 == 0)
							System.gc();
					}
				String path = out_dir + "Residue_" + res_index + "_transformed_PDB_coordinates.txt";
				Matrix_IO.write_Matrix(transformed_subset_PDB_coordinates, path, 9, 3);

				path = out_dir + "Residue_" + res_index + "_transformed_reference_PDB_coordinates.txt";
				Matrix_IO.write_Matrix(transformed_subset_REF_PDB_coordinates, path, 9, 3);

				return transformed_subset_PDB_coordinates;
			}

		/**
		 * @return The transformed conformation RMSDs
		 */
		public List<Double> get_SS_Conformation_RMSDs()
			{
				for (int z = 0; z < COLS; z++)
					{
						Matrix fc_O = subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);
						Matrix fc_T = transformed_subset_PDB_coordinates.getMatrix(0, ROWS - 1, z, z);

						JEDi_Get_Conformation_RMSD gRMSD_O = new JEDi_Get_Conformation_RMSD(subset_REF_PDB_coordinates, fc_O);
						JEDi_Get_Conformation_RMSD gRMSD_T = new JEDi_Get_Conformation_RMSD(transformed_subset_REF_PDB_coordinates, fc_T);

						double rmsd_O = gRMSD_O.get_RMSD();
						original_conformation_rmsds.add(rmsd_O);
						double rmsd_T = gRMSD_T.get_RMSD();
						transformed_conformation_rmsds.add(rmsd_T);
					}
				String path = out_dir + "Residue_" + res_index + "_original_conformation_rmsds.txt";
				List_IO.write_Double_List(original_conformation_rmsds, path, 3);
				path = out_dir + "Residue_" + res_index + "_transformed_conformation_rmsds.txt";
				List_IO.write_Double_List(transformed_conformation_rmsds, path, 3);

				return transformed_conformation_rmsds;
			}


		/**
		 * @return The PDB coordinates matrix with the outliers adjusted to their mean values.
		 */
		public Matrix get_SS_transformed_coordinates_adjusted_ROWS()
			{

				if (z_cut_off > 0)
					{
						Adjust_Outliers_by_Z_Score adr = new Adjust_Outliers_by_Z_Score(transformed_subset_PDB_coordinates);
						adr.set_Z_threshold(z_cut_off);
						adr.adjust_row_data();
						adjusted_PDB_coordinates_rows = adr.get_coorinates_adjusted();
						String path = out_dir + "Residue_" + res_index + "_Z_threshold_" + z_cut_off + "_adjusted_PDB_coordinates_ROWS.txt";
						Matrix_IO.write_Matrix(adjusted_PDB_coordinates_rows, path, 9, 3);
						Matrix var_counts = adr.get_counts();
						path = out_dir + "Residue_" + res_index + "_adjustments_per_variable.txt";
						Matrix_IO.write_Matrix(var_counts, path, 6, 0);
					} else
					{
						adjusted_PDB_coordinates_rows = transformed_subset_PDB_coordinates;
					}
				return adjusted_PDB_coordinates_rows;
			}

		/**
		 * @return The Residue RMSDs (RMSFs)
		 */
		public List<Double> get_SS_Residue_RMSFs()
			{

				JEDi_Get_Residue_RMSFs rrmsd = new JEDi_Get_Residue_RMSFs(transformed_subset_PDB_coordinates);
				transformed_residue_rmsd_list = rrmsd.get_residue_rmsfs();
				var_Z_scores = rrmsd.get_z_scores();

				String path = out_dir + "Residue_" + res_index + "_residue_atom_rmsf.txt";
				List_IO.write_Double_List(transformed_residue_rmsd_list, path, 3);

				path = out_dir + "Residue_" + res_index + "_Variable_Z_Scores.txt";
				Matrix_IO.write_Matrix(var_Z_scores, path, 6, 1);

				return transformed_residue_rmsd_list;
			}

		/* ************************************** GETTERS ******************************************************************************** */

		public Matrix get_Original_reference_coordinates()
			{
				return subset_REF_PDB_coordinates;
			}

		public Matrix get_Transformed_reference_coordinates()
			{
				return transformed_subset_REF_PDB_coordinates;
			}

		public List<Integer> get_conformation_outliers()
			{

				return conformation_outliers;
			}

		public Matrix get_conf_Z_scores()
			{

				return conf_Z_scores;
			}

		public int getNumber_of_residues()
			{

				return number_of_atoms;
			}

		public double get_z_cut()
			{

				return z_cut_off;
			}
	}
