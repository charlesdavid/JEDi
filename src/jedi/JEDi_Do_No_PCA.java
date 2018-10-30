package jedi;

import java.util.List;
import java.util.Vector;

import Jama.Matrix;

public class JEDi_Do_No_PCA
	{

		public int number_of_atoms, number_of_residues;
		public String directory, description;
		public double z_cutoff;
		public List<Double> conformation_rmsds, residue_rmsd_list, Z_Scores;
		public Matrix original_reference_coordinates, transformed_reference_coordinates, original_PDB_coordinates, transformed_subset_PDB_coordinates,
				adjusted_PDB_coordinates_ROWS, conf_Z_scores;
		public Vector<Atom> atoms;

		/* ******************************************* CONSTRUCTOR ******************************************************************************** */

		public JEDi_Do_No_PCA(String dir, String desc, Matrix coordinates, Matrix ref_coords)
			{
				super();
				this.directory = dir;
				this.description = desc;
				this.original_PDB_coordinates = coordinates;
				this.original_reference_coordinates = ref_coords;
			}

		// ********************************************** SETTERS ********************************************************************************* */

		public void set_z_cutoff(double z)
			{
				z_cutoff = z;
			}

		// ********************************************** METHODS ********************************************************************************* */

		public void get_Transformed_Coords()
			{
				JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_coordinates, original_reference_coordinates, directory, description);
					{
						tf_coords.set_z_cutoff(z_cutoff);

						transformed_reference_coordinates = tf_coords.get_Transformed_reference_coordinates();
						transformed_subset_PDB_coordinates = tf_coords.get_SS_Transformed_coords();
						adjusted_PDB_coordinates_ROWS = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
						conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
						residue_rmsd_list = tf_coords.get_SS_RMSF();
						conf_Z_scores = tf_coords.get_conf_Z_scores();
						number_of_atoms = tf_coords.getNumber_of_atoms();
					}
			}

		// ********************************************** GETTERS *********************************************************************************** */

		/**
		 * @return the number_of_residues
		 */
		public int getNumber_of_residues()
			{
				return number_of_residues;
			}

		/**
		 * @return the directory
		 */
		public String getDirectory()
			{
				return directory;
			}

		/**
		 * @return the description
		 */
		public String getDescription()
			{
				return description;
			}

		/**
		 * @return the z_cutoff
		 */
		public double getZ_cutoff()
			{
				return z_cutoff;
			}

		/**
		 * @return the conformation_rmsds
		 */
		public List<Double> getConformation_rmsds()
			{
				return conformation_rmsds;
			}

		/**
		 * @return the residue_rmsd_list
		 */
		public List<Double> getResidue_rmsd_list()
			{
				return residue_rmsd_list;
			}

		/**
		 * @return the z_Scores
		 */
		public List<Double> getZ_Scores()
			{
				return Z_Scores;
			}

		/**
		 * @return the original_PDB_coordinates
		 */
		public Matrix getOriginal_PDB_coordinates()
			{
				return original_PDB_coordinates;
			}

		/**
		 * @return the transformed_subset_PDB_coordinates
		 */
		public Matrix getTransformed_subset_PDB_coordinates()
			{
				return transformed_subset_PDB_coordinates;
			}

		/**
		 * @return the adjusted_PDB_coordinates_ROWS
		 */
		public Matrix getAdjusted_PDB_coordinates_ROWS()
			{
				return adjusted_PDB_coordinates_ROWS;
			}

		/**
		 * @return the conf_Z_scores
		 */
		public Matrix getConf_Z_scores()
			{
				return conf_Z_scores;
			}

		/**
		 * @return the atoms
		 */
		public Vector<Atom> getAtoms()
			{
				return atoms;
			}
	}
