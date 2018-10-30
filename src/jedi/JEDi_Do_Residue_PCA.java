package jedi;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Do_Cartesian: Top class for implementing the All Atom Residue Cartesian analysis. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Do_Residue_PCA
	{

		String directory, out_dir, description, res_index;
		final String Q = "COV", R = "CORR", P = "PCORR";
		int offset, number_of_residues, number_Of_Atoms_In_Residue, residue_Index, residue_Index_adj, number_of_modes_Residue, number_of_conformations;
		double z_cutoff, trace_COV, trace_CORR, trace_PCORR, cond_COV, cond_CORR, cond_PCORR, det_COV, det_CORR, det_PCORR, rank_COV, rank_CORR, rank_PCORR;
		List<Integer> residue_list, numbers_of_atoms_per_residue;
		List<Double> transformed_conformation_rmsds, transformed_residue_rmsd_list, top_cartesian_eigenvalues_COV, top_cartesian_eigenvalues_CORR, top_cartesian_eigenvalues_PCORR,
				eigenvalues_COV, eigenvalues_CORR, eigenvalues_PCORR;
		double[] pca_mode_max_COV, pca_mode_min_COV, pca_mode_max_CORR, pca_mode_min_CORR, pca_mode_max_PCORR, pca_mode_min_PCORR;
		Matrix original_reference_PDB_coordinates, transformed_reference_PDB_coordinates, original_PDB_coordinates, transformed_PDB_coordinates, residue_PDB_coordinates,
				residue_REF_PDB_coordinates, cov, corr, conf_z_scores, var_z_scores;
		Matrix original_Residue_reference_PDB_coordinates, transformed_Residue_reference_coordinates, transformed_Residue_PDB_coordinates, adjusted_Residue_PDB_coordinates_rows,
				delta_vectors, mean_centered_data;
		Matrix top_cartesian_evectors_COV, top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, square_pca_modes_COV, weighted_square_pca_modes_COV, weighted_pca_modes_COV,
				square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR, pca_modes_COV, pca_modes_CORR, square_pca_modes_PCORR,
				weighted_square_pca_modes_PCORR, weighted_pca_modes_PCORR, pca_modes_PCORR;
		Matrix normed_projections_COV, normed_projections_CORR, projections_COV, projections_CORR, normed_projections_PCORR, projections_PCORR, weighted_projections_COV,
				weighted_projections_CORR, weighted_normed_projections_COV, weighted_normed_projections_CORR;
		Matrix trimmmed_PDB_coordinates_COLS, adjusted_PDB_coordinates_rows, Residue_Generalized_Coordinates_COV, Residue_Generalized_Coordinates_CORR,
				Residue_Generalized_Coordinates_PCORR;
		List<Matrix> Residue_Eigenvectors_COV, Residue_Delta_Vectors, Residues_Centered_Data;
		List<List<Double>> Residue_Eigenvalues_COV;
		boolean exist, success;

		/* ********************************************************* CONSTRUCTOR **************************************************************** */

		public JEDi_Do_Residue_PCA(String dir, String desc, List<Integer> residues, List<Integer> residue_atoms, int number_of_modes_Res, Matrix PDB_coordinates,
				Matrix Ref_PDB_coords)
			{

				super();

				this.directory = dir;
				this.description = desc;
				this.residue_list = residues;
				this.number_of_residues = residue_list.size();
				this.numbers_of_atoms_per_residue = residue_atoms;
				this.number_of_modes_Residue = number_of_modes_Res;
				this.original_PDB_coordinates = PDB_coordinates;
				this.original_reference_PDB_coordinates = Ref_PDB_coords;

				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Residue_cPCA" + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();

				this.number_of_conformations = original_PDB_coordinates.getColumnDimension();

				Residue_Generalized_Coordinates_COV = new Matrix(residue_list.size() * number_of_modes_Residue, number_of_conformations);
				Residue_Delta_Vectors = new ArrayList<Matrix>();
				Residue_Eigenvectors_COV = new ArrayList<Matrix>();
				Residues_Centered_Data = new ArrayList<Matrix>();
				Residue_Eigenvalues_COV = new ArrayList<List<Double>>();
			}

		/* ************************************** PUBLIC METHOD ******************************************************* */

		public void do_Cartesian_Residue_Local() // Uses local alignment of each residue
			{
				offset = 0;
				System.out.println("Residue list size " + residue_list.size());

				for (int i = 0; i < residue_list.size(); i++)
					{

						residue_Index = (i);
						residue_Index_adj = residue_list.get(i);
						res_index = String.format("%03d", (residue_Index_adj + 1));
						number_Of_Atoms_In_Residue = numbers_of_atoms_per_residue.get(residue_Index_adj);
						System.out.println("Residue " + res_index + "     Number_Of_Atoms_In_Residue " + number_Of_Atoms_In_Residue);

						JEDi_Get_Residue_Coords ss = new JEDi_Get_Residue_Coords(directory, description, original_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
							{
								residue_PDB_coordinates = ss.get_residue_Coords();
								String name = "Residue_" + res_index + "_PDB_coordinates.txt";
								Matrix_IO.write_Matrix(residue_PDB_coordinates, out_dir, name, 9, 3);
							}

						JEDi_Get_Residue_Coords ref_ss = new JEDi_Get_Residue_Coords(directory, description, original_reference_PDB_coordinates, number_Of_Atoms_In_Residue,
								offset);
							{
								residue_REF_PDB_coordinates = ref_ss.get_residue_Coords();
								String name = "Residue_" + res_index + "_reference_PDB_coordinates.txt";
								Matrix_IO.write_Matrix(residue_REF_PDB_coordinates, out_dir, name, 9, 3);
							}
						System.out.println("\nProcessing Residue  " + res_index);
						get_Local_Transformed_Coords();
						get_Cartesian_Residue_PCA();
						get_Cartesian_Residue_DVPs();
						get_Cartesian_Residue_PCs();

						offset += number_Of_Atoms_In_Residue;
					}
			}

		public void do_Cartesian_Residue_Global() // Uses global alignment of entire protein or selected subset
			{
				offset = 0;
				// System.out.println("Residue list size " + residue_list.size() + "\n");
				// System.out.println("Doing global alignment of all residues:\n");
				get_Global_Transformed_Coords();

				for (int i = 0; i < residue_list.size(); i++)
					{

						residue_Index = (i);
						residue_Index_adj = residue_list.get(i);
						res_index = String.format("%03d", (residue_Index_adj + 1));
						number_Of_Atoms_In_Residue = numbers_of_atoms_per_residue.get(residue_Index_adj);

						// System.out.println("\tResidue " + res_index + " Number_Of_Atoms_In_Residue " + number_Of_Atoms_In_Residue);

						JEDi_Get_Residue_Coords ss = new JEDi_Get_Residue_Coords(directory, description, adjusted_PDB_coordinates_rows, number_Of_Atoms_In_Residue, offset);
							{
								residue_PDB_coordinates = ss.get_residue_Coords();
								String name = "Residue_" + res_index + "_PDB_coordinates.txt";
								Matrix_IO.write_Matrix(residue_PDB_coordinates, out_dir, name, 9, 3);
							}

						JEDi_Get_Residue_Coords ref_ss = new JEDi_Get_Residue_Coords(directory, description, transformed_reference_PDB_coordinates, number_Of_Atoms_In_Residue,
								offset);
							{
								residue_REF_PDB_coordinates = ref_ss.get_residue_Coords();
								String name = "Residue_" + res_index + "_reference_PDB_coordinates.txt";
								Matrix_IO.write_Matrix(residue_REF_PDB_coordinates, out_dir, name, 9, 3);
							}

						System.out.println("\nProcessing Residue  " + res_index);
						get_Cartesian_Residue_PCA();
						get_Cartesian_Residue_DVPs();
						get_Cartesian_Residue_PCs();

						offset += number_Of_Atoms_In_Residue;
					}
			}

		/* ********************************************************** PRIVATE METHODS *********************************************************************** */

		private void get_Local_Transformed_Coords()
			{
				JEDi_Get_Transformed_Residue_Coordinates tf_coords = new JEDi_Get_Transformed_Residue_Coordinates(residue_PDB_coordinates, residue_REF_PDB_coordinates, directory,
						description);
					{
						tf_coords.setRes_index(res_index);
						tf_coords.set_Output_Directory(out_dir);
						tf_coords.set_z_cutoff(z_cutoff);
						original_Residue_reference_PDB_coordinates = tf_coords.get_Original_reference_coordinates();
						transformed_Residue_reference_coordinates = tf_coords.get_Transformed_reference_coordinates();
						transformed_Residue_PDB_coordinates = tf_coords.get_SS_Transformed_coords();
						transformed_conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
						transformed_residue_rmsd_list = tf_coords.get_SS_Residue_RMSFs();
						conf_z_scores = tf_coords.get_conf_Z_scores();
						adjusted_Residue_PDB_coordinates_rows = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
						System.gc();
					}
			}

		private void get_Global_Transformed_Coords()
			{
				JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_coordinates, original_reference_PDB_coordinates, directory,
						description);
					{
						tf_coords.set_Output_Directory(out_dir);
						tf_coords.set_z_cutoff(z_cutoff);

						transformed_reference_PDB_coordinates = tf_coords.get_Transformed_reference_coordinates();
						transformed_PDB_coordinates = tf_coords.get_SS_Transformed_coords();
						transformed_conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
						transformed_residue_rmsd_list = tf_coords.get_SS_RMSF();
						conf_z_scores = tf_coords.get_conf_Z_scores();
						adjusted_PDB_coordinates_rows = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
						System.gc();
					}
			}

		private void get_Cartesian_Residue_PCA()
			{
				JEDi_Get_Residue_PCA cr_pca = new JEDi_Get_Residue_PCA(residue_PDB_coordinates, number_of_modes_Residue, directory, description);
					{
						cr_pca.set_Residue_Index(residue_Index_adj);
						cr_pca.get_Cartesian_PCA();
						mean_centered_data = cr_pca.getCentered_input_data();
						// COVARIANCE METHOD
						cond_COV = cr_pca.get_cond_COV();
						trace_COV = cr_pca.get_trace_COV();
						det_COV = cr_pca.get_det_COV();
						rank_COV = cr_pca.get_rank_COV();
						eigenvalues_COV = cr_pca.getEigenvalues_COV();
						top_cartesian_eigenvalues_COV = cr_pca.getTop_eigenvalues_COV();
						top_cartesian_evectors_COV = cr_pca.getTop_evectors_COV();
						pca_modes_COV = cr_pca.getPca_modes_COV();
						square_pca_modes_COV = cr_pca.getSquare_pca_modes_COV();
						weighted_square_pca_modes_COV = cr_pca.getWeighted_square_pca_modes_COV();
						weighted_pca_modes_COV = cr_pca.getWeighted_pca_modes_COV();
						pca_mode_min_COV = cr_pca.get_pca_mode_COV_min();
						pca_mode_max_COV = cr_pca.get_pca_mode_COV_max();
						// ******************************************************************* // For the hierarchical PCA
						Residue_Eigenvectors_COV.add(top_cartesian_evectors_COV);
						Residue_Eigenvalues_COV.add(top_cartesian_eigenvalues_COV);
						// ******************************************************************* // For the hierarchical PCA
						// CORRELATION METHOD
						trace_CORR = cr_pca.get_trace_CORR();
						eigenvalues_CORR = cr_pca.getEigenvalues_CORR();
						top_cartesian_evectors_CORR = cr_pca.getTop_evectors_CORR();
						top_cartesian_eigenvalues_CORR = cr_pca.getEigenvalues_CORR();
						pca_modes_CORR = cr_pca.getPca_modes_CORR();
						square_pca_modes_CORR = cr_pca.getSquare_pca_modes_CORR();
						weighted_square_pca_modes_CORR = cr_pca.getWeighted_square_pca_modes_CORR();
						weighted_pca_modes_CORR = cr_pca.getWeighted_pca_modes_CORR();
						pca_mode_min_CORR = cr_pca.get_pca_mode_CORR_min();
						pca_mode_max_CORR = cr_pca.get_pca_mode_CORR_max();
						// PARTIAL CORRELATION METHOD
						trace_PCORR = cr_pca.get_trace_PCORR();
						eigenvalues_PCORR = cr_pca.getEigenvalues_PCORR();
						top_cartesian_evectors_PCORR = cr_pca.getTop_evectors_PCORR();
						top_cartesian_eigenvalues_PCORR = cr_pca.getEigenvalues_PCORR();
						pca_modes_PCORR = cr_pca.getPca_modes_PCORR();
						square_pca_modes_PCORR = cr_pca.getSquare_pca_modes_PCORR();
						weighted_square_pca_modes_PCORR = cr_pca.getWeighted_square_pca_modes_PCORR();
						weighted_pca_modes_PCORR = cr_pca.getWeighted_pca_modes_PCORR();
						pca_mode_min_PCORR = cr_pca.get_pca_mode_PCORR_min();
						pca_mode_max_PCORR = cr_pca.get_pca_mode_PCORR_max();
					}
			}

		private void get_Cartesian_Residue_DVPs()
			{
				JEDi_Get_Residue_DVPs dvps_cov = new JEDi_Get_Residue_DVPs(residue_PDB_coordinates, residue_REF_PDB_coordinates, top_cartesian_evectors_COV, eigenvalues_COV,
						directory, description, Q);
					{
						dvps_cov.set_Residue_Index(residue_Index_adj);
						delta_vectors = dvps_cov.get_Cartesian_DV_Series();
						Residue_Delta_Vectors.add(delta_vectors);
						projections_COV = dvps_cov.getProjections();
						normed_projections_COV = dvps_cov.getNormed_projections();
						weighted_projections_COV = dvps_cov.getWeighted_projections();
						weighted_normed_projections_COV = dvps_cov.getWeighted_normed_projections();
					}
				System.gc();
			}

		private void get_Cartesian_Residue_PCs()
			{
				JEDi_Get_Residue_PCs pcs_cov = new JEDi_Get_Residue_PCs(mean_centered_data, top_cartesian_evectors_COV, eigenvalues_COV, directory, description, Q);
					{
						pcs_cov.set_Residue_Index(residue_Index_adj);
						pcs_cov.get_PCs();
						Residues_Centered_Data.add(mean_centered_data);
						projections_COV = pcs_cov.getProjections();
						normed_projections_COV = pcs_cov.getNormed_projections();
						weighted_projections_COV = pcs_cov.getWeighted_projections();
						weighted_normed_projections_COV = pcs_cov.getWeighted_normed_projections();
						// ********************************************************************************************************************************** */
						int row_index1 = residue_Index * number_of_modes_Residue;
						int row_index2 = row_index1 + number_of_modes_Residue - 1;
						int col_index1 = 0;
						int col_index2 = Residue_Generalized_Coordinates_COV.getColumnDimension() - 1;
						// System.out.println("row_index1 " + row_index1 + " row_index2 " + row_index2 + " col_index1 " + col_index1 + " col_index2 " +
						// col_index2); //////
						Residue_Generalized_Coordinates_COV.setMatrix(row_index1, row_index2, col_index1, col_index2, projections_COV.transpose());
						// ********************************************************************************************************************************** */
					}
				System.gc();
			}

		/* ************************************************************ SETTERS ***************************************************************************** */

		public void set_z_cutoff(double z)
			{
				z_cutoff = z;
			}

		/* ************************************************************ GETTERS ***************************************************************************** */

		public double get_z_cut()
			{
				return z_cutoff;
			}

		public double get_Trace_COV()
			{
				return trace_COV;
			}

		public double get_Trace_CORR()
			{
				return trace_CORR;
			}

		public double get_cond_COV()
			{
				return cond_COV;
			}

		public double get_cond_CORR()
			{
				return cond_CORR;
			}

		public double get_Trace_PCORR()
			{
				return trace_PCORR;
			}

		public double get_Cond_COV()
			{
				return cond_COV;
			}

		public double get_Cond_CORR()
			{
				return cond_CORR;
			}

		public double get_Cond_PCORR()
			{
				return cond_PCORR;
			}

		public double get_Det_COV()
			{
				return det_COV;
			}

		public double get_Det_CORR()
			{
				return det_CORR;
			}

		public double get_Det_PCORR()
			{
				return det_PCORR;
			}

		public double get_Rank_COV()
			{
				return rank_COV;
			}

		public double get_Rank_CORR()
			{
				return rank_CORR;
			}

		public double get_Rank_PCORR()
			{
				return rank_PCORR;
			}

		public List<Double> getTop_cartesian_eigenvalues_COV()
			{
				return top_cartesian_eigenvalues_COV;
			}

		public List<Double> getTop_cartesian_eigenvalues_CORR()
			{
				return top_cartesian_eigenvalues_CORR;
			}

		public double[] getPca_mode_max_COV()
			{
				return pca_mode_max_COV;
			}

		public double[] getPca_mode_min_COV()
			{
				return pca_mode_min_COV;
			}

		public double[] getPca_mode_max_CORR()
			{
				return pca_mode_max_CORR;
			}

		public double[] getPca_mode_min_CORR()
			{
				return pca_mode_min_CORR;
			}

		public Matrix getCov()
			{
				return cov;
			}

		public Matrix getTop_cartesian_evectors_COV()
			{
				return top_cartesian_evectors_COV;
			}

		public Matrix getWeighted_pca_modes_COV()
			{
				return weighted_pca_modes_COV;
			}

		public Matrix getSquare_pca_modes_COV()
			{
				return square_pca_modes_COV;
			}

		public Matrix getWeighted_square_pca_modes_COV()
			{
				return weighted_square_pca_modes_COV;
			}

		public Matrix getNormed_projections_COV()
			{
				return normed_projections_COV;
			}

		public Matrix getProjections_COV()
			{
				return projections_COV;
			}

		public Matrix getCorr()
			{
				return corr;
			}

		public Matrix getTop_cartesian_evectors_CORR()
			{
				return top_cartesian_evectors_CORR;
			}

		public Matrix getWeighted_pca_modes_CORR()
			{
				return weighted_pca_modes_CORR;
			}

		public Matrix getSquare_pca_modes_CORR()
			{
				return square_pca_modes_CORR;
			}

		public Matrix getWeighted_square_pca_modes_CORR()
			{
				return weighted_square_pca_modes_CORR;
			}

		public Matrix getNormed_projections_CORR()
			{
				return normed_projections_CORR;
			}

		public Matrix getProjections_CORR()
			{
				return projections_CORR;
			}

		public double getZ_cutoff()
			{
				return z_cutoff;
			}

		public Matrix getNormed_projections_PCORR()
			{
				return normed_projections_PCORR;
			}

		public Matrix getProjections_PCORR()
			{
				return projections_PCORR;
			}

		public Matrix getWeighted_projections_COV()
			{
				return weighted_projections_COV;
			}

		public Matrix getWeighted_projections_CORR()
			{
				return weighted_projections_CORR;
			}

		public Matrix getWeighted_normed_projections_COV()
			{
				return weighted_normed_projections_COV;
			}

		public Matrix getWeighted_normed_projections_CORR()
			{
				return weighted_normed_projections_CORR;
			}

		public Matrix getSubset_PDB_coordinates()
			{
				return residue_PDB_coordinates;
			}

		public Matrix getTransformed_PDB_coordinates()
			{
				return transformed_PDB_coordinates;
			}

		public int getNumber_of_residues()
			{
				return number_of_residues;
			}

		public List<Double> getTransformed_conformation_rmsds()
			{
				return transformed_conformation_rmsds;
			}

		public List<Double> getTransformed_residue_rmsd_list()
			{
				return transformed_residue_rmsd_list;
			}

		public Matrix getPca_modes_COV()
			{
				return pca_modes_COV;
			}

		public Matrix getPca_modes_CORR()
			{
				return pca_modes_CORR;
			}

		public double getTrace_COV()
			{
				return trace_COV;
			}

		public double getTrace_CORR()
			{
				return trace_CORR;
			}

		public double getTrace_PCORR()
			{
				return trace_PCORR;
			}

		public double getDet_COV()
			{
				return det_COV;
			}

		public double getDet_CORR()
			{
				return det_CORR;
			}

		public double getDet_PCORR()
			{
				return det_PCORR;
			}

		public double getRank_COV()
			{
				return rank_COV;
			}

		public double getRank_CORR()
			{
				return rank_CORR;
			}

		public double getRank_PCORR()
			{
				return rank_PCORR;
			}

		public List<Double> getTop_cartesian_eigenvalues_PCORR()
			{
				return top_cartesian_eigenvalues_PCORR;
			}

		public List<Double> getEigenvalues_COV()
			{
				return eigenvalues_COV;
			}

		public List<Double> getEigenvalues_CORR()
			{
				return eigenvalues_CORR;
			}

		public List<Double> getEigenvalues_PCORR()
			{
				return eigenvalues_PCORR;
			}

		public Matrix getTop_cartesian_evectors_PCORR()
			{
				return top_cartesian_evectors_PCORR;
			}

		public double[] getPca_mode_max_PCORR()
			{
				return pca_mode_max_PCORR;
			}

		public double[] getPca_mode_min_PCORR()
			{
				return pca_mode_min_PCORR;
			}

		public Matrix getSquare_pca_modes_PCORR()
			{
				return square_pca_modes_PCORR;
			}

		public Matrix get_Residue_Generalized_Coordinates_COV()
			{
				return Residue_Generalized_Coordinates_COV;
			}

		public Matrix get_Residue_Generalized_Coordinates_CORR()
			{
				return Residue_Generalized_Coordinates_CORR;
			}

		public Matrix get_Residue_Generalized_Coordinates_PCORR()
			{
				return Residue_Generalized_Coordinates_PCORR;
			}

		public List<Matrix> get_Residue_Delta_Vectors()
			{
				return Residue_Delta_Vectors;
			}

		public List<Matrix> get_Residue_Eigenvectors_COV()
			{
				return Residue_Eigenvectors_COV;
			}

		public List<List<Double>> get_Residue_Eigenvalues_COV()
			{
				return Residue_Eigenvalues_COV;
			}

		public List<Matrix> getResidues_Centered_Data()
			{
				return Residues_Centered_Data;
			}
	}
