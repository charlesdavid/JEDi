package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Do_Cartesian: Top class for implementing the Cartesian analysis. Copyright (C) 2012 Dr. Charles David
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

public class JEDi_Do_Cartesian
	{

		public String directory, out_dir, description, type, model, Q = "COV", R = "CORR", P = "PCORR";
		public int number_of_atoms, number_of_residues, number_of_modes_SS, number_of_conformations;
		public double z_cutoff, trace_COV, trace_CORR, trace_PCORR, cond_COV, cond_CORR, cond_PCORR, det_COV, det_CORR, det_PCORR, rank_COV, rank_CORR, rank_PCORR;
		public double[] pca_mode_max_COV, pca_mode_min_COV, pca_mode_max_CORR, pca_mode_min_CORR, pca_mode_max_PCORR, pca_mode_min_PCORR;
		public List<Double> transformed_conformation_rmsds, transformed_residue_rmsd_list, top_cartesian_eigenvalues_COV, top_cartesian_eigenvalues_CORR,
				top_cartesian_eigenvalues_PCORR, eigenvalues_COV, eigenvalues_CORR, eigenvalues_PCORR;
		public Matrix transformed_reference_coordinates, subset_PDB_coordinates, subset_REF_PDB_coordinates, transformed_PDB_coordinates, cov, corr, conf_z_scores, var_z_scores;
		public Matrix top_cartesian_evectors_COV, top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, square_pca_modes_COV, weighted_square_pca_modes_COV,
				weighted_pca_modes_COV, square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR, pca_modes_COV, pca_modes_CORR, square_pca_modes_PCORR,
				weighted_square_pca_modes_PCORR, weighted_pca_modes_PCORR, pca_modes_PCORR;
		public Matrix normed_projections_COV, normed_projections_CORR, projections_COV, projections_CORR, normed_projections_PCORR, projections_PCORR, weighted_projections_COV,
				weighted_projections_CORR, weighted_normed_projections_COV, weighted_normed_projections_CORR;
		public Matrix adjusted_PDB_coordinates_rows;
		boolean exist, success;
		/* ******************************************************************* CONSTRUCTOR **************************************************************** */

		/**
		 * Initiates the Cartesian Subset Analysis:
		 *
		 * @param directory                The working directory
		 * @param                          data.description The job description
		 * @param pdb_ref                  The PDB Reference File
		 * @param rl_SS                    The Cartesian Residue List
		 * @param                          data.number_of_modes_SS The number of Cartesian PCA modes to compute.
		 * @param original_PDB_coordinates The Coordinates Matrix
		 */
		public JEDi_Do_Cartesian(String directory, String desc, int number_of_modes, Matrix PDB_coordinates, Matrix Ref_PDB_coords, String type_pca)
			{

				super();
				this.directory = directory;
				this.description = desc;
				this.number_of_modes_SS = number_of_modes;
				this.subset_PDB_coordinates = PDB_coordinates;
				this.subset_REF_PDB_coordinates = Ref_PDB_coords;
				this.type = type_pca;
				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "type" + File.separatorChar;
				this.number_of_conformations = subset_PDB_coordinates.getColumnDimension();
			}
		/* ************************************************** DRIVER METHODS ******************************************************************************** */

		public void do_Cartesian()
			{
				get_Transformed_Coords();
				get_Cartesian_PCA();
				get_Cartesian_DVPs();
				do_FES();
			}

		/* ****************************************************** METHODS *********************************************************************************** */

		private void get_Transformed_Coords()
			{
				JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates, subset_REF_PDB_coordinates, directory, description);
					{
						System.out.println("Transforming Coordinates:");
						tf_coords.set_Output_Directory(out_dir);
						tf_coords.set_z_cutoff(z_cutoff);

						transformed_reference_coordinates = tf_coords.get_Transformed_reference_coordinates();
						transformed_PDB_coordinates = tf_coords.get_SS_Transformed_coords();

						transformed_conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
						transformed_residue_rmsd_list = tf_coords.get_SS_RMSF();
						conf_z_scores = tf_coords.get_conf_Z_scores();
						adjusted_PDB_coordinates_rows = tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();

						System.out.println("Done:");
						System.gc();
					}
			}

		private void get_Cartesian_PCA()
			{
				JEDi_Get_Cartesian_PCA c_pca = new JEDi_Get_Cartesian_PCA(adjusted_PDB_coordinates_rows, number_of_modes_SS, directory, description, type);
					{
						System.out.println("Doing " + type + " PCA:");
						c_pca.get_Cartesian_PCA();
						System.out.println("Done:");

						// COVARIANCE METHOD
						cond_COV = c_pca.get_cond_COV();
						trace_COV = c_pca.get_trace_COV();
						det_COV = c_pca.get_det_COV();
						rank_COV = c_pca.get_rank_COV();
						eigenvalues_COV = c_pca.getEigenvalues_COV();
						top_cartesian_evectors_COV = c_pca.getTop_evectors_COV();
						top_cartesian_eigenvalues_COV = c_pca.getTop_eigenvalues_COV();
						pca_modes_COV = c_pca.getPca_modes_COV();
						square_pca_modes_COV = c_pca.getSquare_pca_modes_COV();
						weighted_square_pca_modes_COV = c_pca.getWeighted_square_pca_modes_COV();
						weighted_pca_modes_COV = c_pca.getWeighted_pca_modes_COV();
						pca_mode_min_COV = c_pca.get_pca_mode_COV_min();
						pca_mode_max_COV = c_pca.get_pca_mode_COV_max();
						// CORRELATION METHOD
						trace_CORR = c_pca.get_trace_CORR();
						eigenvalues_CORR = c_pca.getEigenvalues_CORR();
						top_cartesian_evectors_CORR = c_pca.getTop_evectors_CORR();
						top_cartesian_eigenvalues_CORR = c_pca.getTop_eigenvalues_CORR();
						pca_modes_CORR = c_pca.getPca_modes_CORR();
						square_pca_modes_CORR = c_pca.getSquare_pca_modes_CORR();
						weighted_square_pca_modes_CORR = c_pca.getWeighted_square_pca_modes_CORR();
						weighted_pca_modes_CORR = c_pca.getWeighted_pca_modes_CORR();
						pca_mode_min_CORR = c_pca.get_pca_mode_CORR_min();
						pca_mode_max_CORR = c_pca.get_pca_mode_CORR_max();
						// PARTIAL CORRELATION METHOD
						trace_PCORR = c_pca.get_trace_PCORR();
						eigenvalues_PCORR = c_pca.getEigenvalues_PCORR();
						top_cartesian_evectors_PCORR = c_pca.getTop_evectors_PCORR();
						top_cartesian_eigenvalues_PCORR = c_pca.getEigenvalues_PCORR();
						pca_modes_PCORR = c_pca.getPca_modes_PCORR();
						square_pca_modes_PCORR = c_pca.getSquare_pca_modes_PCORR();
						weighted_square_pca_modes_PCORR = c_pca.getWeighted_square_pca_modes_PCORR();
						weighted_pca_modes_PCORR = c_pca.getWeighted_pca_modes_PCORR();
						pca_mode_min_PCORR = c_pca.get_pca_mode_PCORR_min();
						pca_mode_max_PCORR = c_pca.get_pca_mode_PCORR_max();
					}
			}

		private void get_Cartesian_DVPs()
			{
				JEDi_Get_Cartesian_DVPs pcs_cov = new JEDi_Get_Cartesian_DVPs(transformed_PDB_coordinates, transformed_reference_coordinates, top_cartesian_evectors_COV,
						eigenvalues_COV, directory, description, type, Q);
					{
						System.out.println("Doing " + type + " DVPs:");
						pcs_cov.get_Cartesian_DV_Series();
						projections_COV = pcs_cov.getProjections();
						normed_projections_COV = pcs_cov.getNormed_projections();
						weighted_projections_COV = pcs_cov.getWeighted_projections();
						weighted_normed_projections_COV = pcs_cov.getWeighted_normed_projections();

					}
				JEDi_Get_Cartesian_DVPs pcs_corr = new JEDi_Get_Cartesian_DVPs(transformed_PDB_coordinates, transformed_reference_coordinates, top_cartesian_evectors_CORR,
						eigenvalues_CORR, directory, description, type, R);
					{
						pcs_corr.get_Cartesian_DV_Series();
						projections_CORR = pcs_corr.getProjections();
						normed_projections_CORR = pcs_corr.getNormed_projections();
						weighted_projections_CORR = pcs_corr.getWeighted_projections();
						weighted_normed_projections_CORR = pcs_corr.getWeighted_normed_projections();

					}
				JEDi_Get_Cartesian_DVPs pcs_pcorr = new JEDi_Get_Cartesian_DVPs(transformed_PDB_coordinates, transformed_reference_coordinates, top_cartesian_evectors_PCORR,
						eigenvalues_PCORR, directory, description, type, P);
					{
						pcs_pcorr.get_Cartesian_DV_Series();
						projections_PCORR = pcs_pcorr.getProjections();
						normed_projections_PCORR = pcs_pcorr.getNormed_projections();
					}
				System.out.println("Done:");
				transformed_PDB_coordinates = null;
				System.gc();
			}

		private void do_FES()
			{
				System.out.println("Calculating the Free Energy Surface for the 3 PCA models... ");

				JEDi_Get_FES fes = new JEDi_Get_FES(weighted_projections_CORR, 0, 1, number_of_conformations, 0, 0, directory, description, type, "CORR");
				fes.get_FES();
				fes.write_FES_Log();

				fes = new JEDi_Get_FES(weighted_projections_COV, 0, 1, number_of_conformations, 0, 0, directory, description, type, "COV");
				fes.get_FES();
				fes.write_FES_Log();

				fes = new JEDi_Get_FES(projections_PCORR, 0, 1, number_of_conformations, 0, 0, directory, description, type, "PCORR");
				fes.get_FES();
				fes.write_FES_Log();

				System.out.println("Done.");
			}


		/* ************************************************************ SETTERS ********************************************************************* */

		/**
		 * Sets the Z cutoff for outlier processing
		 *
		 * @param z
		 */
		public void set_z_cutoff(double z)
			{
				z_cutoff = z;
			}

		public void set_Out_dir(String out_dir)
			{
				this.out_dir = out_dir;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();
			}

		/* ************************************************************** GETTERS ******************************************************************* */

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

		/**
		 * @return the z_cutoff
		 */
		public double getZ_cutoff()
			{
				return z_cutoff;
			}

		/**
		 * @return the normed_projections_PCORR
		 */
		public Matrix getNormed_projections_PCORR()
			{
				return normed_projections_PCORR;
			}

		/**
		 * @return the projections_PCORR
		 */
		public Matrix getProjections_PCORR()
			{
				return projections_PCORR;
			}

		/**
		 * @return the weighted_projections_COV
		 */
		public Matrix getWeighted_projections_COV()
			{
				return weighted_projections_COV;
			}

		/**
		 * @return the weighted_projections_CORR
		 */
		public Matrix getWeighted_projections_CORR()
			{
				return weighted_projections_CORR;
			}

		/**
		 * @return the weighted_normed_projections_COV
		 */
		public Matrix getWeighted_normed_projections_COV()
			{
				return weighted_normed_projections_COV;
			}

		/**
		 * @return the weighted_normed_projections_CORR
		 */
		public Matrix getWeighted_normed_projections_CORR()
			{
				return weighted_normed_projections_CORR;
			}

		public Matrix getSubset_PDB_coordinates()
			{

				return subset_PDB_coordinates;
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

		/**
		 * @return the trace_COV
		 */
		public double getTrace_COV()
			{
				return trace_COV;
			}

		/**
		 * @return the trace_CORR
		 */
		public double getTrace_CORR()
			{
				return trace_CORR;
			}

		/**
		 * @return the trace_PCORR
		 */
		public double getTrace_PCORR()
			{
				return trace_PCORR;
			}

		/**
		 * @return the det_COV
		 */
		public double getDet_COV()
			{
				return det_COV;
			}

		/**
		 * @return the det_CORR
		 */
		public double getDet_CORR()
			{
				return det_CORR;
			}

		/**
		 * @return the det_PCORR
		 */
		public double getDet_PCORR()
			{
				return det_PCORR;
			}

		/**
		 * @return the rank_COV
		 */
		public double getRank_COV()
			{
				return rank_COV;
			}

		/**
		 * @return the rank_CORR
		 */
		public double getRank_CORR()
			{
				return rank_CORR;
			}

		/**
		 * @return the rank_PCORR
		 */
		public double getRank_PCORR()
			{
				return rank_PCORR;
			}

		/**
		 * @return the top_cartesian_eigenvalues_PCORR
		 */
		public List<Double> getTop_cartesian_eigenvalues_PCORR()
			{
				return top_cartesian_eigenvalues_PCORR;
			}

		/**
		 * @return the eigenvalues_COV
		 */
		public List<Double> getEigenvalues_COV()
			{
				return eigenvalues_COV;
			}

		/**
		 * @return the eigenvalues_CORR
		 */
		public List<Double> getEigenvalues_CORR()
			{
				return eigenvalues_CORR;
			}

		/**
		 * @return the eigenvalues_PCORR
		 */
		public List<Double> getEigenvalues_PCORR()
			{
				return eigenvalues_PCORR;
			}

		/**
		 * @return the top_cartesian_evectors_PCORR
		 */
		public Matrix getTop_cartesian_evectors_PCORR()
			{
				return top_cartesian_evectors_PCORR;
			}

		/**
		 * @return the pca_mode_max_PCORR
		 */
		public double[] getPca_mode_max_PCORR()
			{
				return pca_mode_max_PCORR;
			}

		/**
		 * @return the pca_mode_min_PCORR
		 */
		public double[] getPca_mode_min_PCORR()
			{
				return pca_mode_min_PCORR;
			}

		/**
		 * @return the square_pca_modes_PCORR
		 */
		public Matrix getSquare_pca_modes_PCORR()
			{
				return square_pca_modes_PCORR;
			}
	}
