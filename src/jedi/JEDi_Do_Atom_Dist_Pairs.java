package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;

import Jama.Matrix;

public class JEDi_Do_Atom_Dist_Pairs
	{

		int number_of_Chains, number_of_pairs, number_of_modes_dist_pairs, number_of_conformations;
		double trace_dist_COV, trace_dist_CORR, trace_dist_PCORR, cond_cov, cond_corr, cond_pcorr, det_cov, det_corr, det_pcorr, rank_cov, rank_corr, rank_pcorr, z_cutoff;
		String directory, out_dir, description, path, line, type, model;
		List<Integer> atom_list_dist1, atom_list_dist_orig1, atom_list_dist2, atom_list_dist_orig2;
		List<String> chain_idents1, chain_idents2;
		List<Double> top_distance_eigenvalues_COV, top_distance_eigenvalues_CORR, top_distance_eigenvalues_PCORR;
		double[] residue_distance_means, residue_distance_std_devs, dist_pca_mode_COV_min, dist_pca_mode_COV_max, dist_pca_mode_CORR_min, dist_pca_mode_CORR_max;
		Matrix reference_pdb_coordinates, reference_distances, distance_matrix, distance_matrix_cond, original_PDB_coordinates, subset_PDB_coordinates_dist,
				transformed_subset_PDB_coordinates_dist, cov_dist, corr_dist, top_distance_evectors_COV, top_distance_evectors_CORR, top_distance_evectors_PCORR, Z_scores, counts;
		Matrix normed_projections_dist_COV, projections_dist_COV, projections_dist_CORR, normed_projections_dist_CORR, weighted_normed_projections_dist_COV,
				weighted_normed_projections_dist_CORR, weighted_projections_dist_COV, weighted_projections_dist_CORR, normed_projections_dist_PCORR, projections_dist_PCORR;
		NumberFormat nf;
		RoundingMode rm;
		boolean exist, success;

		/* ***************************************************** CONSTRUCTOR ********************************************************************** */

		public JEDi_Do_Atom_Dist_Pairs(String dir, String desc, int modes_dist_pairs, Matrix PDB_coordinates, Matrix ref_coords, List<Integer> atm_list1,
				List<Integer> atm_list1_orig, List<Integer> atm_list2, List<Integer> atm_list2_orig, List<String> chain_ids1, List<String> chain_ids2)
			{
				super();
				this.directory = dir;
				this.description = desc;
				this.number_of_modes_dist_pairs = modes_dist_pairs;
				this.original_PDB_coordinates = PDB_coordinates;
				this.reference_pdb_coordinates = ref_coords;
				this.atom_list_dist1 = atm_list1;
				this.atom_list_dist_orig1 = atm_list1_orig;
				this.atom_list_dist2 = atm_list2;
				this.atom_list_dist_orig2 = atm_list2_orig;
				this.chain_idents1 = chain_ids1;
				this.chain_idents2 = chain_ids2;

				this.number_of_pairs = atom_list_dist1.size();
				this.number_of_conformations = original_PDB_coordinates.getColumnDimension();
				this.type = "dpPCA";

				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + type + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist) success = (new File(out_dir)).mkdirs();

				nf = NumberFormat.getInstance();
				rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
			}

		/* ************************************************** DRIVER METHOD ************************************************************************* ******* */

		/**
		 * Method for performing dPCA, processing Single Chain PDBs with no chain IDs.
		 */
		public void do_Dist()
			{
				get_Reference_Atom_Pair_Distances();
				get_Atom_Pair_Distances();
				remove_Outliers_Z_Score();
				get_Distance_Pair_PCA();
				get_Distance_Pair_DVPs();
				write_Distance_Pair_Stats();
				do_FES();
			}

		/**
		 * Method for performing dPCA, processing Multi chain PDBs with chain IDs.
		 */
		public void do_Dist_Multi()
			{
				get_Reference_Atom_Pair_Distances();
				get_Atom_Pair_Distances();
				remove_Outliers_Z_Score();
				get_Distance_Pair_PCA();
				get_Distance_Pair_DVPs();
				write_Distance_Pair_Stats_Multi();
				do_FES();
			}

		/* ******************************************************** METHODS ************************************************************************* ******* */

		private void get_Distance_Pair_DVPs()
			{
				JEDi_Get_Distance_Pair_DVPs SS_dvp_COV = new JEDi_Get_Distance_Pair_DVPs(distance_matrix, reference_distances, top_distance_evectors_COV,
						top_distance_eigenvalues_COV, directory, description, "COV", number_of_modes_dist_pairs);
					{
						// COVARIANCE METHOD
						SS_dvp_COV.get_Distance_Pair_DV_Series();
						projections_dist_COV = SS_dvp_COV.getProjections();
						normed_projections_dist_COV = SS_dvp_COV.getNormed_projections();
						weighted_projections_dist_COV = SS_dvp_COV.getWeighted_projections();
						weighted_normed_projections_dist_COV = SS_dvp_COV.getWeighted_normed_projections();

					}
				JEDi_Get_Distance_Pair_DVPs SS_dvp_CORR = new JEDi_Get_Distance_Pair_DVPs(distance_matrix, reference_distances, top_distance_evectors_CORR,
						top_distance_eigenvalues_CORR, directory, description, "CORR", number_of_modes_dist_pairs);
					{
						// CORRELATION METHOD
						SS_dvp_CORR.get_Distance_Pair_DV_Series();
						projections_dist_CORR = SS_dvp_CORR.getProjections();
						normed_projections_dist_CORR = SS_dvp_CORR.getNormed_projections();
						weighted_projections_dist_CORR = SS_dvp_CORR.getWeighted_projections();
						weighted_normed_projections_dist_CORR = SS_dvp_CORR.getWeighted_normed_projections();
					}
				JEDi_Get_Distance_Pair_DVPs SS_dvp_PCORR = new JEDi_Get_Distance_Pair_DVPs(distance_matrix, reference_distances, top_distance_evectors_CORR,
						top_distance_eigenvalues_CORR, directory, description, "PCORR", number_of_modes_dist_pairs);
					{
						// PARTIAL CORRELATION METHOD
						SS_dvp_PCORR.get_Distance_Pair_DV_Series();
						projections_dist_PCORR = SS_dvp_PCORR.getProjections();
						normed_projections_dist_PCORR = SS_dvp_PCORR.getNormed_projections();
					}
			}

		private void get_Distance_Pair_PCA()
			{
				JEDi_Get_Distance_Pair_PCA dp_pca = new JEDi_Get_Distance_Pair_PCA(distance_matrix_cond, directory, description, number_of_modes_dist_pairs, number_of_pairs);
					{
						dp_pca.get_Distance_Pair_PCA();

						residue_distance_means = dp_pca.getResidue_means().getColumnPackedCopy();
						residue_distance_std_devs = dp_pca.getResidue_std_devs().getColumnPackedCopy();

						// COVARIANCE METHOD
						trace_dist_COV = dp_pca.getTrace_COV();
						cond_cov = dp_pca.getCond_COV();
						det_cov = dp_pca.getDet_COV();
						rank_cov = dp_pca.getRank_COV();
						top_distance_evectors_COV = dp_pca.getTop_evectors_dist_COV();
						top_distance_eigenvalues_COV = dp_pca.getTop_eigenvalues_COV();

						// CORRELATION METHOD
						trace_dist_CORR = dp_pca.getTrace_CORR();
						top_distance_evectors_CORR = dp_pca.getTop_evectors_dist_CORR();
						top_distance_eigenvalues_CORR = dp_pca.getTop_eigenvalues_CORR();

						// PARTIAL CORRELATION METHOD
						trace_dist_PCORR = dp_pca.getTrace_PCORR();
						top_distance_evectors_PCORR = dp_pca.getTop_evectors_dist_PCORR();
						top_distance_eigenvalues_PCORR = dp_pca.getTop_eigenvalues_PCORR();
					}
			}

		private void remove_Outliers_Z_Score()
			{
				if (z_cutoff > 0)
					{
						Adjust_Outliers_by_Z_Score cv = new Adjust_Outliers_by_Z_Score(distance_matrix);
							{
								cv.set_Z_threshold(z_cutoff);
								cv.adjust_row_data();
								distance_matrix_cond = cv.get_coorinates_adjusted();
								Z_scores = cv.get_z_scores();
								Matrix_IO.write_Matrix(Z_scores, out_dir + "ss_" + number_of_pairs + "_Atom_Pairs_Distance_Z_scores.txt", 6, 1);
								counts = cv.get_counts();
								Matrix_IO.write_Matrix(counts, out_dir + "ss_" + number_of_pairs + "_Atom_Pairs_Outliers_Per_Variable.txt", 6, 0);
							}
					}
				if (z_cutoff == 0) distance_matrix_cond = distance_matrix;
			}

		private void write_Distance_Pair_Stats()
			{
				try
					{
						path = out_dir + "ss_" + number_of_pairs + "_Atom_Pair_Distance_Stats.txt";
						File d_stats = new File(path);
						BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
						d_stats_writer.write("MEANs and STANDARD DEVIATIONS for the Atom Pair Distances: " + "\n");
						d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Atom1", "Atom2", "Mean", "Std_Dev"));
						for (int i = 0; i < number_of_pairs; i++)
							{
									{
										d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", atom_list_dist_orig1.get(i), atom_list_dist_orig2.get(i),
												nf.format(residue_distance_means[i]), nf.format(residue_distance_std_devs[i])));
									}
							}
						d_stats_writer.close();

					}
				catch (IOException io)
					{
						System.err.println("IOException thrown. Could not write the file: " + path);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}
			}

		private void get_Reference_Atom_Pair_Distances()
			{
				JEDi_Get_Distances_for_Atom_Pairs ref_dist_pairs = new JEDi_Get_Distances_for_Atom_Pairs(directory, description, reference_pdb_coordinates, atom_list_dist1,
						atom_list_dist2);
					{
						reference_distances = ref_dist_pairs.Get_Ref_Distances();
					}
			}

		private void get_Atom_Pair_Distances()
			{
				JEDi_Get_Distances_for_Atom_Pairs dist_pairs = new JEDi_Get_Distances_for_Atom_Pairs(directory, description, original_PDB_coordinates, atom_list_dist1,
						atom_list_dist2);
					{
						distance_matrix = dist_pairs.get_Distances();
					}
			}

		private void write_Distance_Pair_Stats_Multi()
			{
				try
					{
						path = out_dir + "ss_" + number_of_pairs + "_Atom_Pair_Distance_Stats.txt";
						File d_stats = new File(path);
						BufferedWriter d_stats_writer = new BufferedWriter(new FileWriter(d_stats));
						d_stats_writer.write("MEANs and STANDARD DEVIATIONS for the Atom Pair Distances: " + "\n");
						d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Atom1", "Atom2", "Mean", "Std_Dev"));
						for (int i = 0; i < number_of_pairs; i++)
							{
									{
										d_stats_writer.write(String.format("%-12s%-16s%-16s%-16s%n", chain_idents1.get(i) + atom_list_dist_orig1.get(i),
												chain_idents2.get(i) + atom_list_dist_orig2.get(i), nf.format(residue_distance_means[i]), nf.format(residue_distance_std_devs[i])));
									}
							}
						d_stats_writer.close();

					}
				catch (IOException io)
					{
						System.out.println("IOException thrown. Could not write the file: " + path);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}
			}

		private void do_FES()
			{
				System.out.println("Calculating the Free Energy Surface for the 3 PCA models... ");

				JEDi_Get_FES fes = new JEDi_Get_FES(weighted_projections_dist_CORR, 0, 1, number_of_conformations, 0, 0, directory, description, type, "CORR");
				fes.get_FES();
				fes.write_FES_Log();

				fes = new JEDi_Get_FES(weighted_projections_dist_COV, 0, 1, number_of_conformations, 0, 0, directory, description, type, "COV");
				fes.get_FES();
				fes.write_FES_Log();

				fes = new JEDi_Get_FES(projections_dist_PCORR, 0, 1, number_of_conformations, 0, 0, directory, description, type, "PCORR");
				fes.get_FES();
				fes.write_FES_Log();

				System.out.println("Done.");
			}


		/* ************************************* SETTERS ************************************************************************* ** */

		/**
		 * Sets the Z cutoff for outlier processing.
		 *
		 * @param z The Z cutoff
		 */
		public void set_z_cutoff(double z)
			{

				z_cutoff = z;
			}

		/* ************************************* GETTERS ************************************************************************* ** */

		/**
		 * @return the number_of_Chains
		 */
		public int getNumber_of_Chains()
			{
				return number_of_Chains;
			}

		/**
		 * @return the number_of_pairs
		 */
		public int getNumber_of_pairs()
			{
				return number_of_pairs;
			}

		/**
		 * @return the number_of_modes_dist_pairs
		 */
		public int getNumber_of_modes_dist_pairs()
			{
				return number_of_modes_dist_pairs;
			}

		/**
		 * @return the directory
		 */
		public String getDirectory()
			{
				return directory;
			}

		/**
		 * @return the out_dir
		 */
		public String getOut_dir()
			{
				return out_dir;
			}

		/**
		 * @return the description
		 */
		public String getDescription()
			{
				return description;
			}

		/**
		 * @return the path
		 */
		public String getPath()
			{
				return path;
			}

		/**
		 * @return the line
		 */
		public String getLine()
			{
				return line;
			}

		/**
		 * @return the trace_dist_COV
		 */
		public double getTrace_dist_COV()
			{
				return trace_dist_COV;
			}

		/**
		 * @return the trace_dist_CORR
		 */
		public double getTrace_dist_CORR()
			{
				return trace_dist_CORR;
			}

		/**
		 * @return the trace_dist_PCORR
		 */
		public double getTrace_dist_PCORR()
			{
				return trace_dist_PCORR;
			}

		/**
		 * @return the cond_cov
		 */
		public double getCond_cov()
			{
				return cond_cov;
			}

		/**
		 * @return the cond_corr
		 */
		public double getCond_corr()
			{
				return cond_corr;
			}

		/**
		 * @return the cond_pcorr
		 */
		public double getCond_pcorr()
			{
				return cond_pcorr;
			}

		/**
		 * @return the det_cov
		 */
		public double getDet_cov()
			{
				return det_cov;
			}

		/**
		 * @return the det_corr
		 */
		public double getDet_corr()
			{
				return det_corr;
			}

		/**
		 * @return the det_pcorr
		 */
		public double getDet_pcorr()
			{
				return det_pcorr;
			}

		/**
		 * @return the rank_cov
		 */
		public double getRank_cov()
			{
				return rank_cov;
			}

		/**
		 * @return the rank_corr
		 */
		public double getRank_corr()
			{
				return rank_corr;
			}

		/**
		 * @return the rank_pcorr
		 */
		public double getRank_pcorr()
			{
				return rank_pcorr;
			}

		/**
		 * @return the z_cutoff
		 */
		public double getZ_cutoff()
			{
				return z_cutoff;
			}

		/**
		 * @return the residue_list_dist1
		 */
		public List<Integer> getResidue_list_dist1()
			{
				return atom_list_dist1;
			}

		/**
		 * @return the residue_list_dist_orig1
		 */
		public List<Integer> getResidue_list_dist_orig1()
			{
				return atom_list_dist_orig1;
			}

		/**
		 * @return the residue_list_dist2
		 */
		public List<Integer> getResidue_list_dist2()
			{
				return atom_list_dist2;
			}

		/**
		 * @return the residue_list_dist_orig2
		 */
		public List<Integer> getResidue_list_dist_orig2()
			{
				return atom_list_dist_orig2;
			}

		/**
		 * @return the chain_idents1
		 */
		public List<String> getChain_idents1()
			{
				return chain_idents1;
			}

		/**
		 * @return the chain_idents2
		 */
		public List<String> getChain_idents2()
			{
				return chain_idents2;
			}

		/**
		 * @return the top_distance_eigenvalues_COV
		 */
		public List<Double> getTop_distance_eigenvalues_COV()
			{
				return top_distance_eigenvalues_COV;
			}

		/**
		 * @return the top_distance_eigenvalues_CORR
		 */
		public List<Double> getTop_distance_eigenvalues_CORR()
			{
				return top_distance_eigenvalues_CORR;
			}

		/**
		 * @return the top_distance_eigenvalues_PCORR
		 */
		public List<Double> getTop_distance_eigenvalues_PCORR()
			{
				return top_distance_eigenvalues_PCORR;
			}

		/**
		 * @return the residue_distance_means
		 */
		public double[] getResidue_distance_means()
			{
				return residue_distance_means;
			}

		/**
		 * @return the residue_distance_std_devs
		 */
		public double[] getResidue_distance_std_devs()
			{
				return residue_distance_std_devs;
			}

		/**
		 * @return the dist_pca_mode_COV_min
		 */
		public double[] getDist_pca_mode_COV_min()
			{
				return dist_pca_mode_COV_min;
			}

		/**
		 * @return the dist_pca_mode_COV_max
		 */
		public double[] getDist_pca_mode_COV_max()
			{
				return dist_pca_mode_COV_max;
			}

		/**
		 * @return the dist_pca_mode_CORR_min
		 */
		public double[] getDist_pca_mode_CORR_min()
			{
				return dist_pca_mode_CORR_min;
			}

		/**
		 * @return the dist_pca_mode_CORR_max
		 */
		public double[] getDist_pca_mode_CORR_max()
			{
				return dist_pca_mode_CORR_max;
			}

		/**
		 * @return the distance_matrix
		 */
		public Matrix getDistance_matrix()
			{
				return distance_matrix;
			}

		/**
		 * @return the distance_matrix_cond
		 */
		public Matrix getDistance_matrix_cond()
			{
				return distance_matrix_cond;
			}

		/**
		 * @return the original_PDB_coordinates
		 */
		public Matrix getOriginal_PDB_coordinates()
			{
				return original_PDB_coordinates;
			}

		/**
		 * @return the subset_PDB_coordinates_dist
		 */
		public Matrix getSubset_PDB_coordinates_dist()
			{
				return subset_PDB_coordinates_dist;
			}

		/**
		 * @return the transformed_subset_PDB_coordinates_dist
		 */
		public Matrix getTransformed_subset_PDB_coordinates_dist()
			{
				return transformed_subset_PDB_coordinates_dist;
			}

		/**
		 * @return the cov_dist
		 */
		public Matrix getCov_dist()
			{
				return cov_dist;
			}

		/**
		 * @return the corr_dist
		 */
		public Matrix getCorr_dist()
			{
				return corr_dist;
			}

		/**
		 * @return the top_distance_evectors_COV
		 */
		public Matrix getTop_distance_evectors_COV()
			{
				return top_distance_evectors_COV;
			}

		/**
		 * @return the top_distance_evectors_CORR
		 */
		public Matrix getTop_distance_evectors_CORR()
			{
				return top_distance_evectors_CORR;
			}

		/**
		 * @return the top_distance_evectors_PCORR
		 */
		public Matrix getTop_distance_evectors_PCORR()
			{
				return top_distance_evectors_PCORR;
			}

		/**
		 * @return the normed_projections_dist_COV
		 */
		public Matrix getNormed_projections_dist_COV()
			{
				return normed_projections_dist_COV;
			}

		/**
		 * @return the projections_dist_COV
		 */
		public Matrix getProjections_dist_COV()
			{
				return projections_dist_COV;
			}

		/**
		 * @return the projections_dist_CORR
		 */
		public Matrix getProjections_dist_CORR()
			{
				return projections_dist_CORR;
			}

		/**
		 * @return the normed_projections_dist_CORR
		 */
		public Matrix getNormed_projections_dist_CORR()
			{
				return normed_projections_dist_CORR;
			}

		/**
		 * @return the weighted_normed_projections_dist_COV
		 */
		public Matrix getWeighted_normed_projections_dist_COV()
			{
				return weighted_normed_projections_dist_COV;
			}

		/**
		 * @return the weighted_normed_projections_dist_CORR
		 */
		public Matrix getWeighted_normed_projections_dist_CORR()
			{
				return weighted_normed_projections_dist_CORR;
			}

		/**
		 * @return the weighted_projections_dist_COV
		 */
		public Matrix getWeighted_projections_dist_COV()
			{
				return weighted_projections_dist_COV;
			}

		/**
		 * @return the weighted_projections_dist_CORR
		 */
		public Matrix getWeighted_projections_dist_CORR()
			{
				return weighted_projections_dist_CORR;
			}

		/**
		 * @return the normed_projections_dist_PCORR
		 */
		public Matrix getNormed_projections_dist_PCORR()
			{
				return normed_projections_dist_PCORR;
			}

		/**
		 * @return the projections_dist_PCORR
		 */
		public Matrix getProjections_dist_PCORR()
			{
				return projections_dist_PCORR;
			}

		/**
		 * @return the z_scores
		 */
		public Matrix getZ_scores()
			{
				return Z_scores;
			}

		/**
		 * @return the counts
		 */
		public Matrix getCounts()
			{
				return counts;
			}

		/**
		 * @return the number format
		 */
		public NumberFormat getNf()
			{
				return nf;
			}

		/**
		 * @return the rounding
		 */
		public RoundingMode getRm()
			{
				return rm;
			}
	}
