package jedi;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Do_Generalized_Cartesian: Top class for implementing the Generalized Cartesian analysis. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Do_Hierarchical_PCA
	{

		boolean exist, success;
		int number_of_residues, number_modes_RGC, number_of_atoms, number_of_modes_SS, number_of_modes_Residues, number_of_conformations;
		double trace, cond, det, rank;
		String directory, out_dir, description;
		double[] pca_mode_max, pca_mode_min;
		Matrix big_G, G_Square_Modes, big_DV, DVPs, gc_coordinates, cov, centered_data;
		Matrix top_cartesian_evectors, square_pca_modes, weighted_square_pca_modes, weighted_pca_modes, pca_modes;
		Matrix normed_projections_rgc, projections_rgc, weighted_projections_rgc, weighted_normed_projections_rgc, normed_projections_pc, projections_pc, weighted_projections_pc,
				weighted_normed_projections_pc;
		List<Double> top_eigenvalues, eigenvalues;
		List<Matrix> residue_Gs, residue_Eigenvectors, residue_DVs;
		JEDi_Get_FES fes;

		/* ****************************************************** CONSTRUCTOR ************************************************************************************** */

		public JEDi_Do_Hierarchical_PCA(String directory, String desc, int number_of_modes_SS, Matrix gc_cov, List<Matrix> Residue_Evects, List<Matrix> Residue_DVs, int modes_res)
			{
				super();
				this.directory = directory;
				this.description = desc;
				this.number_of_modes_SS = number_of_modes_SS;
				this.gc_coordinates = gc_cov;
				this.residue_Eigenvectors = Residue_Evects;
				this.residue_DVs = Residue_DVs;
				this.number_of_modes_Residues = modes_res;
				this.number_of_residues = (gc_cov.getRowDimension() / number_of_modes_Residues);

				residue_Gs = new ArrayList<Matrix>();

				out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Hierarchical_PCA" + File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();
			}

		/* *********************************************************** PUBLIC DRIVER METHOD *************************************************** */

		public void do_Hierarchical_PCA()
			{
				get_Hierarchical_PCA();
				get_Hierarchical_PCs();
				get_Big_G();
				get_big_DVs();
				get_big_DVPs();
				do_FES();
			}

		/* ************************************************************** PRIVATE METHODS *************************************************** */

		private void get_Hierarchical_PCA()
			{
				JEDi_Get_Hierarchical_PCA gcov_pca = new JEDi_Get_Hierarchical_PCA(gc_coordinates, number_of_modes_SS, directory, description, number_of_modes_Residues);
					{
						gcov_pca.get_Hierarchical_PCA();

						// COVARIANCE METHOD ONLY
						cond = gcov_pca.get_cond_COV();
						trace = gcov_pca.get_trace_COV();
						det = gcov_pca.get_det_COV();
						rank = gcov_pca.get_rank_COV();
						eigenvalues = gcov_pca.getEigenvalues_COV();
						top_cartesian_evectors = gcov_pca.getTop_evectors_COV();
						top_eigenvalues = gcov_pca.getTop_eigenvalues_COV();
						pca_modes = gcov_pca.getPca_modes_COV();
						square_pca_modes = gcov_pca.getSquare_pca_modes_COV();
						weighted_square_pca_modes = gcov_pca.getWeighted_square_pca_modes_COV();
						weighted_pca_modes = gcov_pca.getWeighted_pca_modes_COV();
						pca_mode_min = gcov_pca.get_pca_mode_COV_min();
						pca_mode_max = gcov_pca.get_pca_mode_COV_max();
						centered_data = gcov_pca.getCentered_input_data();
					}
			}

		private void get_Hierarchical_PCs()
			{
				JEDi_Get_Hierarchical_PCs pcs_cov = new JEDi_Get_Hierarchical_PCs(centered_data, top_cartesian_evectors, eigenvalues, directory,
						description, number_of_modes_Residues);
					{
						pcs_cov.get_PCs();
						projections_pc = pcs_cov.getPCs();
						normed_projections_pc = pcs_cov.getNormed_PCs();
						weighted_projections_pc = pcs_cov.getWeighted_PCs();
						weighted_normed_projections_pc = pcs_cov.getWeighted_normed_PCs();
					}
			}

		private void get_Big_G()
			{
				JEDi_Get_Big_G getBG = new JEDi_Get_Big_G(directory, description, number_of_modes_SS, number_of_modes_Residues, top_cartesian_evectors, residue_Eigenvectors);
				big_G = getBG.get_big_G_Col();
				number_of_atoms = getBG.getNumber_of_atoms();

				int rowsG = big_G.getRowDimension();
				int colsG = big_G.getColumnDimension();

				G_Square_Modes = new Matrix(rowsG, colsG);

				for (int i = 0; i < rowsG; i++)
					{
						for (int j = 0; j < colsG; j++)
							{
								double val = big_G.get(i, j);
								double val_SQ = val * val;
								G_Square_Modes.set(i, j, val_SQ);
							}
					}

				String path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_big_G.txt";
				Matrix_IO.write_Matrix(big_G, path, 12, 9);

				path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_big_G_Square_Modes.txt";
				Matrix_IO.write_Matrix(G_Square_Modes, path, 12, 9);
			}

		private void get_big_DVs()
			{
				number_of_conformations = residue_DVs.get(0).getColumnDimension();
				big_DV = new Matrix(number_of_atoms * 3, number_of_conformations);

				int big_DV_offset = 0;
				for (int i = 0; i < number_of_residues; i++)
					{
						Matrix dv = residue_DVs.get(i);
						int num_of_atoms_Res = dv.getRowDimension() / 3;
						for (int j = 0; j < number_of_conformations; j++)
							{
								for (int k = 0; k < num_of_atoms_Res; k++)
									{
										double X = dv.get(k, j);
										double Y = dv.get(k + num_of_atoms_Res, j);
										double Z = dv.get(k + 2 * num_of_atoms_Res, j);
										big_DV.set(k + big_DV_offset, j, X);
										big_DV.set(k + big_DV_offset + number_of_atoms, j, Y);
										big_DV.set(k + big_DV_offset + 2 * number_of_atoms, j, Z);
									}
							}
					}
				String path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_big_DVs.txt";
				Matrix_IO.write_Matrix(big_DV, path, 12, 9);
			}

		private void get_big_DVPs()
			{
				projections_rgc = new Matrix(number_of_conformations, number_of_modes_Residues);
				weighted_projections_rgc = new Matrix(number_of_conformations, number_of_modes_Residues);
				normed_projections_rgc = new Matrix(number_of_conformations, number_of_modes_Residues);
				weighted_normed_projections_rgc = new Matrix(number_of_conformations, number_of_modes_Residues);

				for (int outer = 0; outer < number_of_modes_Residues; outer++)
					{
						for (int inner = 0; inner < number_of_conformations; inner++)
							{
								int row_index_1 = 0;
								int row_index_2 = 3 * number_of_atoms - 1;

								Matrix data1 = big_G.getMatrix(row_index_1, row_index_2, outer, outer);
								Matrix data2 = big_DV.getMatrix(row_index_1, row_index_2, inner, inner);
								Matrix vector1 = Projector.get_Normed_array(data1);
								Matrix vector2 = Projector.get_Normed_array(data2);
								double weight = eigenvalues.get(outer);
								if (weight < 0)
									weight = 0; // Fixes the -0.00000000 problem.
								double sqrt_weight = Math.sqrt(weight);
								double dp = Projector.get_InnerProduct(data1, data2);
								double normed_dp = Projector.get_InnerProduct(vector1, vector2);
								if (dp == Double.NaN)
									dp = 0.000;
								if (normed_dp == Double.NaN)
									normed_dp = 0.000;
								double w_dp = sqrt_weight * dp;
								double weighted_normed_dp = sqrt_weight * normed_dp;
								projections_rgc.set(inner, outer, dp);
								weighted_projections_rgc.set(inner, outer, w_dp);
								normed_projections_rgc.set(inner, outer, normed_dp);
								weighted_normed_projections_rgc.set(inner, outer, weighted_normed_dp);
							}
					}
				String path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_DVPs.txt";
				Matrix_IO.write_Matrix(projections_rgc, path, 12, 9);

				path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_normed_DVPs.txt";
				Matrix_IO.write_Matrix(normed_projections_rgc, path, 12, 9);

				path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_weighted_DVPs.txt";
				Matrix_IO.write_Matrix(weighted_projections_rgc, path, 9, 6);

				path = out_dir + "ss_" + number_of_residues + "_top_" + number_of_modes_Residues + "_weighted_normed_DVPs.txt";
				Matrix_IO.write_Matrix(weighted_normed_projections_rgc, path, 9, 6);
			}

		private void do_FES()
			{
				System.out.println("Calculating the Free Energy Surface... ");

				fes = new JEDi_Get_FES(weighted_projections_rgc, 0, 1, number_of_conformations, 0, 0, directory, description, "Hierarchical_PCA", "COV");
				fes.get_FES();
				fes.write_FES_Log();

				System.out.println("Done.");
			}


		/* ************************************************************** GETTERS ******************************************************************* */

		public double get_Trace()
			{

				return trace;
			}


		public double get_cond()
			{

				return cond;
			}

		public double get_Cond()
			{
				return cond;
			}

		public double get_Det()
			{
				return det;
			}

		public double get_Rank()
			{
				return rank;
			}

		public List<Double> getTop_cartesian_eigenvalues()
			{

				return top_eigenvalues;
			}

		public double[] getPca_mode_max()
			{

				return pca_mode_max;
			}

		public double[] getPca_mode_min()
			{

				return pca_mode_min;
			}

		public Matrix getCov()
			{

				return cov;
			}

		public Matrix getTop_cartesian_evectors()
			{

				return top_cartesian_evectors;
			}

		public Matrix getWeighted_pca_modes()
			{

				return weighted_pca_modes;
			}

		public Matrix getSquare_pca_modes()
			{

				return square_pca_modes;
			}

		public Matrix getWeighted_square_pca_modes()
			{

				return weighted_square_pca_modes;
			}

		public Matrix getNormed_projections()
			{

				return normed_projections_rgc;
			}

		public Matrix getProjections()
			{

				return projections_rgc;
			}

		public Matrix getWeighted_projections()
			{
				return weighted_projections_rgc;
			}

		public Matrix getWeighted_normed_projections()
			{
				return weighted_normed_projections_rgc;
			}

		public int getNumber_of_residues()
			{

				return number_of_residues;
			}

		public Matrix getPca_modes()
			{

				return pca_modes;
			}

		public double getTrace()
			{
				return trace;
			}

		public double getDet()
			{
				return det;
			}

		public double getRank()
			{
				return rank;
			}

		public List<Double> getEigenvalues()
			{
				return eigenvalues;
			}

		public Matrix getBig_G()
			{
				return big_G;
			}

		public Matrix getBig_DV()
			{
				return big_DV;
			}

		public Matrix getDVPs()
			{
				return DVPs;
			}

		public Matrix getNormed_projections_rgc()
			{
				return normed_projections_rgc;
			}

		public Matrix getProjections_rgc()
			{
				return projections_rgc;
			}

		public Matrix getWeighted_projections_rgc()
			{
				return weighted_projections_rgc;
			}

		public Matrix getWeighted_normed_projections_rgc()
			{
				return weighted_normed_projections_rgc;
			}

		public Matrix getNormed_projections_pc()
			{
				return normed_projections_pc;
			}

		public Matrix getProjections_pc()
			{
				return projections_pc;
			}

		public Matrix getWeighted_projections_pc()
			{
				return weighted_projections_pc;
			}

		public Matrix getWeighted_normed_projections_pc()
			{
				return weighted_normed_projections_pc;
			}

		public List<Matrix> getResidue_Gs()
			{
				return residue_Gs;
			}

		public Matrix getG_Square_Modes()
			{
				return G_Square_Modes;
			}
	}
