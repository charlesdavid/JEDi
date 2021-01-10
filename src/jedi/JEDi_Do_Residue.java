package jedi;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import support.Atom;
import support.Residue_ID_Pair;
import support.Row_Center_Data;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;
import supportPlot.HeatMap;
import supportPlot.MODES_Plot;
import supportPlot.RMSD_Plot;
import supportPlot.STATS_Plot;

/**
 * JED class JED_Do_Residue: Top class for implementing the Residue Cartesian analysis
 * 
 * Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Do_Residue
{
	boolean exist, success, verbose;
	boolean doCORR, doPCORR, doSPARSIFY, doReduce, doFES, doKPCA, do_ESSENTIAL_VIZ, do_INDIVIDUAL_VIZ, doLOCAL, doGLOBAL, doHeavyAtoms, doAllAtoms;

	String local_out_dir, global_out_dir, pairs_out_dir, out_dir_pca, out_dir_ssa, out_dir_fes, out_dir_kpca, out_dir_viz;

	String res_index, res_index1, res_index2, Residue_Pair_Index, type_alignment;

	final String cRes_AA = "Residue_All_Atom_PCA", Q = "COV", R = "CORR", P = "PCORR";

	int offset, number_Of_Atoms_In_Residue, number_Of_Atoms_In_Residue1, number_Of_Atoms_In_Residue2, residue_Index, residue_Index_Orig, residue_Index1_Orig, residue_Index2_Orig;
	final int number_of_conformations, number_modes_Residue, number_of_modes_viz, number_residues;

	double trace_COV, trace_CORR, trace_CORR_SPARSE, trace_PCORR, trace_PCORR_SPARSE, cond_COV, det_COV, rank_COV;
	final double MAX_COUPLING_SCORE = 100;

	List<Double> transformed_conformation_rmsds, transformed_residue_rmsd_list;
	List<Double> top_eigenvalues_COV, top_eigenvalues_CORR, top_eigenvalues_PCORR, eigenvalues_COV, eigenvalues_CORR, eigenvalues_PCORR;
	List<Double> top_eigenvalues_CORR_SPARSE, top_eigenvalues_PCORR_SPARSE, eigenvalues_CORR_SPARSE, eigenvalues_PCORR_SPARSE;

	double[] pca_mode_maxes_COV, pca_mode_mins_COV, pca_mode_maxes_CORR, pca_mode_mins_CORR, pca_mode_maxes_PCORR, pca_mode_mins_PCORR;
	double[] pca_mode_maxes_CORR_SPARSE, pca_mode_mins_CORR_SPARSE, pca_mode_maxes_PCORR_SPARSE, pca_mode_mins_PCORR_SPARSE;

	Matrix cov, coupling_scores;
	Matrix residue_PDB_coordinates, adjusted_residue_PDB_coordinates, residue_PDB_coordinates1, residue_PDB_coordinates2, residue_REF_PDB_coordinates, residue_REF_PDB_coordinates1,
			residue_REF_PDB_coordinates2;
	Matrix original_Residue_REF_coords, aligned_Residue_REF_coords, aligned_Residue_coords, aligned_Residue_coords_Outliers_REMOVED, aligned_Residue_coords_Outliers_SELECTED,
			delta_vectors, mean_centered_data;
	Matrix evectors_COV, evectors_CORR, evectors_PCORR, evectors_CORR_SPARSE, evectors_PCORR_SPARSE;

	Matrix pca_modes_COV, square_pca_modes_COV, weighted_square_pca_modes_COV, weighted_pca_modes_COV;
	Matrix pca_modes_CORR, square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR;
	Matrix pca_modes_CORR_SPARSE, square_pca_modes_CORR_SPARSE, weighted_square_pca_modes_CORR_SPARSE, weighted_pca_modes_CORR_SPARSE;
	Matrix pca_modes_PCORR, square_pca_modes_PCORR, weighted_square_pca_modes_PCORR, weighted_pca_modes_PCORR;
	Matrix pca_modes_PCORR_SPARSE, square_pca_modes_PCORR_SPARSE, weighted_square_pca_modes_PCORR_SPARSE, weighted_pca_modes_PCORR_SPARSE;

	Matrix projections_COV, normed_projections_COV, weighted_projections_COV, weighted_normed_projections_COV;
	Matrix projections_CORR, normed_projections_CORR, weighted_projections_CORR, weighted_normed_projections_CORR;
	Matrix projections_CORR_SPARSE, normed_projections_CORR_SPARSE, weighted_projections_CORR_SPARSE, weighted_normed_projections_CORR_SPARSE;
	Matrix projections_PCORR, normed_projections_PCORR, weighted_projections_PCORR, weighted_normed_projections_PCORR;
	Matrix projections_PCORR_SPARSE, normed_projections_PCORR_SPARSE, weighted_projections_PCORR_SPARSE, weighted_normed_projections_PCORR_SPARSE;

	Matrix Residue_Generalized_Coordinates_Pair, Residue_Generalized_Coordinates_Global;

	final Matrix subset_PDB_coordinates, subset_reference_PDB_coordinates, adjusted_subset_PDB_coordinates;

	List<Matrix> Residue_Eigenvectors_Global, Residues_Centered_Data_Global;

	final List<List<Atom>> Residue_Atoms_List;
	List<Atom> residueAtoms, residueAtoms1, residueAtoms2;
	List<Residue_ID_Pair> residue_list;

	/* ********************************************************* CONSTRUCTORS **************************************************************** */
	/**
	 * Constructor for INDIVIDUAL Residue Analysis
	 */
	public JEDi_Do_Residue(List<Residue_ID_Pair> id_pairs, Matrix PDB_coordinates, Matrix Ref_PDB_coords, Matrix Adj_PDB_coords, List<List<Atom>> residueAtomsList)
	{
		this.residue_list = id_pairs;
		this.subset_PDB_coordinates = PDB_coordinates;
		this.subset_reference_PDB_coordinates = Ref_PDB_coords;
		this.adjusted_subset_PDB_coordinates = Adj_PDB_coords;
		this.Residue_Atoms_List = residueAtomsList;

		this.number_residues = residue_list.size();
		this.number_of_conformations = subset_PDB_coordinates.getColumnDimension();

		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doSPARSIFY = Input_Parameters.doSPARSIFY;
		this.doReduce = Input_Parameters.doREDUCE;
		this.doFES = Input_Parameters.doFES;
		this.doKPCA = Input_Parameters.doKPCA;
		this.do_ESSENTIAL_VIZ = Input_Parameters.doEssentialViz;
		this.do_INDIVIDUAL_VIZ = Input_Parameters.doModeViz;
		this.number_modes_Residue = Input_Parameters.MODES_RESIDUE_INDIVIDUAL;
		this.number_of_modes_viz = Input_Parameters.MODES_VIZ;
		this.verbose = Input_Parameters.verbose;
	}

	/**
	 * Constructor for Residue PAIRS Analysis
	 */
	public JEDi_Do_Residue(List<Residue_ID_Pair> id_pairs, Matrix PDB_coordinates, Matrix Ref_PDB_coords, List<List<Atom>> residueAtomsList)
	{
		this.residue_list = id_pairs;
		this.subset_PDB_coordinates = PDB_coordinates;
		this.subset_reference_PDB_coordinates = Ref_PDB_coords;
		this.Residue_Atoms_List = residueAtomsList;
		this.adjusted_subset_PDB_coordinates = null;

		this.number_residues = residue_list.size();
		this.number_of_conformations = subset_PDB_coordinates.getColumnDimension();

		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doSPARSIFY = Input_Parameters.doSPARSIFY;
		this.doReduce = Input_Parameters.doREDUCE;
		this.doFES = Input_Parameters.doFES;
		this.doKPCA = Input_Parameters.doKPCA;
		this.do_ESSENTIAL_VIZ = Input_Parameters.doEssentialViz;
		this.do_INDIVIDUAL_VIZ = Input_Parameters.doModeViz;
		this.number_modes_Residue = Input_Parameters.MODES_EIGEN_RESIDUE_PAIRS;
		this.number_of_modes_viz = Input_Parameters.MODES_VIZ;
		this.verbose = Input_Parameters.verbose;
	}

	/**
	 * Constructor for GLOBAL Residue Analysis --> Basis for Hierarchical PCA
	 */
	public JEDi_Do_Residue(List<Residue_ID_Pair> id_pairs, int number_of_modes, Matrix PDB_coordinates, Matrix Ref_PDB_coords, Matrix Adj_PDB_coords,
			List<List<Atom>> residueAtomsList, int modes_viz)
	{
		this.residue_list = id_pairs;
		this.number_residues = residue_list.size();
		this.number_modes_Residue = number_of_modes;
		this.subset_PDB_coordinates = PDB_coordinates;
		this.subset_reference_PDB_coordinates = Ref_PDB_coords;
		this.adjusted_subset_PDB_coordinates = Adj_PDB_coords;

		this.number_of_conformations = subset_PDB_coordinates.getColumnDimension();
		this.Residue_Atoms_List = residueAtomsList;
		this.number_of_modes_viz = modes_viz;

		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doReduce = Input_Parameters.doREDUCE;
		this.doFES = Input_Parameters.doFES;
		this.doKPCA = Input_Parameters.doKPCA;
		this.do_ESSENTIAL_VIZ = Input_Parameters.doEssentialViz;
		this.do_INDIVIDUAL_VIZ = Input_Parameters.doModeViz;
		this.verbose = Input_Parameters.verbose;

		Residue_Generalized_Coordinates_Global = new Matrix(residue_list.size() * number_modes_Residue, number_of_conformations);
		Residue_Eigenvectors_Global = new ArrayList<Matrix>();
		Residues_Centered_Data_Global = new ArrayList<Matrix>();
	}

	/* ******************************************************** PUBLIC METHODS *************************************************************** */

	public void do_Individual_Residue_Analysis() // Uses local alignment of each residue individually
	{
		type_alignment = "Local";
		offset = 0;
		String out;

		for (int i = 0; i < number_residues; i++)
			{
				String chainID = residue_list.get(i).getChain_ID();
				residue_Index_Orig = residue_list.get(i).getResidue_Number();
				res_index = String.format("%s_%03d", chainID, residue_Index_Orig);
				residueAtoms = Residue_Atoms_List.get(i);
				number_Of_Atoms_In_Residue = residueAtoms.size();

				local_out_dir = out_dir_pca + ("Residue_" + res_index) + File.separatorChar;
				create_Directory(local_out_dir);

				JEDi_Get_Residue_Coords ss = new JEDi_Get_Residue_Coords(subset_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
					{
						residue_PDB_coordinates = ss.get_residue_Coords();
					}

				JEDi_Get_Residue_Coords ref_ss = new JEDi_Get_Residue_Coords(subset_reference_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
					{
						residue_REF_PDB_coordinates = ref_ss.get_residue_Coords();
					}

				get_Local_Transformed_Coords();
				get_Local_Residue_PCA();
				get_Local_Residue_DVPs();


				out = out_dir_ssa + ("Residue_" + res_index) + File.separatorChar;
				create_Directory(out);

				do_SSA(out);

				if (doFES)
					{
						out = out_dir_fes + "Residue_" + res_index + File.separatorChar;
						create_Directory(out);
						do_FES(out);
					}

				if (doKPCA)
					{
						out = out_dir_kpca + "Residue_" + res_index + File.separatorChar;
						create_Directory(out);
						do_KPCA(out);
					}

				if (do_ESSENTIAL_VIZ)
					{
						out = out_dir_viz + "Residue_" + res_index + File.separatorChar;
						create_Directory(out);
						get_Local_Residue_Mode_Visualization(out);
					}

				offset += number_Of_Atoms_In_Residue;
			}
	}

	public void do_Global_Residue_Analysis() // Uses global alignment of selected subset. Needed for consistency in the Hierarchical PCA Method
	{
		type_alignment = "Global";
		offset = 0;

		for (int i = 0; i < number_residues; i++)
			{
				residue_Index = (i);
				String chainID = residue_list.get(i).getChain_ID();
				residue_Index_Orig = residue_list.get(i).getResidue_Number();
				res_index = String.format("%s_%03d", chainID, residue_Index_Orig);
				residueAtoms = Residue_Atoms_List.get(residue_Index);
				number_Of_Atoms_In_Residue = residueAtoms.size();

				JEDi_Get_Residue_Coords ss = new JEDi_Get_Residue_Coords(subset_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
					{
						residue_PDB_coordinates = ss.get_residue_Coords();
					}

				JEDi_Get_Residue_Coords ref_ss = new JEDi_Get_Residue_Coords(subset_reference_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
					{
						residue_REF_PDB_coordinates = ref_ss.get_residue_Coords();
					}

				JEDi_Get_Residue_Coords adj_ss = new JEDi_Get_Residue_Coords(adjusted_subset_PDB_coordinates, number_Of_Atoms_In_Residue, offset);
					{
						adjusted_residue_PDB_coordinates = adj_ss.get_residue_Coords();
					}

				get_Global_Residue_PCA(residue_PDB_coordinates, adjusted_residue_PDB_coordinates);
				get_Global_Residue_PCs();
				offset += number_Of_Atoms_In_Residue;
			}
	}

	public void do_Residue_Pairs(List<Matrix> ref_residue_coords, List<Matrix> residue_coords) // Uses global alignment of each residue PAIR
	{
		type_alignment = "Pairs";
		coupling_scores = new Matrix(number_residues, number_residues);

		for (int i = 0; i < number_residues; i++) // Residue 1 index
			{
				String chainID1 = residue_list.get(i).getChain_ID();
				residue_Index1_Orig = residue_list.get(i).getResidue_Number();
				res_index1 = String.format("%s_%03d", chainID1, residue_Index1_Orig);
				residueAtoms1 = Residue_Atoms_List.get(i);
				number_Of_Atoms_In_Residue1 = residueAtoms1.size();

				Matrix RC1 = residue_coords.get(i);
				Matrix RC1_ref = ref_residue_coords.get(i);

				for (int j = i + 1; j < number_residues; j++) // Residue 2 index
					{
						String chainID2 = residue_list.get(j).getChain_ID();
						residue_Index2_Orig = residue_list.get(j).getResidue_Number();
						res_index2 = String.format("%s_%03d", chainID2, residue_Index2_Orig);
						residueAtoms2 = Residue_Atoms_List.get(j);
						number_Of_Atoms_In_Residue2 = residueAtoms2.size();

						Matrix RC2 = residue_coords.get(j);
						Matrix RC2_ref = ref_residue_coords.get(j);

						// Create the reference coordinates matrix for the pair:
						int number_of_atoms_in_pair = number_Of_Atoms_In_Residue1 + number_Of_Atoms_In_Residue2;

						Matrix pair_ref = new Matrix(3 * number_of_atoms_in_pair, 1);
						Matrix X1r = RC1_ref.getMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, 0);
						Matrix Y1r = RC1_ref.getMatrix(number_Of_Atoms_In_Residue1, 2 * number_Of_Atoms_In_Residue1 - 1, 0, 0);
						Matrix Z1r = RC1_ref.getMatrix(2 * number_Of_Atoms_In_Residue1, 3 * number_Of_Atoms_In_Residue1 - 1, 0, 0);
						Matrix X2r = RC2_ref.getMatrix(0, number_Of_Atoms_In_Residue2 - 1, 0, 0);
						Matrix Y2r = RC2_ref.getMatrix(number_Of_Atoms_In_Residue2, 2 * number_Of_Atoms_In_Residue2 - 1, 0, 0);
						Matrix Z2r = RC2_ref.getMatrix(2 * number_Of_Atoms_In_Residue2, 3 * number_Of_Atoms_In_Residue2 - 1, 0, 0);

						pair_ref.setMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, 0, X1r);
						pair_ref.setMatrix(number_of_atoms_in_pair, number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0, 0, Y1r);
						pair_ref.setMatrix(2 * number_of_atoms_in_pair, 2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0, 0, Z1r);
						pair_ref.setMatrix(number_Of_Atoms_In_Residue1, number_of_atoms_in_pair - 1, 0, 0, X2r);
						pair_ref.setMatrix(number_of_atoms_in_pair + number_Of_Atoms_In_Residue1, 2 * number_of_atoms_in_pair - 1, 0, 0, Y2r);
						pair_ref.setMatrix(2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue1, 3 * number_of_atoms_in_pair - 1, 0, 0, Z2r);

						// Create the coordinates matrix for the pair:
						Matrix pair = new Matrix(3 * number_of_atoms_in_pair, number_of_conformations);
						Matrix X1 = RC1.getMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1);
						Matrix Y1 = RC1.getMatrix(number_Of_Atoms_In_Residue1, 2 * number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1);
						Matrix Z1 = RC1.getMatrix(2 * number_Of_Atoms_In_Residue1, 3 * number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1);
						Matrix X2 = RC2.getMatrix(0, number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1);
						Matrix Y2 = RC2.getMatrix(number_Of_Atoms_In_Residue2, 2 * number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1);
						Matrix Z2 = RC2.getMatrix(2 * number_Of_Atoms_In_Residue2, 3 * number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1);

						pair.setMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, X1);
						pair.setMatrix(number_of_atoms_in_pair, number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, Y1);
						pair.setMatrix(2 * number_of_atoms_in_pair, 2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, Z1);
						pair.setMatrix(number_Of_Atoms_In_Residue1, number_of_atoms_in_pair - 1, 0, number_of_conformations - 1, X2);
						pair.setMatrix(number_of_atoms_in_pair + number_Of_Atoms_In_Residue1, 2 * number_of_atoms_in_pair - 1, 0, number_of_conformations - 1, Y2);
						pair.setMatrix(2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue1, 3 * number_of_atoms_in_pair - 1, 0, number_of_conformations - 1, Z2);

						// Align the pair as a subset:
						JEDi_Get_Transformed_Residue_Coordinates tf_coords = new JEDi_Get_Transformed_Residue_Coordinates(pair, pair_ref);
						aligned_Residue_REF_coords = tf_coords.get_Transformed_reference_coordinates();
						aligned_Residue_coords = tf_coords.get_SS_Transformed_coords();

						// Get the Pair-Aligned Coordinates for each residue in the pair:
						// ---> Residue1:
						Matrix residue_1_pair_aligned_coords = new Matrix(3 * number_Of_Atoms_In_Residue1, number_of_conformations);

						X1 = aligned_Residue_coords.getMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1);
						Y1 = aligned_Residue_coords.getMatrix(number_of_atoms_in_pair, number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1);
						Z1 = aligned_Residue_coords.getMatrix(2 * number_of_atoms_in_pair, 2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue1 - 1, 0,
								number_of_conformations - 1);

						residue_1_pair_aligned_coords.setMatrix(0, number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, X1);
						residue_1_pair_aligned_coords.setMatrix(number_Of_Atoms_In_Residue1, 2 * number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, X1);
						residue_1_pair_aligned_coords.setMatrix(2 * number_Of_Atoms_In_Residue1, 3 * number_Of_Atoms_In_Residue1 - 1, 0, number_of_conformations - 1, X1);

						// ---> Residue2:
						Matrix residue_2_pair_aligned_coords = new Matrix(3 * number_Of_Atoms_In_Residue2, number_of_conformations);

						X1 = aligned_Residue_coords.getMatrix(0, number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1);
						Y1 = aligned_Residue_coords.getMatrix(number_of_atoms_in_pair, number_of_atoms_in_pair + number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1);
						Z1 = aligned_Residue_coords.getMatrix(2 * number_of_atoms_in_pair, 2 * number_of_atoms_in_pair + number_Of_Atoms_In_Residue2 - 1, 0,
								number_of_conformations - 1);

						residue_2_pair_aligned_coords.setMatrix(0, number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1, X1);
						residue_2_pair_aligned_coords.setMatrix(number_Of_Atoms_In_Residue2, 2 * number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1, X1);
						residue_2_pair_aligned_coords.setMatrix(2 * number_Of_Atoms_In_Residue2, 3 * number_Of_Atoms_In_Residue2 - 1, 0, number_of_conformations - 1, X1);

						// Get EigenResidue1
						JEDi_Get_Residue_Pair_PCA PCA1 = new JEDi_Get_Residue_Pair_PCA(residue_1_pair_aligned_coords, number_modes_Residue);
						PCA1.get_Residue_Pair_PCA();
						mean_centered_data = PCA1.getCentered_input_data();
						eigenvalues_COV = PCA1.getEigenvalues_COV();
						evectors_COV = PCA1.getTop_evectors_COV();

						// Need the PCs for each residue in the pair ---> to get the Q-COV matrix...
						Residue_Generalized_Coordinates_Pair = new Matrix(2 * number_modes_Residue, number_of_conformations);

						// Get PCs for residue 1
						JEDi_Get_Residue_Pair_PCs PCs1 = new JEDi_Get_Residue_Pair_PCs(mean_centered_data, evectors_COV, eigenvalues_COV, res_index);
						projections_COV = PCs1.get_PCs();
						Residue_Generalized_Coordinates_Pair.setMatrix(0, number_modes_Residue - 1, 0, number_of_conformations - 1, projections_COV.transpose());

						// Get EigenResidue2
						JEDi_Get_Residue_Pair_PCA PCA2 = new JEDi_Get_Residue_Pair_PCA(residue_2_pair_aligned_coords, number_modes_Residue);
						PCA2.get_Residue_Pair_PCA();
						mean_centered_data = PCA2.getCentered_input_data();
						eigenvalues_COV = PCA2.getEigenvalues_COV();
						evectors_COV = PCA2.getTop_evectors_COV();

						// Get the PCs for Residue 2
						JEDi_Get_Residue_Pair_PCs PCs2 = new JEDi_Get_Residue_Pair_PCs(mean_centered_data, evectors_COV, eigenvalues_COV, res_index);
						projections_COV = PCs2.get_PCs();
						Residue_Generalized_Coordinates_Pair.setMatrix(number_modes_Residue, 2 * number_modes_Residue - 1, 0, number_of_conformations - 1,
								projections_COV.transpose());

						// Get the U eigenvector matrix to compute the coupling score for the (i,j) Residue Pair
						JEDi_Get_Residue_Pair_Coupling rpc = new JEDi_Get_Residue_Pair_Coupling(Residue_Generalized_Coordinates_Pair, Residue_Pair_Index);
						rpc.get_Pairwise_Couping_Analysis();
						double score = rpc.getCoupling_score();
						// coupling_scores.set(i, j, score);
						coupling_scores.set(j, i, score); // this makes a lower triangular matrix
						coupling_scores.set(i, i, MAX_COUPLING_SCORE); // Set the diagonal elements to max score... (currently, it is 100)
					}
			}
		coupling_scores.set(number_residues - 1, number_residues - 1, MAX_COUPLING_SCORE); // Set the last diagonal element to max score... currently, it is 100

		String title = "Residue_Pair_Coupling_Score_Matrix_" + number_modes_Residue + "_DOF";
		String name = "Residue_Pair_Coupling_Score_Matrix_" + number_modes_Residue + "_DOF.txt";

		Matrix_IO.write_Matrix(coupling_scores, out_dir_pca + name, 9, 2);
		HeatMap hm1 = new HeatMap(coupling_scores, out_dir_pca, title, "Residue 1", "Residue 2");
		hm1.createCoupingScorePlot();
	}

	/* ******************************************************** PRIVATE METHODS ************************************************************** */

	private void get_Local_Transformed_Coords()
	{
		JEDi_Get_Transformed_Residue_Coordinates tf_coords = new JEDi_Get_Transformed_Residue_Coordinates(residue_PDB_coordinates, residue_REF_PDB_coordinates);

		aligned_Residue_REF_coords = tf_coords.get_Transformed_reference_coordinates();
		aligned_Residue_coords = tf_coords.get_SS_Transformed_coords();
		transformed_conformation_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		transformed_residue_rmsd_list = tf_coords.get_SS_Residue_RMSFs();

		Matrix stats = tf_coords.get_SS_coordinate_STATS();
		String name = "Residue_Coordinate_Stats.txt.bz2";
		String path = local_out_dir + name;

		if (verbose) Matrix_IO.write_BZ2_Matrix(stats, path, 12, 6);

		STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Residue_Variable_Means", "Atom Index", "Mean", stats, residueAtoms, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Residue_Variable_Variances", "Atom Index", "Variance", stats, residueAtoms, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Residue_Variable_Skews", "Atom Index", "Skew", stats, residueAtoms, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Residue_Variable_Kurtosis", "Atom Index", "Kurtosis", stats, residueAtoms, 4);

		if (Input_Parameters.doOutlierProcessing)
			{
				Matrix adjustments_per_variable;
				tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();


				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable = tf_coords.getAdjustments_per_variable_REMOVE_Outliers();
						name = "Adjustments_per_Variable_REMOVE_Outliers_MAD.txt.bz2";
						path = local_out_dir + name;
						if (verbose) Matrix_IO.write_BZ2_Matrix(adjustments_per_variable, path, 12, 0);
						STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Index", "Counts", adjustments_per_variable,
								residueAtoms, 1);

						adjustments_per_variable = tf_coords.getAdjustments_per_variable_SELECT_Outliers();
						name = "Adjustments_per_Variable_SELECT_Outliers_MAD.txt.bz2";
						path = local_out_dir + name;
						if (verbose) Matrix_IO.write_BZ2_Matrix(adjustments_per_variable, path, 12, 0);
						STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Index", "Counts", adjustments_per_variable,
								residueAtoms, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable = tf_coords.getAdjustments_per_variable_REMOVE_Outliers();
						name = "Adjustments_per_Variable_REMOVE_Outliers_Z.txt.bz2";
						path = local_out_dir + name;
						if (verbose) Matrix_IO.write_BZ2_Matrix(adjustments_per_variable, path, 12, 0);
						STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Index", "Counts", adjustments_per_variable,
								residueAtoms, 1);
						adjustments_per_variable = tf_coords.getAdjustments_per_variable_SELECT_Outliers();
						name = "Adjustments_per_Variable_SELECT_Outliers_Z.txt.bz2";
						path = local_out_dir + name;
						if (verbose) Matrix_IO.write_BZ2_Matrix(adjustments_per_variable, path, 12, 0);
						STATS_Plot.create_Variables_Stat_Line_Chart(local_out_dir, "Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Index", "Counts", adjustments_per_variable,
								residueAtoms, 1);
					}
			}
		RMSD_Plot.create_RMSF_Line_Chart(local_out_dir, "Residue_Atomic_RMSFs", "Atom Number", transformed_residue_rmsd_list, residueAtoms);
		RMSD_Plot.create_RMSD_XY_Chart(local_out_dir, "Residue_Conformation_RMSDs", "Conformation Number", transformed_conformation_rmsds);
	}

	private void get_Local_Residue_PCA()
	{
		JEDi_Get_Residue_PCA cr_pca = new JEDi_Get_Residue_PCA(aligned_Residue_coords, number_modes_Residue, type_alignment);

		cr_pca.set_Residue_Index(residue_Index_Orig);
		cr_pca.set_Output_Directory(local_out_dir);
		cr_pca.Do_COV_PCA();

		mean_centered_data = cr_pca.getCentered_input_data();

		cond_COV = cr_pca.get_cond_COV();
		trace_COV = cr_pca.get_trace_COV();
		det_COV = cr_pca.get_det_COV();
		rank_COV = cr_pca.get_rank_COV();
		eigenvalues_COV = cr_pca.get_Eigenvalues_COV();
		top_eigenvalues_COV = cr_pca.getTop_eigenvalues_COV();
		evectors_COV = cr_pca.get_Top_evectors_COV();
		pca_modes_COV = cr_pca.getPca_modes_COV();
		square_pca_modes_COV = cr_pca.getSquare_pca_modes_COV();
		weighted_square_pca_modes_COV = cr_pca.getWeighted_square_pca_modes_COV();
		weighted_pca_modes_COV = cr_pca.getWeighted_pca_modes_COV();
		pca_mode_mins_COV = cr_pca.get_pca_mode_mins_COV();
		pca_mode_maxes_COV = cr_pca.get_pca_mode_maxes_COV();

		MODES_Plot.createLineChart2Series(local_out_dir + "COV" + File.separatorChar, "Residue_" + res_index + "_Top_2_Weighted_Square_PCA_Modes_COV",
				weighted_square_pca_modes_COV, residueAtoms);

		if (doCORR)
			{
				cr_pca.Do_CORR_PCA();

				trace_CORR = cr_pca.get_trace_CORR();
				eigenvalues_CORR = cr_pca.get_Eigenvalues_CORR();
				evectors_CORR = cr_pca.get_Top_evectors_CORR();
				top_eigenvalues_CORR = cr_pca.get_Eigenvalues_CORR();
				pca_modes_CORR = cr_pca.getPca_modes_CORR();
				square_pca_modes_CORR = cr_pca.getSquare_pca_modes_CORR();
				weighted_square_pca_modes_CORR = cr_pca.getWeighted_square_pca_modes_CORR();
				weighted_pca_modes_CORR = cr_pca.getWeighted_pca_modes_CORR();
				pca_mode_mins_CORR = cr_pca.get_pca_mode_mins_CORR();
				pca_mode_maxes_CORR = cr_pca.get_pca_mode_maxes_CORR();

				MODES_Plot.createLineChart2Series(local_out_dir + "CORR" + File.separatorChar, "Residue_" + res_index + "_Top_2_Weighted_Square_PCA_Modes_CORR",
						weighted_square_pca_modes_CORR, residueAtoms);

				if (doSPARSIFY)
					{
						cr_pca.Do_CORR_SPARSE_PCA();

						trace_CORR_SPARSE = cr_pca.getTrace_CORR_SPARSE();
						eigenvalues_CORR_SPARSE = cr_pca.getEigenvalues_CORR_SPARSE();
						evectors_CORR_SPARSE = cr_pca.get_Top_evectors_CORR_SPARSE();
						top_eigenvalues_CORR_SPARSE = cr_pca.getTop_eigenvalues_CORR_SPARSE();
						pca_modes_CORR_SPARSE = cr_pca.getPca_modes_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE = cr_pca.getSquare_pca_modes_CORR_SPARSE();
						weighted_square_pca_modes_CORR_SPARSE = cr_pca.getWeighted_square_pca_modes_CORR_SPARSE();
						weighted_pca_modes_CORR_SPARSE = cr_pca.getWeighted_pca_modes_CORR_SPARSE();
						pca_mode_mins_CORR_SPARSE = cr_pca.get_pca_mode_mins_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE = cr_pca.get_pca_mode_maxes_CORR_SPARSE();

						String out = local_out_dir + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						create_Directory(out);
						MODES_Plot.createLineChart2Series(out, "Residue_" + res_index + "_Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", weighted_square_pca_modes_CORR_SPARSE,
								residueAtoms);
					}
			}
		if (doPCORR)
			{
				cr_pca.Do_PCORR_PCA();

				trace_PCORR = cr_pca.get_trace_PCORR();
				eigenvalues_PCORR = cr_pca.get_Eigenvalues_PCORR();
				evectors_PCORR = cr_pca.getTop_evectors_PCORR();
				top_eigenvalues_PCORR = cr_pca.get_Eigenvalues_PCORR();
				pca_modes_PCORR = cr_pca.getPca_modes_PCORR();
				square_pca_modes_PCORR = cr_pca.getSquare_pca_modes_PCORR();
				weighted_square_pca_modes_PCORR = cr_pca.getWeighted_square_pca_modes_PCORR();
				weighted_pca_modes_PCORR = cr_pca.getWeighted_pca_modes_PCORR();
				pca_mode_mins_PCORR = cr_pca.get_pca_mode_mins_PCORR();
				pca_mode_maxes_PCORR = cr_pca.get_pca_mode_maxes_PCORR();

				MODES_Plot.createLineChart2Series(local_out_dir + "PCORR" + File.separatorChar, "Residue_" + res_index + "_Top_2_Weighted_Square_PCA_Modes_PCORR",
						weighted_square_pca_modes_PCORR, residueAtoms);

				if (doSPARSIFY)
					{
						cr_pca.Do_PCORR_SPARSE_PCA();

						trace_PCORR_SPARSE = cr_pca.getTrace_PCORR_SPARSE();
						eigenvalues_PCORR_SPARSE = cr_pca.getEigenvalues_PCORR_SPARSE();
						evectors_PCORR_SPARSE = cr_pca.get_Top_evectors_PCORR_SPARSE();
						top_eigenvalues_PCORR_SPARSE = cr_pca.getTop_eigenvalues_PCORR_SPARSE();
						pca_modes_PCORR_SPARSE = cr_pca.getPca_modes_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE = cr_pca.getSquare_pca_modes_PCORR_SPARSE();
						weighted_square_pca_modes_PCORR_SPARSE = cr_pca.getWeighted_square_pca_modes_PCORR_SPARSE();
						weighted_pca_modes_PCORR_SPARSE = cr_pca.getWeighted_pca_modes_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE = cr_pca.get_pca_mode_mins_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE = cr_pca.get_pca_mode_maxes_PCORR_SPARSE();

						String out = local_out_dir + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						create_Directory(out);
						MODES_Plot.createLineChart2Series(out, "Residue_" + res_index + "_Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", weighted_square_pca_modes_PCORR_SPARSE,
								residueAtoms);
					}
			}
	}

	private void get_Global_Residue_PCA(Matrix data, Matrix adj_data)
	{
		JEDi_Get_Residue_PCA cr_pca = new JEDi_Get_Residue_PCA(adj_data, number_modes_Residue, type_alignment); // Using data with Outlier Processing for the PCA...
			{
				cr_pca.set_Residue_Index(residue_Index_Orig);
				cr_pca.Do_COV_PCA();

				Row_Center_Data cd = new Row_Center_Data(data); // Need to project the original aligned data for the PCs...
				mean_centered_data = cd.get_row_centered_data();

				eigenvalues_COV = cr_pca.get_Eigenvalues_COV();
				evectors_COV = cr_pca.get_Top_evectors_COV();

				/* ********** For the Hierarchical PCA **************** */
				Residue_Eigenvectors_Global.add(evectors_COV);
			}
	}

	private void get_Local_Residue_DVPs()
	{
		JEDi_Get_Residue_DVPs dvps_cov = new JEDi_Get_Residue_DVPs(aligned_Residue_coords, aligned_Residue_REF_coords, evectors_COV, eigenvalues_COV, Q);

		dvps_cov.set_Residue_Index(residue_Index_Orig);
		dvps_cov.set_Output_Directory(local_out_dir + Q + File.separatorChar);

		delta_vectors = dvps_cov.get_Cartesian_DV_Series();
		projections_COV = dvps_cov.getProjections();
		normed_projections_COV = dvps_cov.getNormed_projections();
		weighted_projections_COV = dvps_cov.getWeighted_projections();
		weighted_normed_projections_COV = dvps_cov.getWeighted_normed_projections();

		if (doCORR)
			{
				JEDi_Get_Residue_DVPs dvps_corr = new JEDi_Get_Residue_DVPs(aligned_Residue_coords, aligned_Residue_REF_coords, evectors_CORR, eigenvalues_CORR, R);

				dvps_corr.set_Residue_Index(residue_Index_Orig);
				dvps_corr.set_Output_Directory(local_out_dir + R + File.separatorChar);

				delta_vectors = dvps_corr.get_Cartesian_DV_Series();
				projections_CORR = dvps_corr.getProjections();
				normed_projections_CORR = dvps_corr.getNormed_projections();
				weighted_projections_CORR = dvps_corr.getWeighted_projections();
				weighted_normed_projections_CORR = dvps_corr.getWeighted_normed_projections();

				if (doSPARSIFY)
					{
						String out = local_out_dir + R + File.separatorChar + "sparse" + File.separatorChar;
						dvps_corr = new JEDi_Get_Residue_DVPs(aligned_Residue_coords, aligned_Residue_REF_coords, evectors_CORR_SPARSE, eigenvalues_CORR_SPARSE, R);

						dvps_corr.set_Residue_Index(residue_Index_Orig);
						dvps_corr.set_Output_Directory(out);

						delta_vectors = dvps_corr.get_Cartesian_DV_Series();
						projections_CORR_SPARSE = dvps_corr.getProjections();
						normed_projections_CORR_SPARSE = dvps_corr.getNormed_projections();
						weighted_projections_CORR_SPARSE = dvps_corr.getWeighted_projections();
						weighted_normed_projections_CORR_SPARSE = dvps_corr.getWeighted_normed_projections();
					}
			}

		if (doPCORR)
			{
				JEDi_Get_Residue_DVPs dvps_pcorr = new JEDi_Get_Residue_DVPs(aligned_Residue_coords, aligned_Residue_REF_coords, evectors_PCORR, eigenvalues_PCORR, P);

				dvps_pcorr.set_Residue_Index(residue_Index_Orig);
				dvps_pcorr.set_Output_Directory(local_out_dir + P + File.separatorChar);

				delta_vectors = dvps_pcorr.get_Cartesian_DV_Series();
				projections_PCORR = dvps_pcorr.getProjections();
				normed_projections_PCORR = dvps_pcorr.getNormed_projections();
				weighted_projections_PCORR = dvps_pcorr.getWeighted_projections();
				weighted_normed_projections_PCORR = dvps_pcorr.getWeighted_normed_projections();

				if (doSPARSIFY)
					{
						String out = local_out_dir + P + File.separatorChar + "sparse" + File.separatorChar;
						dvps_pcorr = new JEDi_Get_Residue_DVPs(aligned_Residue_coords, aligned_Residue_REF_coords, evectors_PCORR_SPARSE, eigenvalues_PCORR_SPARSE, P);

						dvps_pcorr.set_Residue_Index(residue_Index_Orig);
						dvps_pcorr.set_Output_Directory(out);

						delta_vectors = dvps_pcorr.get_Cartesian_DV_Series();
						projections_PCORR_SPARSE = dvps_pcorr.getProjections();
						normed_projections_PCORR_SPARSE = dvps_pcorr.getNormed_projections();
						weighted_projections_PCORR_SPARSE = dvps_pcorr.getWeighted_projections();
						weighted_normed_projections_PCORR_SPARSE = dvps_pcorr.getWeighted_normed_projections();
					}
			}
	}

	private void get_Global_Residue_PCs()
	{
		JEDi_Get_Residue_PCs pcs_cov = new JEDi_Get_Residue_PCs(mean_centered_data, evectors_COV, eigenvalues_COV);
		pcs_cov.get_PCs();

		Residues_Centered_Data_Global.add(mean_centered_data);
		projections_COV = pcs_cov.getProjections();
		normed_projections_COV = pcs_cov.getNormed_projections();
		weighted_projections_COV = pcs_cov.getWeighted_projections();
		weighted_normed_projections_COV = pcs_cov.getWeighted_normed_projections();

		int row_index1 = residue_Index * number_modes_Residue;
		int row_index2 = row_index1 + number_modes_Residue - 1;
		int col_index1 = 0;
		int col_index2 = projections_COV.transpose().getColumnDimension() - 1;

		Residue_Generalized_Coordinates_Global.setMatrix(row_index1, row_index2, col_index1, col_index2, projections_COV.transpose());
	}

	private void get_Local_Residue_Mode_Visualization(String outdir)
	{
		JEDi_Get_PCA_Mode_Vizualization mvQ = new JEDi_Get_PCA_Mode_Vizualization(residueAtoms, eigenvalues_COV, evectors_COV, square_pca_modes_COV, pca_mode_maxes_COV,
				pca_mode_mins_COV, number_of_modes_viz, ("Residue_PCA" + File.separatorChar + type_alignment), Q);

		mvQ.set_Output_Directory(outdir + Q + File.separatorChar);
		mvQ.get_Essential_Visualization_All_Atom();
		if (do_INDIVIDUAL_VIZ) mvQ.get_Mode_Visualizations_All_Atom();
		mvQ = null;

		if (doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization mvR = new JEDi_Get_PCA_Mode_Vizualization(residueAtoms, eigenvalues_CORR, evectors_CORR, square_pca_modes_CORR, pca_mode_maxes_CORR,
						pca_mode_mins_CORR, number_of_modes_viz, ("Residue_PCA" + File.separatorChar + type_alignment), R);

				mvR.set_Output_Directory(outdir + R + File.separatorChar);

				mvR.get_Essential_Visualization_All_Atom();
				if (do_INDIVIDUAL_VIZ) mvR.get_Mode_Visualizations_All_Atom();
				mvR = null;

				if (doSPARSIFY)
					{
						mvR = new JEDi_Get_PCA_Mode_Vizualization(residueAtoms, eigenvalues_CORR_SPARSE, evectors_CORR_SPARSE, square_pca_modes_CORR_SPARSE,
								pca_mode_maxes_CORR_SPARSE, pca_mode_mins_CORR_SPARSE, number_of_modes_viz, ("Residue_PCA" + File.separatorChar + type_alignment), R);

						mvR.set_Output_Directory(outdir + R + File.separatorChar + "sparse" + File.separatorChar);

						mvR.get_Essential_Visualization_All_Atom();
						if (do_INDIVIDUAL_VIZ) mvR.get_Mode_Visualizations_All_Atom();
						mvR = null;
					}
			}

		if (doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization mvP = new JEDi_Get_PCA_Mode_Vizualization(residueAtoms, eigenvalues_PCORR, evectors_PCORR, square_pca_modes_PCORR,
						pca_mode_maxes_PCORR, pca_mode_mins_PCORR, number_of_modes_viz, ("Residue_PCA" + File.separatorChar + type_alignment), P);

				mvP.set_Output_Directory(outdir + P + File.separatorChar);

				mvP.get_Essential_Visualization_All_Atom();
				if (do_INDIVIDUAL_VIZ) mvP.get_Mode_Visualizations_All_Atom();
				mvP = null;

				if (doSPARSIFY)
					{
						mvP = new JEDi_Get_PCA_Mode_Vizualization(residueAtoms, eigenvalues_PCORR_SPARSE, evectors_PCORR_SPARSE, square_pca_modes_PCORR_SPARSE,
								pca_mode_maxes_PCORR_SPARSE, pca_mode_mins_PCORR_SPARSE, number_of_modes_viz, ("Residue_PCA" + File.separatorChar + type_alignment), P);

						mvP.set_Output_Directory(outdir + P + File.separatorChar + "sparse" + File.separatorChar);

						mvP.get_Essential_Visualization_All_Atom();
						if (do_INDIVIDUAL_VIZ) mvP.get_Mode_Visualizations_All_Atom();
						mvP = null;
					}
			}
	}

	private void do_FES(String outdir)
	{
		JEDi_Get_FES fes = new JEDi_Get_FES(normed_projections_COV, 0, 1, number_of_conformations, 0, 0);
		fes.set_Out_dir(outdir + Q + File.separatorChar);
		fes.get_FES();
		fes.write_FES_Log();

		if (doCORR)
			{
				fes = new JEDi_Get_FES(normed_projections_CORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(outdir + R + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();

				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_CORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(outdir + R + File.separatorChar + "SPARSE" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}
		if (doPCORR)
			{

				fes = new JEDi_Get_FES(normed_projections_PCORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(outdir + P + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();

				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_PCORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(outdir + P + File.separatorChar + "SPARSE" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}
	}

	private void do_KPCA(String outdir)
	{
		JEDi_Do_Kernel_PCA kpca = new JEDi_Do_Kernel_PCA(normed_projections_COV, outdir + Q + File.separatorChar);
		kpca.kPCA_Driver();

		if (doCORR)
			{
				kpca = new JEDi_Do_Kernel_PCA(normed_projections_CORR, outdir + R + File.separatorChar);
				kpca.kPCA_Driver();

				if (doSPARSIFY)
					{
						kpca = new JEDi_Do_Kernel_PCA(normed_projections_CORR_SPARSE, outdir + R + File.separatorChar + "SPARSE" + File.separatorChar);
						kpca.kPCA_Driver();
					}
			}

		if (doPCORR)
			{
				kpca = new JEDi_Do_Kernel_PCA(normed_projections_PCORR, outdir + P + File.separatorChar);
				kpca.kPCA_Driver();

				if (doSPARSIFY)
					{
						kpca = new JEDi_Do_Kernel_PCA(normed_projections_PCORR_SPARSE, outdir + P + File.separatorChar + "SPARSE" + File.separatorChar);
						kpca.kPCA_Driver();
					}
			}
	}

	private void do_SSA(String outdir)
	{

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_Subspace_Analysis jssa1_C_AA = new JEDi_Get_Subspace_Analysis(cRes_AA, "COV_vs_CORR", evectors_CORR, evectors_COV);

				jssa1_C_AA.setOut_dir(outdir);
				jssa1_C_AA.get_SSA();
				jssa1_C_AA.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_Subspace_Analysis jssa2_C_AA = new JEDi_Get_Subspace_Analysis(cRes_AA, "COV_vs_PCORR", evectors_PCORR, evectors_COV);

				jssa2_C_AA.setOut_dir(outdir);
				jssa2_C_AA.get_SSA();
				jssa2_C_AA.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				JEDi_Get_Subspace_Analysis jssa3_C_AA = new JEDi_Get_Subspace_Analysis(cRes_AA, "CORR_vs_PCORR", evectors_PCORR, evectors_CORR);

				jssa3_C_AA.setOut_dir(outdir);
				jssa3_C_AA.get_SSA();
				jssa3_C_AA.get_FSSA_Iterated();
			}

		if (doSPARSIFY)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cRes_AA, "CORR_vs_CORR_SPARSE", evectors_CORR_SPARSE, evectors_CORR);

				ssa.setOut_dir(outdir);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();

				ssa = new JEDi_Get_Subspace_Analysis(cRes_AA, "PCORR_vs_PCORR_SPARSE", evectors_PCORR_SPARSE, evectors_PCORR);

				ssa.setOut_dir(outdir);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}


	}

	private void create_Directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
	}

	/* ************************************************************ SETTERS ***************************************************************************** */

	public void set_out_dir_pca(String out)
	{
		this.out_dir_pca = out;
		create_Directory(out_dir_pca);
	}

	public void set_out_dir_ssa(String out)
	{
		this.out_dir_ssa = out;
		create_Directory(out_dir_ssa);
	}

	public void set_out_dir_fes(String out)
	{
		this.out_dir_fes = out;
		create_Directory(out_dir_fes);
	}

	public void set_out_dir_kpca(String out)
	{
		this.out_dir_kpca = out;
		create_Directory(out_dir_kpca);
	}

	public void set_out_dir_viz(String out)
	{
		this.out_dir_viz = out;
		create_Directory(out_dir_viz);
	}

//---------------------------------------------------------------------------------------------------------------------------------

	public void setDo_CORR(boolean doCORR)
	{
		this.doCORR = doCORR;
	}

	public void setDo_PCORR(boolean doPCORR)
	{
		this.doPCORR = doPCORR;
	}

	public void setDoReduce(boolean doReduce)
	{
		this.doReduce = doReduce;
	}

	public void setDo_ESSENTIAL_VIZ(boolean do_essential)
	{
		this.do_ESSENTIAL_VIZ = do_essential;
	}

	public void setDo_INDIVIDUAL_VIZ(boolean do_individual)
	{
		this.do_INDIVIDUAL_VIZ = do_individual;
	}

	public void setDoHeavyAtoms(boolean doHeavyAtoms)
	{
		this.doHeavyAtoms = doHeavyAtoms;
	}

	public void setDoAllAtoms(boolean doAllAtoms)
	{
		this.doAllAtoms = doAllAtoms;
	}

	public void setDoFES(boolean doFES)
	{
		this.doFES = doFES;
	}

	public void setDoKPCA(boolean doKPCA)
	{
		this.doKPCA = doKPCA;
	}

	/* ************************************************************ GETTERS ***************************************************************************** */

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

	public double get_Trace_PCORR()
	{
		return trace_PCORR;
	}

	public double get_Cond_COV()
	{
		return cond_COV;
	}

	public double get_Det_COV()
	{
		return det_COV;
	}

	public double get_Rank_COV()
	{
		return rank_COV;
	}

	public List<Double> getTop_cartesian_eigenvalues_COV()
	{
		return top_eigenvalues_COV;
	}

	public List<Double> getTop_cartesian_eigenvalues_CORR()
	{
		return top_eigenvalues_CORR;
	}

	public double[] getPca_mode_max_COV()
	{
		return pca_mode_maxes_COV;
	}

	public double[] getPca_mode_min_COV()
	{
		return pca_mode_mins_COV;
	}

	public double[] getPca_mode_max_CORR()
	{
		return pca_mode_maxes_CORR;
	}

	public double[] getPca_mode_min_CORR()
	{
		return pca_mode_mins_CORR;
	}

	public Matrix getCov()
	{
		return cov;
	}

	public Matrix getTop_cartesian_evectors_COV()
	{
		return evectors_COV;
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

	public Matrix getTop_cartesian_evectors_CORR()
	{
		return evectors_CORR;
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

	public int getNumber_of_residues()
	{
		return number_residues;
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

	public double getRank_COV()
	{
		return rank_COV;
	}

	public List<Double> getTop_cartesian_eigenvalues_PCORR()
	{
		return top_eigenvalues_PCORR;
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
		return evectors_PCORR;
	}

	public double[] getPca_mode_max_PCORR()
	{
		return pca_mode_maxes_PCORR;
	}

	public double[] getPca_mode_min_PCORR()
	{
		return pca_mode_mins_PCORR;
	}

	public Matrix getSquare_pca_modes_PCORR()
	{
		return square_pca_modes_PCORR;
	}

	public Matrix get_Residue_Generalized_Coordinates_Global()
	{
		return Residue_Generalized_Coordinates_Global;
	}

	public List<Matrix> get_Residue_Eigenvectors()
	{
		return Residue_Eigenvectors_Global;
	}

	public List<Matrix> getResidues_Centered_Data()
	{
		return Residues_Centered_Data_Global;
	}
}
