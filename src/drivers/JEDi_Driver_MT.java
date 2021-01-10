package drivers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import Jama.Matrix;
import jedi.JEDi_Do_Cartesian;
import jedi.JEDi_Do_Dist_Pairs;
import jedi.JEDi_Do_Hierarchical;
import jedi.JEDi_Do_PDB_Processing;
import jedi.JEDi_Do_Residue;
import jedi.JEDi_Get_Coordinates_from_Matrix;
import jedi.JEDi_Get_Distances_for_Atom_Pairs;
import jedi.JEDi_Get_PCA_Mode_Vizualization;
import jedi.JEDi_Get_Subspace_Analysis;
import jedi.JEDi_Get_Transformed_Coordinates;
import support.Atom;
import support.Atom_ID_Pair;
import support.KMO_MSA;
import support.Outlier_Processing;
import support.Residue_ID_Pair;
import support.Select_Variables_by_Statistical_Threshold;
import supportIO.DateUtils;
import supportIO.Input_Parameters;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportIO.Output_Control;
import supportIO.PDB_IO;
import supportIO.Write_DP_Stats;
import supportPlot.MODES_Plot;
import supportPlot.PC_Plot;
import supportPlot.Plot_Line_Chart;
import supportPlot.RMSD_Plot;
import supportPlot.STATS_Plot;

/**
 * JEDi class JEDi_Driver_MT: Driver program for running JEDi Software with Multiple Threads.
 * 
 * The default Input file is "JEDi_Parameters.txt"
 * 
 * Optionally, the first passed program argument can be used to specify a custom input file.
 * 
 * Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */
public class JEDi_Driver_MT
{
	static int number_of_residues_REF, number_of_residues_AA_SS, number_of_residues_BB_SS, number_of_residues_HA_SS, number_of_residues_CA_SS, number_of_residues_local_SS,
			number_of_residues_pairs_SS, number_of_residues_HAA_SS, number_of_residues_HHA_SS;

	static int number_of_atoms_REF, number_of_heavy_atoms_REF, number_of_atoms_local_SS, number_of_atoms_residue_pairs_SS, number_of_atoms_hierarchical_AA_SS,
			number_of_atoms_hierarchical_HA_SS, number_of_atoms_AA_SS, number_of_atoms_BB_SS, number_of_atoms_HA_SS, number_of_atoms_AL_SS, number_of_atom_pairs_DP_SS;

	static int number_conformations, ROWS_AA, ROWS_BB, ROWS_CA, ROWS_HA, COLS, number_Of_Input_Lines;
	static int rank_AA, rank_HA, rank_BB, rank_CA, rank_HAA, rank_HHA, rank_d, rank_AL;

	static double trace_d, cond_d, det_d, shrinkage_d;
	static double shrinkage_AA, trace_AA, cond_AA, det_AA;
	static double shrinkage_AL, trace_AL, cond_AL, det_AL;
	static double shrinkage_HA, trace_HA, cond_HA, det_HA;
	static double shrinkage_BB, trace_BB, cond_BB, det_BB;
	static double shrinkage_CA, trace_CA, cond_CA, det_CA;
	static double shrinkage_HAA, trace_HAA, cond_HAA, det_HAA;
	static double shrinkage_HHA, trace_HHA, cond_HHA, det_HHA;
	static double KMO_AA, KMO_HA, KMO_BB, KMO_CA, KMO_DP, KMO_AL;

	static String key, value, path, name, line, date, type, model;

	static final String input_file = "JEDi_Parameters.txt", delim = "=", Q = "COV", R = "CORR", P = "PCORR", S = "sparse", PP = "Pre_Process", AA = "All_Atom", HA = "Heavy_Atom",
			BB = "BackBone", CA = "Alpha_Carbon", cAL = "Atom_List_PCA", cAA = "All_Atom_PCA", cBB = "Backbone_PCA", cAC = "Alpha_Carbon_PCA", cHA = "Heavy_Atom_PCA",
			Res_AA = "Residue_All_Atom_PCA", Res_PAIRS = "Residue_Pair_Analysis", HAA = "Hierarchical_All_Atom_PCA", HHA = "Hierarchical_Heavy_Atom_PCA", DP = "Distance_Pair_PCA";

	static double[] atomic_distance_means, atomic_distance_std_devs;
	static double[] convoluted_pca_mode_mins_AA, convoluted_pca_mode_maxes_AA, convoluted_pca_mode_mins_HA, convoluted_pca_mode_maxes_HA;

	static double[] pca_mode_mins_COV_AA, pca_mode_maxes_COV_AA, pca_mode_mins_CORR_AA, pca_mode_maxes_CORR_AA, pca_mode_mins_PCORR_AA, pca_mode_maxes_PCORR_AA;
	static double[] pca_mode_mins_COV_AL, pca_mode_maxes_COV_AL, pca_mode_mins_CORR_AL, pca_mode_maxes_CORR_AL, pca_mode_mins_PCORR_AL, pca_mode_maxes_PCORR_AL;
	static double[] pca_mode_mins_COV_BB, pca_mode_maxes_COV_BB, pca_mode_mins_CORR_BB, pca_mode_maxes_CORR_BB, pca_mode_mins_PCORR_BB, pca_mode_maxes_PCORR_BB;
	static double[] pca_mode_mins_COV_CA, pca_mode_maxes_COV_CA, pca_mode_mins_CORR_CA, pca_mode_maxes_CORR_CA, pca_mode_mins_PCORR_CA, pca_mode_maxes_PCORR_CA;
	static double[] pca_mode_mins_COV_HA, pca_mode_maxes_COV_HA, pca_mode_mins_CORR_HA, pca_mode_maxes_CORR_HA, pca_mode_mins_PCORR_HA, pca_mode_maxes_PCORR_HA;

	static double[] pca_mode_mins_CORR_SPARSE_AA, pca_mode_maxes_CORR_SPARSE_AA, pca_mode_mins_PCORR_SPARSE_AA, pca_mode_maxes_PCORR_SPARSE_AA;
	static double[] pca_mode_mins_CORR_SPARSE_AL, pca_mode_maxes_CORR_SPARSE_AL, pca_mode_mins_PCORR_SPARSE_AL, pca_mode_maxes_PCORR_SPARSE_AL;
	static double[] pca_mode_mins_CORR_SPARSE_BB, pca_mode_maxes_CORR_SPARSE_BB, pca_mode_mins_PCORR_SPARSE_BB, pca_mode_maxes_PCORR_SPARSE_BB;
	static double[] pca_mode_mins_CORR_SPARSE_CA, pca_mode_maxes_CORR_SPARSE_CA, pca_mode_mins_PCORR_SPARSE_CA, pca_mode_maxes_PCORR_SPARSE_CA;
	static double[] pca_mode_mins_CORR_SPARSE_HA, pca_mode_maxes_CORR_SPARSE_HA, pca_mode_mins_PCORR_SPARSE_HA, pca_mode_maxes_PCORR_SPARSE_HA;

	static List<Integer> atom_list, atom_list_original, atom_list_dp1, atom_list_dp2, atom_list_dp_original1, atom_list_dp_original2, atoms_read, heavy_atoms_read,
			alpha_carbons_read, numbers_Of_Atoms_in_Residues, numbers_Of_Heavy_Atoms_in_Residues, residues_read, SS_Atom_Numbers_AA, SS_Atom_Numbers_BB, SS_Atom_Numbers_HA,
			SS_Atom_Numbers_CA, SS_Atom_Numbers_local, SS_Atom_Numbers_pairs, SS_Atom_Numbers_Hierarchical_AA, SS_Atom_Numbers_Hierarchical_HA;

	static List<Double> aligned_atomic_RMSFs_AL, aligned_atomic_RMSFs_AA, aligned_atomic_RMSFs_HA, aligned_atomic_RMSFs_BB, aligned_atomic_RMSFs_CA,
			aligned_atomic_RMSFs_Hierarchical_AA, aligned_atomic_RMSFs_Hierarchical_HA, aligned_atomic_RMSFs_local;

	static List<Double> top_eigenvalues_COV_AA, top_eigenvalues_CORR_AA, top_eigenvalues_PCORR_AA, top_eigenvalues_CORR_SPARSE_AA, top_eigenvalues_PCORR_SPARSE_AA;
	static List<Double> top_eigenvalues_COV_HA, top_eigenvalues_CORR_HA, top_eigenvalues_PCORR_HA, top_eigenvalues_CORR_SPARSE_HA, top_eigenvalues_PCORR_SPARSE_HA;
	static List<Double> top_eigenvalues_COV_CA, top_eigenvalues_CORR_CA, top_eigenvalues_PCORR_CA, top_eigenvalues_CORR_SPARSE_CA, top_eigenvalues_PCORR_SPARSE_CA;
	static List<Double> top_eigenvalues_COV_BB, top_eigenvalues_CORR_BB, top_eigenvalues_PCORR_BB, top_eigenvalues_CORR_SPARSE_BB, top_eigenvalues_PCORR_SPARSE_BB;
	static List<Double> top_eigenvalues_COV_AL, top_eigenvalues_CORR_AL, top_eigenvalues_PCORR_AL, top_eigenvalues_CORR_SPARSE_AL, top_eigenvalues_PCORR_SPARSE_AL;
	static List<Double> top_eigenvalues_COV_DP, top_eigenvalues_CORR_DP, top_eigenvalues_PCORR_DP;
	static List<Double> top_eigenvalues_HAA, top_eigenvalues_HHA;

	static List<String> lines, pdb_file_names, chain_idents0, chain_idents1, chain_idents2, chain_ids_read;

	static List<Matrix> Residue_Delta_Vectors_Global_AA, Residue_Delta_Vectors_Global_HA, Residue_Centered_Coordinates_Global_AA, Residue_Centered_Coordinates_Global_HA,
			Residue_Eigenvectors_Global_AA, Residue_Eigenvectors_Global_HA, Residue_PCs_Global_AA, Residue_PCs_Global_HA;

	static Matrix original_reference_PDB_coordinates_AA, original_reference_PDB_coordinates_BB, original_reference_PDB_coordinates_CA, original_reference_PDB_coordinates_HA,
			original_PDB_coordinates_AA, aligned_subset_REF_PDB_coordinates_hierarchical_AA, aligned_subset_REF_PDB_coordinates_hierarchical_HA;

	static Matrix subset_PDB_coordinates_AA, subset_PDB_coordinates_BB, subset_PDB_coordinates_CA, subset_PDB_coordinates_HA, subset_PDB_coordinates_local,
			subset_PDB_coordinates_pairs, subset_PDB_coordinates_hierarchical_AA, subset_PDB_coordinates_hierarchical_HA, subset_PDB_coordinates_AL;

	static Matrix aligned_subset_PDB_coordinates_AA, aligned_subset_PDB_coordinates_BB, aligned_subset_PDB_coordinates_CA, aligned_subset_PDB_coordinates_HA,
			aligned_subset_PDB_coordinates_hierarchical_AA, aligned_subset_PDB_coordinates_hierarchical_HA, aligned_subset_PDB_coordinates_AL, aligned_subset_PDB_coordinates_local;

	static Matrix subset_PDB_coordinates_Outliers_REMOVED_AA, subset_PDB_coordinates_Outliers_REMOVED_BB, subset_PDB_coordinates_Outliers_REMOVED_CA,
			subset_PDB_coordinates_Outliers_REMOVED_HA, subset_PDB_coordinates_hierarchical_Outliers_REMOVED_AA, subset_PDB_coordinates_hierarchical_Outliers_REMOVED_HA,
			subset_PDB_coordinates_Outliers_REMOVED_AL, aligned_subset_PDB_coordinates_Outliers_REMOVED_local;

	static Matrix subset_PDB_coordinates_Outliers_SELECTED_AA, subset_PDB_coordinates_Outliers_SELECTED_BB, subset_PDB_coordinates_Outliers_SELECTED_CA,
			subset_PDB_coordinates_Outliers_SELECTED_HA, subset_PDB_coordinates_hierarchical_Outliers_SELECTED_AA, subset_PDB_coordinates_hierarchical_Outliers_SELECTED_HA,
			subset_PDB_coordinates_Outliers_SELECTED_AL, aligned_subset_PDB_coordinates_Outliers_SELECTED_local;

	static Matrix aligned_subset_REF_PDB_coordinates_AA, aligned_subset_REF_PDB_coordinates_BB, aligned_subset_REF_PDB_coordinates_CA, aligned_subset_REF_PDB_coordinates_HA,
			aligned_subset_REF_PDB_coordinates_AL, aligned_subset_REF_PDB_coordinates_local;

	static Matrix adjustments_per_variable_REMOVE_AA, adjustments_per_variable_REMOVE_HA, adjustments_per_variable_REMOVE_BB, adjustments_per_variable_REMOVE_CA,
			adjustments_per_variable_REMOVE_HAA, adjustments_per_variable_REMOVE_HHA, adjustments_per_variable_REMOVE_AL;

	static Matrix adjustments_per_variable_SELECT_AA, adjustments_per_variable_SELECT_HA, adjustments_per_variable_SELECT_BB, adjustments_per_variable_SELECT_CA,
			adjustments_per_variable_SELECT_HAA, adjustments_per_variable_SELECT_HHA, adjustments_per_variable_SELECT_AL;

	static Matrix adjustments_per_variable_AA, adjustments_per_variable_HA, adjustments_per_variable_BB, adjustments_per_variable_CA, adjustments_per_variable_HAA,
			adjustments_per_variable_HHA, adjustments_per_variable_AL;

	static Matrix subset_reference_PDB_coordinates_CA, subset_reference_PDB_coordinates_AA, subset_reference_PDB_coordinates_BB, subset_reference_PDB_coordinates_HA,
			subset_reference_PDB_coordinates_local, subset_reference_PDB_coordinates_pairs, subset_reference_PDB_coordinates_HAA, subset_reference_PDB_coordinates_HHA,
			subset_reference_PDB_coordinates_AL;

	static Matrix distances, distances_Outliers_REMOVED, distances_Outliers_SELECTED, reference_distances;

	static Matrix square_pca_modes_COV_AA, square_pca_modes_CORR_AA, square_pca_modes_PCORR_AA, square_pca_modes_CORR_SPARSE_AA, square_pca_modes_PCORR_SPARSE_AA;
	static Matrix square_pca_modes_COV_HA, square_pca_modes_CORR_HA, square_pca_modes_PCORR_HA, square_pca_modes_CORR_SPARSE_HA, square_pca_modes_PCORR_SPARSE_HA;
	static Matrix square_pca_modes_COV_BB, square_pca_modes_CORR_BB, square_pca_modes_PCORR_BB, square_pca_modes_CORR_SPARSE_BB, square_pca_modes_PCORR_SPARSE_BB;
	static Matrix square_pca_modes_COV_CA, square_pca_modes_CORR_CA, square_pca_modes_PCORR_CA, square_pca_modes_CORR_SPARSE_CA, square_pca_modes_PCORR_SPARSE_CA;
	static Matrix square_pca_modes_COV_AL, square_pca_modes_CORR_AL, square_pca_modes_PCORR_AL, square_pca_modes_CORR_SPARSE_AL, square_pca_modes_PCORR_SPARSE_AL;
	static Matrix square_pca_modes_HAA, square_pca_modes_HHA;
	static Matrix weighted_square_pca_modes_COV_dist, weighted_square_pca_modes_CORR_dist, weighted_square_pca_modes_PCORR_dist;

	static Matrix top_evectors_COV_AA, top_evectors_CORR_AA, top_evectors_PCORR_AA, top_evectors_CORR_SPARSE_AA, top_evectors_PCORR_SPARSE_AA;
	static Matrix top_evectors_COV_AL, top_evectors_CORR_AL, top_evectors_PCORR_AL, top_evectors_CORR_SPARSE_AL, top_evectors_PCORR_SPARSE_AL;
	static Matrix top_evectors_COV_HA, top_evectors_CORR_HA, top_evectors_PCORR_HA, top_evectors_CORR_SPARSE_HA, top_evectors_PCORR_SPARSE_HA;
	static Matrix top_evectors_COV_BB, top_evectors_CORR_BB, top_evectors_PCORR_BB, top_evectors_CORR_SPARSE_BB, top_evectors_PCORR_SPARSE_BB;
	static Matrix top_evectors_COV_CA, top_evectors_CORR_CA, top_evectors_PCORR_CA, top_evectors_CORR_SPARSE_CA, top_evectors_PCORR_SPARSE_CA;

	static Matrix top_evectors_HAA, top_evectors_HHA;

	static Matrix top_distance_evectors_COV, top_distance_evectors_CORR, top_distance_evectors_PCORR;

	static Matrix HAA_Evects, HAA_Evects_OR, HAA_Evects_OS, HHA_Evects, HHA_Evects_OR, HHA_Evects_OS, AA_Evects, AA_Evects_OR, AA_Evects_OS, HA_Evects, HA_Evects_OR, HA_Evects_OS;

	static Matrix projections_COV, normed_projections_COV, weighted_projections_COV, weighted_normed_projections_COV;

	static Matrix projections_dist_COV, normed_projections_dist_COV, projections_dist_CORR, normed_projections_dist_CORR, weighted_normed_projections_dist_COV,
			weighted_normed_projections_dist_CORR, weighted_projections_dist_CORR, weighted_projections_dist_COV, projections_dist_PCORR, normed_projections_dist_PCORR,
			weighted_projections_dist_PCORR, weighted_normed_projections_dist_PCORR;

	static Matrix Residue_Generalized_Coordinates_Global_AA_Cart, Residue_Generalized_Coordinates_Global_HA_Cart, Convoluted_Eigenvectors_AA_Cart, Convoluted_Square_Modes_AA_Cart,
			Convoluted_DVs_AA_Cart, Convoluted_Eigenvectors_HA_Cart, Convoluted_Square_Modes_HA_Cart, Convoluted_DVs_HA_Cart;

	static List<Atom> ref_atoms_AA, ref_atoms_BB, ref_atoms_CA, ref_atoms_HA;

	static List<Atom> ref_subset_atoms_AA, ref_subset_atoms_BB, ref_subset_atoms_CA, ref_subset_atoms_HA, ref_subset_atoms_AL, ref_subset_atoms_DP1, ref_subset_atoms_DP2,
			ref_subset_atoms_local, ref_subset_atoms_pairs, ref_subset_atoms_HAA, ref_subset_atoms_HHA, amended_subset_atoms_CA_BB_VIZ;

	static List<List<Atom>> Residue_Atoms_List_AA, Residue_Atoms_List_HA, Residue_Atoms_List_Local, Residue_Atoms_List_pairs;

	static List<Matrix> residue_coordinates_list, residue_reference_coordinates_list;

	static List<Atom_ID_Pair> atom_ID_pairs_read;

	static List<Residue_ID_Pair> residue_ID_pairs_read;

	static List<Residue_ID_Pair> residue_list_AA, residue_list_AA_orig, residue_list_BB, residue_list_BB_orig, residue_list_HA, residue_list_HA_orig, residue_list_CA,
			residue_list_CA_orig, residue_list_local, residue_list_local_orig, residue_list_pairs, residue_list_pairs_orig, residue_list_HAA, residue_list_HAA_orig,
			residue_list_HHA, residue_list_HHA_orig;

	static boolean doInterSSA_AA, doInterSSA_HA, verbose, doStatThresholds, outputAlignedCoords, doInterSSA_OR_OS, exist, success;
	static BufferedReader input_reader;
	static DateUtils now;
	static StringTokenizer sToken;
	static NumberFormat nf0, nf3, nf6, df;
	static RoundingMode rm;
	static long jediStart, jediEnd, startTime, endTime, totalTime;
	static Hashtable<String, String> parameters, parameters_read;
	static JEDi_Do_PDB_Processing refPDB;
	static JEDi_Driver_MT jedi;
	static Input_Parameters ip;
	static Output_Control oc;
	static ExecutorService service = Executors.newCachedThreadPool();
	@SuppressWarnings("rawtypes")
	static Future future0, future1, future2, future3, future4, future5, future6, future7, future8, future9;

	// ****************************************** CONSTRUCTOR *********************************************************************************** //

	public JEDi_Driver_MT()
	{
		rm = RoundingMode.HALF_UP;
		df = new DecimalFormat("0.###E0");
		df.setRoundingMode(rm);

		nf0 = NumberFormat.getInstance();
		nf0.setRoundingMode(rm);
		nf0.setMaximumFractionDigits(0);
		nf0.setMinimumFractionDigits(0);

		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);

		nf6 = NumberFormat.getInstance();
		nf6.setMaximumFractionDigits(6);
		nf6.setMinimumFractionDigits(6);
		nf6.setRoundingMode(rm);

		parameters = new Hashtable<>();
		lines = new ArrayList<String>();
	}

	// **************************************** PRIVATE METHODS ********************************************************************************* //

	private static void create_directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
	}

	private static void read_input_file()
	{
		number_Of_Input_Lines = 0;
		System.out.println("Reading file: " + input_file);
		try
			{
				while ((line = input_reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						if (!(line.startsWith("#") || line.startsWith("-")))
							{
								sToken = new StringTokenizer(line);
								key = sToken.nextToken(delim);
								value = sToken.nextToken(delim);
								parameters.put(key, value);
								number_Of_Input_Lines++;
							}
					}
				input_reader.close();
				System.out.println("\tThe total number of KEY=VALUE Pairs in the input file is: " + number_Of_Input_Lines + "\n");
			}
		catch (IOException e)
			{
				System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void process_Reference_PDB()
	{
		refPDB = new JEDi_Do_PDB_Processing();
		refPDB.parse_Reference_PDB();

		ref_atoms_AA = refPDB.get_ref_Atoms_AA();
		ref_atoms_BB = refPDB.get_ref_Atoms_BB();
		ref_atoms_CA = refPDB.get_ref_Atoms_CA();
		ref_atoms_HA = refPDB.get_ref_Atoms_HA();

		number_of_atoms_REF = refPDB.getNumber_of_atoms_AA();
		number_of_residues_REF = refPDB.getNumber_of_residues();
		number_of_heavy_atoms_REF = ref_atoms_HA.size();

		original_reference_PDB_coordinates_AA = refPDB.getReference_PDB_coordinates_AA();
		original_reference_PDB_coordinates_BB = refPDB.getReference_PDB_coordinates_BB();
		original_reference_PDB_coordinates_CA = refPDB.getReference_PDB_coordinates_CA();
		original_reference_PDB_coordinates_HA = refPDB.getReference_PDB_coordinates_HA();

		chain_ids_read = refPDB.get_Chain_IDs_Read();
		residue_ID_pairs_read = refPDB.getResidue_ID_pairs_read();
		atom_ID_pairs_read = refPDB.getAtom_ID_pairs_read();

		atoms_read = refPDB.get_Atoms_Read();
		heavy_atoms_read = refPDB.get_Heavy_Atoms_Read();
		alpha_carbons_read = refPDB.get_Alpha_Carbons_Read();

		residues_read = refPDB.get_Residue_Numbers_Read();

		numbers_Of_Atoms_in_Residues = refPDB.get_Number_of_Atoms_in_Residues();
		numbers_Of_Heavy_Atoms_in_Residues = refPDB.getNumber_of_Heavy_Atoms_in_Residues();

		ROWS_AA = original_reference_PDB_coordinates_AA.getRowDimension();
		ROWS_BB = original_reference_PDB_coordinates_BB.getRowDimension();
		ROWS_CA = original_reference_PDB_coordinates_CA.getRowDimension();
		ROWS_HA = original_reference_PDB_coordinates_HA.getRowDimension();

		name = "numbers_Of_Atoms_in_Residues.txt";
		path = Input_Parameters.OUT_DIR + name;
		ArrayList<String> NAR = new ArrayList<String>(numbers_Of_Atoms_in_Residues.size());
		for (int num = 0; num < numbers_Of_Atoms_in_Residues.size(); num++)
			{
				String cat = residue_ID_pairs_read.get(num).toShortString() + "\t" + numbers_Of_Atoms_in_Residues.get(num);
				NAR.add(cat);
			}

		if (Input_Parameters.verbose) List_IO.write_String_List(NAR, path);

		name = "numbers_Of_Heavy_Atoms_in_Residues.txt";
		path = Input_Parameters.OUT_DIR + name;
		ArrayList<String> NHR = new ArrayList<String>(numbers_Of_Heavy_Atoms_in_Residues.size());
		for (int num = 0; num < numbers_Of_Heavy_Atoms_in_Residues.size(); num++)
			{
				String cat = residue_ID_pairs_read.get(num).toShortString() + "\t" + numbers_Of_Heavy_Atoms_in_Residues.get(num);
				NHR.add(cat);
			}
		if (Input_Parameters.verbose) List_IO.write_String_List(NHR, path);

		name = "All_PDB_Residues_JEDi.txt";
		path = Input_Parameters.OUT_DIR + name;
		List_IO.write_ResidueID_Pair_List(residue_ID_pairs_read, path);

		oc.write_Reference_PDB_Log(number_of_residues_REF, number_of_atoms_REF, number_of_heavy_atoms_REF);
	}

	private static void read_PDB_Files()
	{
		if (Input_Parameters.doREAD_ARCHIVE)
			{
				if (Input_Parameters.ARCHIVE_NAME.endsWith("zip"))
					{
						original_PDB_coordinates_AA = refPDB.get_Matrix_of_Original_PDB_Coords_from_ZIP_Archive();
					}
				if (Input_Parameters.ARCHIVE_NAME.endsWith("tar.bz2") || Input_Parameters.ARCHIVE_NAME.endsWith("tar.gz"))
					{
						original_PDB_coordinates_AA = refPDB.get_Matrix_of_Original_PDB_Coords_from_TAR_Archive();
					}
			}
		else
			original_PDB_coordinates_AA = refPDB.get_Matrix_of_Original_PDB_Coords_from_PDB_Files();

		COLS = original_PDB_coordinates_AA.getColumnDimension();
		pdb_file_names = refPDB.get_PDB_file_names();
		number_conformations = COLS;

		name = "original_PDB_Coordinates_AA.txt.bz2";
		path = Input_Parameters.OUT_DIR + name;
		Matrix_IO.write_BZ2_Matrix(original_PDB_coordinates_AA, path, 9, 3);

		name = "PDB_Read_Log.txt";
		path = Input_Parameters.OUT_DIR + name;
		List_IO.write_String_List(pdb_file_names, path);
	}

	private static void read_All_Atom_Coordinate_File()
	{
		JEDi_Get_Coordinates_from_Matrix mat_coords = new JEDi_Get_Coordinates_from_Matrix(Input_Parameters.DIRECTORY, Input_Parameters.ORIGINAL_PDB_COORDS);
		original_PDB_coordinates_AA = mat_coords.get_Original_PDB_coordinates();
		number_conformations = original_PDB_coordinates_AA.getColumnDimension();

		oc.write_Coordinate_File_Log(ROWS_AA, number_conformations);
	}

	private static void do_Down_Sample()
	{
		int rows = original_PDB_coordinates_AA.getRowDimension();
		int cols = original_PDB_coordinates_AA.getColumnDimension();
		int size = (cols / Input_Parameters.DOWNSAMPLE);
		System.out.println("\tDown-Sampling the coordinates matrix:");
		System.out.println("\tThe number of frames in the original matrix is: " + cols);
		System.out.println("\tThe number of frames in the reduced matrix will be: " + size);
		int[] column_indices = new int[size];
		for (int i = 0; i < size; i++)
			{
				column_indices[i] = i * Input_Parameters.DOWNSAMPLE;
			}
		Matrix Reduced_Matrix = original_PDB_coordinates_AA.getMatrix(0, rows - 1, column_indices);
		original_PDB_coordinates_AA = Reduced_Matrix;

		oc.write_Down_Sample_Log(cols, size);
	}

	private static void do_Frame_Select()
	{
		int rows = original_PDB_coordinates_AA.getRowDimension();
		int cols = original_PDB_coordinates_AA.getColumnDimension();

		System.out.println("Selecting frames from the coordinates matrix:");
		System.out.println("\tFirst frame is: " + Input_Parameters.FRAME_START + " The last frame is: " + Input_Parameters.FRAME_END);
		System.out.println("\tThe number of frames before selection: " + (cols));
		System.out.println("\tThe number of frames in the final matrix will be: " + (Input_Parameters.FRAME_END - Input_Parameters.FRAME_START + 1));

		Matrix Reduced_Matrix = original_PDB_coordinates_AA.getMatrix(0, rows - 1, Input_Parameters.FRAME_START - 1, Input_Parameters.FRAME_END - 1);
		original_PDB_coordinates_AA = Reduced_Matrix;

		oc.write_Frame_Select_Log(cols);
	}

	private static void get_Subsets()
	{
		if (Input_Parameters.doAA)
			{
				get_All_Atom_Subset();
			}

		if (Input_Parameters.doBB)
			{
				get_BackBone_Subset();
			}

		if (Input_Parameters.doHA)
			{
				get_Heavy_Atom_Subset();
			}

		if (Input_Parameters.doCA)
			{
				get_Alpha_Carbon_Subset();
			}

		if (Input_Parameters.doRESIDUE_INDIVIDUAL)
			{
				get_Individual_Residue_Subset();
			}

		if (Input_Parameters.doRESIDUE_PAIRS)
			{
				get_Residue_Pairs_Subset();
			}
		if (Input_Parameters.doHIERARCHICAL_AA)
			{
				get_Hierarchical_All_Atom_Subset();
			}

		if (Input_Parameters.doHIERARCHICAL_HA)
			{
				get_Hierarchical_Heavy_Atom_Subset();
			}

		if (Input_Parameters.doATOM_LIST)
			{
				get_Atom_List_Subset();
			}

		if (Input_Parameters.doDIST_PAIRS)
			{
				get_Distance_Pairs_Subset();
			}
	}

	private static void get_Distance_Pairs_Subset()
	{
		System.out.println("\n\tGetting DISTANCE PAIR Subset: ");

		refPDB.read_Atom_Pairs(Input_Parameters.ATOM_PAIRS_LIST);

		atom_list_dp1 = refPDB.getAtom_list1();
		atom_list_dp2 = refPDB.getAtom_list2();
		atom_list_dp_original1 = refPDB.getAtom_list1_original();
		atom_list_dp_original2 = refPDB.getAtom_list2_original();
		chain_idents1 = refPDB.getChain_ids1();
		chain_idents2 = refPDB.getChain_ids2();
		number_of_atom_pairs_DP_SS = atom_list_dp1.size();
		ref_subset_atoms_DP1 = new ArrayList<Atom>();
		ref_subset_atoms_DP2 = new ArrayList<Atom>();

		HashMap<Integer, Atom> DP_Pairs = new HashMap<Integer, Atom>();
		for (Atom atm : ref_atoms_AA)
			{
				DP_Pairs.put(atm.atom_number, atm);
			}

		for (Integer x : atom_list_dp_original1)
			{
				if (DP_Pairs.containsKey(x))
					{
						Atom dp = new Atom(DP_Pairs.get(x));
						ref_subset_atoms_DP1.add(dp);
					}
			}
		for (Integer x : atom_list_dp_original2)
			{
				if (DP_Pairs.containsKey(x))
					{
						Atom dp = new Atom(DP_Pairs.get(x));
						ref_subset_atoms_DP2.add(dp);
					}
			}

		JEDi_Get_Distances_for_Atom_Pairs ref_dist_pairs = new JEDi_Get_Distances_for_Atom_Pairs(original_reference_PDB_coordinates_AA, atom_list_dp1, atom_list_dp2);
		reference_distances = ref_dist_pairs.Get_Ref_Distances();

		JEDi_Get_Distances_for_Atom_Pairs dist_pairs = new JEDi_Get_Distances_for_Atom_Pairs(original_PDB_coordinates_AA, atom_list_dp1, atom_list_dp2);
		distances = dist_pairs.get_Distances();


		// Mode check:
		if (Input_Parameters.MODES_DISTANCE_PAIRS > number_of_atom_pairs_DP_SS)
			{
				Input_Parameters.MODES_DISTANCE_PAIRS = number_of_atom_pairs_DP_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_DISTANCE_PAIRS reset to the number of available modes: " + number_of_atom_pairs_DP_SS);
			}

		oc.write_Subset_Log_DP(number_of_atom_pairs_DP_SS, ref_subset_atoms_DP1, ref_subset_atoms_DP2);
	}

	private static void get_Atom_List_Subset()
	{
		System.out.println("\n\tGetting ATOM LIST Subset: ");

		refPDB.read_Atoms_List(Input_Parameters.ATOMS_LIST);
		atom_list = refPDB.getAtom_list();
		atom_list_original = refPDB.getAtom_list_original();
		chain_idents0 = refPDB.getChain_ids0();
		number_of_atoms_AL_SS = atom_list.size();
		ref_subset_atoms_AL = new ArrayList<Atom>();

		HashMap<Integer, Atom> Atom_List = new HashMap<Integer, Atom>();
		for (Atom atm : ref_atoms_AA)
			{
				Atom_List.put(atm.atom_number, atm);
			}

		for (Integer x : atom_list_original)
			{
				if (Atom_List.containsKey(x))
					{
						Atom dp = new Atom(Atom_List.get(x));
						ref_subset_atoms_AL.add(dp);
					}
			}

		subset_reference_PDB_coordinates_AL = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_AL);
		subset_PDB_coordinates_AL = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, atom_list);

		// Mode check:
		if (Input_Parameters.MODES_ATOMS_LIST > (3 * number_of_atoms_AL_SS))
			{
				Input_Parameters.MODES_ATOMS_LIST = (3 * number_of_atoms_AL_SS);
				System.err.println("WARNING: NUMBER_OF_MODES_ATOMS_LIST reset to the number of available modes: " + (3 * number_of_atoms_AL_SS));
			}

		oc.write_Subset_Log(cAL, number_of_atoms_AL_SS, number_of_atoms_AL_SS, subset_reference_PDB_coordinates_AL, subset_PDB_coordinates_AL, ref_subset_atoms_AL);
	}

	private static void get_Hierarchical_Heavy_Atom_Subset()
	{
		System.out.println("\n\tGetting the HIERARCHICAL HEAVY ATOM Subset: \n");

		residue_list_HHA = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_HIERARCHICAL_HA);
		residue_list_HHA_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_HHA = refPDB.get_Reference_Atom_Subset_List(ref_atoms_HA, residue_list_HHA_orig);
		SS_Atom_Numbers_Hierarchical_HA = refPDB.get_SS_Atom_Numbers();

		Residue_Atoms_List_HA = refPDB.get_Residue_Atoms_List(ref_subset_atoms_HHA, residue_list_HHA_orig);

		number_of_atoms_hierarchical_HA_SS = ref_subset_atoms_HHA.size();
		number_of_residues_HHA_SS = residue_list_HHA_orig.size();

		subset_reference_PDB_coordinates_HHA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_HHA);
		subset_PDB_coordinates_hierarchical_HA = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_Hierarchical_HA);

		// Mode check:
		if (Input_Parameters.MODES_HIERARCHICAL_HA > Input_Parameters.MODES_EIGEN_RESIDUE_HA * number_of_residues_HHA_SS)
			{
				Input_Parameters.MODES_HIERARCHICAL_HA = Input_Parameters.MODES_EIGEN_RESIDUE_HA * number_of_residues_HHA_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_HIERARCHICAL_HA reset to the number of available modes: "
						+ (Input_Parameters.MODES_EIGEN_RESIDUE_HA * number_of_residues_HHA_SS));
			}

		oc.write_Subset_Log(HHA, number_of_residues_HHA_SS, number_of_atoms_hierarchical_HA_SS, subset_reference_PDB_coordinates_HHA, subset_PDB_coordinates_hierarchical_HA,
				ref_subset_atoms_HHA);
	}

	private static void get_Hierarchical_All_Atom_Subset()
	{
		System.out.println("\n\tGetting the HIERARCHICAL ALL ATOM Subset: \n");

		residue_list_HAA = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_HIERARCHICAL_AA);
		residue_list_HAA_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_HAA = refPDB.get_Reference_Atom_Subset_List(ref_atoms_AA, residue_list_HAA_orig);
		SS_Atom_Numbers_Hierarchical_AA = refPDB.get_SS_Atom_Numbers();

		number_of_atoms_hierarchical_AA_SS = ref_subset_atoms_HAA.size();
		number_of_residues_HAA_SS = residue_list_HAA_orig.size();

		subset_reference_PDB_coordinates_HAA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_HAA);
		subset_PDB_coordinates_hierarchical_AA = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_Hierarchical_AA);

		Residue_Atoms_List_AA = refPDB.get_Residue_Atoms_List(ref_subset_atoms_HAA, residue_list_HAA_orig);

		// Mode check:
		if (Input_Parameters.MODES_HIERARCHICAL_AA > Input_Parameters.MODES_EIGEN_RESIDUE_AA * number_of_residues_HAA_SS)
			{
				Input_Parameters.MODES_HIERARCHICAL_AA = Input_Parameters.MODES_EIGEN_RESIDUE_AA * number_of_residues_HAA_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_HIERARCHICAL_AA reset to the number of available modes: "
						+ (Input_Parameters.MODES_EIGEN_RESIDUE_AA * number_of_residues_HAA_SS));
			}

		oc.write_Subset_Log(HAA, number_of_residues_HAA_SS, number_of_atoms_hierarchical_AA_SS, subset_reference_PDB_coordinates_HAA, subset_PDB_coordinates_hierarchical_AA,
				ref_subset_atoms_HAA);
	}

	private static void get_Residue_Pairs_Subset()
	{
		System.out.println("\n\tGetting RESIDUE PAIRS Subset: \n");

		List<Integer> num_atms_in_residue_pairs_subset = new ArrayList<Integer>();
		residue_list_pairs_orig = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_PAIRS);

		ref_subset_atoms_pairs = refPDB.get_Reference_Atom_Subset_List(ref_atoms_AA, residue_list_pairs_orig);
		SS_Atom_Numbers_pairs = refPDB.get_SS_Atom_Numbers();

		Residue_Atoms_List_pairs = refPDB.get_Residue_Atoms_List(ref_subset_atoms_pairs, residue_list_pairs_orig);

		for (Residue_ID_Pair pair : residue_list_pairs_orig)
			{
				int numAtomsAA = 0;
				String id = pair.getChain_ID();
				int num = pair.getResidue_Number();
				for (Atom a : ref_subset_atoms_pairs)
					{
						if (a.chainID.equals(id) && a.res_number == num)
							{
								numAtomsAA++;
							}
					}
				num_atms_in_residue_pairs_subset.add(numAtomsAA);
			}

		number_of_atoms_residue_pairs_SS = ref_subset_atoms_pairs.size();
		number_of_residues_pairs_SS = ref_subset_atoms_pairs.size();

		subset_reference_PDB_coordinates_pairs = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_pairs);
		subset_PDB_coordinates_pairs = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_pairs);

		residue_reference_coordinates_list = refPDB.get_Residue_Coords(subset_reference_PDB_coordinates_pairs, num_atms_in_residue_pairs_subset);
		residue_coordinates_list = refPDB.get_Residue_Coords(subset_PDB_coordinates_pairs, num_atms_in_residue_pairs_subset);

		oc.write_Subset_Log(Res_PAIRS, number_of_atoms_residue_pairs_SS, number_of_atoms_residue_pairs_SS, subset_reference_PDB_coordinates_pairs, subset_PDB_coordinates_pairs,
				ref_subset_atoms_pairs);
	}

	private static void get_Individual_Residue_Subset()
	{
		System.out.println("\n\tGetting INDIVIDUAL RESIDUE Subset: \n");

		residue_list_local = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_INDIVIDUAL);
		residue_list_local_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_local = refPDB.get_Reference_Atom_Subset_List(ref_atoms_AA, residue_list_local_orig);
		SS_Atom_Numbers_local = refPDB.get_SS_Atom_Numbers();

		Residue_Atoms_List_Local = refPDB.get_Residue_Atoms_List(ref_subset_atoms_local, residue_list_local_orig);

		number_of_atoms_local_SS = ref_subset_atoms_local.size();
		number_of_residues_local_SS = residue_list_local.size();

		subset_reference_PDB_coordinates_local = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_local);
		subset_PDB_coordinates_local = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_local);

		oc.write_Subset_Log(Res_AA, number_of_residues_local_SS, number_of_atoms_local_SS, subset_reference_PDB_coordinates_local, subset_PDB_coordinates_local,
				ref_subset_atoms_local);
	}

	private static void get_Alpha_Carbon_Subset()
	{
		System.out.println("\n\tGetting the ALPHA CARBON Subset: \n");

		residue_list_CA = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_ALPHA_CARBON);
		residue_list_CA_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_CA = refPDB.get_Reference_Atom_Subset_List(ref_atoms_CA, residue_list_CA_orig);
		SS_Atom_Numbers_CA = refPDB.get_SS_Atom_Numbers();

		if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
			{
				System.out.println("\n\tGetting Alpha-Carbon Amended Subset for VIZ: \n");
				amended_subset_atoms_CA_BB_VIZ = refPDB.get_Reference_Atom_Subset_List(ref_atoms_BB, residue_list_CA_orig);
			}
		number_of_residues_CA_SS = ref_subset_atoms_CA.size();

		subset_reference_PDB_coordinates_CA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_CA);
		subset_PDB_coordinates_CA = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_CA);

		// Mode check:
		if (Input_Parameters.MODES_ALPHA_CARBON > 3 * number_of_residues_CA_SS)
			{
				Input_Parameters.MODES_ALPHA_CARBON = 3 * number_of_residues_CA_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_ALPHA_CARBON reset to the number of available modes: " + (3 * number_of_residues_CA_SS));
			}

		oc.write_Subset_Log(cAC, number_of_residues_CA_SS, number_of_residues_CA_SS, subset_reference_PDB_coordinates_CA, subset_PDB_coordinates_CA, ref_subset_atoms_CA);
	}

	private static void get_Heavy_Atom_Subset()
	{
		System.out.println("\n\tGetting the HEAVY ATOM Subset: \n");

		residue_list_HA = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_HEAVY_ATOM);
		residue_list_HA_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_HA = refPDB.get_Reference_Atom_Subset_List(ref_atoms_HA, residue_list_HA_orig);
		SS_Atom_Numbers_HA = refPDB.get_SS_Atom_Numbers();

		number_of_atoms_HA_SS = ref_subset_atoms_HA.size();
		number_of_residues_HA_SS = residue_list_HA.size();

		subset_reference_PDB_coordinates_HA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_HA);
		subset_PDB_coordinates_HA = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_HA);

		// Mode check:
		if (Input_Parameters.MODES_HEAVY_ATOM > (3 * number_of_atoms_HA_SS))
			{
				Input_Parameters.MODES_HEAVY_ATOM = 3 * number_of_atoms_HA_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_HEAVY_ATOM reset to the number of available modes: " + (3 * number_of_atoms_HA_SS));
			}

		oc.write_Subset_Log(cHA, number_of_residues_HA_SS, number_of_atoms_HA_SS, subset_reference_PDB_coordinates_HA, subset_PDB_coordinates_HA, ref_subset_atoms_HA);
	}

	private static void get_BackBone_Subset()
	{
		System.out.println("\n\tGetting the BACKBONE Subset: \n");

		residue_list_BB = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_BACKBONE);
		residue_list_BB_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_BB = refPDB.get_Reference_Atom_Subset_List(ref_atoms_BB, residue_list_BB_orig);
		SS_Atom_Numbers_BB = refPDB.get_SS_Atom_Numbers();

		number_of_atoms_BB_SS = ref_subset_atoms_BB.size();
		number_of_residues_BB_SS = residue_list_BB.size();

		subset_reference_PDB_coordinates_BB = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_BB);
		subset_PDB_coordinates_BB = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_BB);

		// Mode check:
		if (Input_Parameters.MODES_BACKBONE > 3 * number_of_atoms_BB_SS)
			{
				Input_Parameters.MODES_BACKBONE = 3 * number_of_atoms_BB_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_BACKBONE reset to the number of available modes: " + (3 * number_of_atoms_BB_SS));
			}
		oc.write_Subset_Log(cBB, number_of_residues_BB_SS, number_of_atoms_BB_SS, subset_reference_PDB_coordinates_BB, subset_PDB_coordinates_BB, ref_subset_atoms_BB);
	}

	private static void get_All_Atom_Subset()
	{
		System.out.println("\n\tGetting the ALL ATOM Subset: \n");

		residue_list_AA = refPDB.read_Residue_List(Input_Parameters.RESIDUE_LIST_ALL_ATOM);
		residue_list_AA_orig = refPDB.get_Residue_List_Original();

		ref_subset_atoms_AA = refPDB.get_Reference_Atom_Subset_List(ref_atoms_AA, residue_list_AA_orig);
		SS_Atom_Numbers_AA = refPDB.get_SS_Atom_Numbers();

		number_of_atoms_AA_SS = ref_subset_atoms_AA.size();
		number_of_residues_AA_SS = residue_list_AA.size();

		subset_reference_PDB_coordinates_AA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_AA);
		subset_PDB_coordinates_AA = refPDB.get_Subset_Coords(original_PDB_coordinates_AA, SS_Atom_Numbers_AA);

		// Mode check:
		if (Input_Parameters.MODES_ALL_ATOM > 3 * number_of_atoms_AA_SS)
			{
				Input_Parameters.MODES_ALL_ATOM = 3 * number_of_atoms_AA_SS;
				System.err.println("WARNING: NUMBER_OF_MODES_ALL_ATOM reset to the number of available modes: " + (3 * number_of_atoms_AA_SS));
			}

		oc.write_Subset_Log(cAA, number_of_residues_AA_SS, number_of_atoms_AA_SS, subset_reference_PDB_coordinates_AA, subset_PDB_coordinates_AA, ref_subset_atoms_AA);
	}

	private static void align_Subsets()
	{
		verbose = Input_Parameters.verbose;
		doStatThresholds = Input_Parameters.do_StatThresholds;
		outputAlignedCoords = Input_Parameters.doOutputCoordinates;

		if (Input_Parameters.doAA)
			{
				align_All_Atom_Subset();
			}

		if (Input_Parameters.doBB)
			{
				align_BackBone_Atom_Subset();
			}

		if (Input_Parameters.doHA)
			{
				align_Heavy_Atom_Subset();
			}

		if (Input_Parameters.doCA)
			{
				align_Alpha_Carbon_Subset();
			}

		if (Input_Parameters.do_atom_list)
			{
				align_Atom_List_Subset();
			}

		if (Input_Parameters.do_residue_individual)
			{
				align_All_Atom_Residue_Subset();
			}

		if (Input_Parameters.doHIERARCHICAL_AA)
			{
				align_Hierarchical_All_Atom_Subset();
			}

		if (Input_Parameters.doHIERARCHICAL_HA)
			{
				align_Hierarchical_Heavy_Atom_Subset();
			}

		if (Input_Parameters.doDIST_PAIRS)
			{
				process_Distance_Pairs_Subset();
			}
	}

	private static void align_Hierarchical_All_Atom_Subset()
	{
		System.out.println("\n\tProcessing the HIERARCHICAL ALL ATOM Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_hierarchical_AA, subset_reference_PDB_coordinates_HAA);

		aligned_subset_PDB_coordinates_hierarchical_AA = tf_coords.get_SS_Transformed_coords();
		aligned_subset_REF_PDB_coordinates_hierarchical_AA = tf_coords.get_Aligned_Reference_Coordinates();

		String out = Input_Parameters.OUT_DIR + HAA + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_HAA_SS + "_" + number_of_atoms_hierarchical_AA_SS + "_Hierarchical_All_Atom_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_hierarchical_AA, path, 9, 3);
			}

		aligned_atomic_RMSFs_Hierarchical_AA = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_HAA, aligned_atomic_RMSFs_Hierarchical_AA);
		name = "ss_" + number_of_residues_HAA_SS + "_" + number_of_atoms_hierarchical_AA_SS + "_Hierarchical_AA_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_Hierarchical_AA, ref_subset_atoms_HAA);
		RMSD_Plot.create_RMSF_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_Hierarchical_AA);
		RMSD_Plot.create_RMSD_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "HIERARCHICAL_ALL_ATOM_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "HIERARCHICAL_ALL_ATOM_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_Hierarchical_AA, path, 12);
				name = "HIERARCHICAL_ALL_ATOM_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_HAA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_HAA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_HAA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_HAA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_HAA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_HAA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_HAA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_HAA, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(HAA, ref_subset_atoms_HAA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(HAA);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_hierarchical_Outliers_REMOVED_AA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_hierarchical_Outliers_SELECTED_AA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HAA = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_HAA = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HAA, path, 12, 0);

								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HAA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HAA, ref_subset_atoms_HAA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HAA = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_HAA = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HAA, path, 12, 0);

								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HAA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HAA, ref_subset_atoms_HAA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HAA, ref_subset_atoms_HAA, 1);
					}
			}
	}

	private static void align_Atom_List_Subset()
	{
		System.out.println("\n\tProcessing the ATOM LIST Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_AL, subset_reference_PDB_coordinates_AL);

		aligned_subset_REF_PDB_coordinates_AL = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_AL = tf_coords.get_SS_Transformed_coords();

		String out = Input_Parameters.OUT_DIR + cAL + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_atoms_AL_SS + "_" + number_of_atoms_AL_SS + "_Atom_List_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_AL, path, 9, 3);
			}

		aligned_atomic_RMSFs_AL = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_AL, aligned_atomic_RMSFs_AL);
		name = "ss_" + number_of_atoms_AL_SS + "_AL_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "ATOM_LIST_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_AL, ref_subset_atoms_AL);
		RMSD_Plot.create_RMSF_XY_Chart(out, "ATOM_LIST_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_AL);
		RMSD_Plot.create_RMSD_XY_Chart(out, "ATOM_LIST_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "ATOM_LIST_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "ATOM_LIST_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_AL, path, 12);
				name = "ATOM_LIST_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_AL, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_AL, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_AL, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_AL, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_AL, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_AL, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_AL, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_AL, 4);

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_Outliers_REMOVED_AL = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_Outliers_SELECTED_AL = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_AL = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_AL = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_AL, path, 12, 0);

								name = "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_AL, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_AL, ref_subset_atoms_AL, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_AL = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_AL = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_AL, path, 12, 0);

								name = "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_AL, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_AL, ref_subset_atoms_AL, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ATOM_LIST_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_AL, ref_subset_atoms_AL, 1);
					}
			}
	}

	private static void align_All_Atom_Residue_Subset()
	{
		System.out.println("\n\tProcessing the RESIDUE AA Subset: \n");

		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_local, subset_reference_PDB_coordinates_local);

		aligned_subset_REF_PDB_coordinates_local = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_local = tf_coords.get_SS_Transformed_coords();

		String out = Input_Parameters.OUT_DIR + Res_AA + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_local_SS + "_" + number_of_atoms_local_SS + "_Individual_Residues_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_local, path, 9, 3);
			}

		aligned_atomic_RMSFs_local = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_local, aligned_atomic_RMSFs_local);
		name = "ss_" + number_of_atoms_local_SS + "_Individual_Residue_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "AA_RESIDUE_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_local, ref_subset_atoms_local);
		RMSD_Plot.create_RMSF_XY_Chart(out, "AA_RESIDUE_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_local);
		RMSD_Plot.create_RMSD_XY_Chart(out, "AA_RESIDUE_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "AA_RESIDUE_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "AA_RESIDUE_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_local, path, 12);
				name = "AA_RESIDUE_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_local, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_local, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_local, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_local, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_local, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_local, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_local, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_local, 4);

		if (Input_Parameters.doOutlierProcessing)
			{
				aligned_subset_PDB_coordinates_Outliers_REMOVED_local = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				aligned_subset_PDB_coordinates_Outliers_SELECTED_local = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						Matrix adjustments_per_variable_REMOVE_local = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						Matrix adjustments_per_variable_SELECT_local = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_local, path, 12, 0);

								name = "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_local, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_local, ref_subset_atoms_local, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						Matrix adjustments_per_variable_REMOVE_local = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						Matrix adjustments_per_variable_SELECT_local = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_local, path, 12, 0);

								name = "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_local, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_local, ref_subset_atoms_local, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "AA_RESIDUE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_local, ref_subset_atoms_local, 1);
					}
			}
	}

	private static void align_Alpha_Carbon_Subset()
	{
		System.out.println("\n\tProcessing the ALPHA CARBON Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_CA, subset_reference_PDB_coordinates_CA);

		aligned_subset_REF_PDB_coordinates_CA = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_CA = tf_coords.get_SS_Transformed_coords();

		String out = Input_Parameters.OUT_DIR + cAC + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_CA_SS + "_" + number_of_residues_CA_SS + "_Alpha_Carbon_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_CA, path, 9, 3);
			}

		aligned_atomic_RMSFs_CA = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_CA, aligned_atomic_RMSFs_CA);
		name = "ss_" + number_of_residues_CA_SS + "_" + number_of_residues_CA_SS + "_CA_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "ALPHA_CARBON_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_CA, ref_subset_atoms_CA);
		RMSD_Plot.create_RMSF_XY_Chart(out, "ALPHA_CARBON_Subset_Atomic_RMSFs_XY", "Atom Number", aligned_atomic_RMSFs_CA);
		RMSD_Plot.create_RMSD_XY_Chart(out, "ALPHA_CARBON_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "ALPHA_CARBON_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "ALPHA_CARBON_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_CA, path, 12);
				name = "ALPHA_CARBON_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_CA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_CA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_CA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_CA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_CA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_CA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_CA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_CA, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(CA, ref_subset_atoms_CA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(CA);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_Outliers_REMOVED_CA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_Outliers_SELECTED_CA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_CA = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_CA = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_CA, path, 12, 0);

								name = "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_CA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_CA, ref_subset_atoms_CA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_CA = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_CA = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_CA, path, 12, 0);

								name = "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_CA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_CA, ref_subset_atoms_CA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALPHA_CARBON_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_CA, ref_subset_atoms_CA, 1);
					}
			}
	}

	private static void align_Heavy_Atom_Subset()
	{
		System.out.println("\n\tProcessing the HEAVY ATOM Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_HA, subset_reference_PDB_coordinates_HA);

		aligned_subset_REF_PDB_coordinates_HA = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_HA = tf_coords.get_SS_Transformed_coords();

		String out = Input_Parameters.OUT_DIR + cHA + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_HA_SS + "_" + number_of_atoms_HA_SS + "_Heavy_Atom_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_HA, path, 9, 3);
			}

		aligned_atomic_RMSFs_HA = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_HA, aligned_atomic_RMSFs_HA);
		name = "ss_" + number_of_residues_HA_SS + "_" + number_of_atoms_HA_SS + "_HA_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "HEAVY_ATOM_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_HA, ref_subset_atoms_HA);
		RMSD_Plot.create_RMSF_XY_Chart(out, "HEAVY_ATOM_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_HA);
		RMSD_Plot.create_RMSD_XY_Chart(out, "HEAVY_ATOM_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "HEAVY_ATOM_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "HEAVY_ATOM_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_HA, path, 12);
				name = "HEAVY_ATOM_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_HA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_HA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_HA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_HA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_HA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_HA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_HA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_HA, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(HA, ref_subset_atoms_HA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(HA);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_Outliers_REMOVED_HA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_Outliers_SELECTED_HA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HA = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_HA = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HA, path, 12, 0);

								name = "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HA, ref_subset_atoms_HA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HA = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_HA = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HA, path, 12, 0);

								name = "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HA, ref_subset_atoms_HA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HA, ref_subset_atoms_HA, 1);
					}
			}
	}

	private static void align_BackBone_Atom_Subset()
	{
		System.out.println("\n\tProcessing the BACKBONE Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_BB, subset_reference_PDB_coordinates_BB);

		String out = Input_Parameters.OUT_DIR + cBB + File.separatorChar;
		create_directory(out);

		aligned_subset_REF_PDB_coordinates_BB = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_BB = tf_coords.get_SS_Transformed_coords();

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_BB_SS + "_" + number_of_atoms_BB_SS + "_BackBone_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_BB, path, 9, 3);
			}

		aligned_atomic_RMSFs_BB = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_BB, aligned_atomic_RMSFs_BB);
		name = "ss_" + number_of_residues_BB_SS + "_" + number_of_atoms_BB_SS + "_BB_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "BACKBONE_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_BB, ref_subset_atoms_BB);
		RMSD_Plot.create_RMSF_XY_Chart(out, "BACKBONE_Subset_Atomic_RMSFs_XY", "Atom Number", aligned_atomic_RMSFs_BB);
		RMSD_Plot.create_RMSD_XY_Chart(out, "BACKBONE_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "BACKBONE_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "BACKBONE_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_BB, path, 12);
				name = "BACKBONE_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_BB, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_BB, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_BB, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_BB, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_BB, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_BB, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_BB, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_BB, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(BB, ref_subset_atoms_BB, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(BB);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_Outliers_REMOVED_BB = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_Outliers_SELECTED_BB = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_BB = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_BB = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_AA, path, 12, 0);

								name = "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_AA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_BB = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_BB = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_BB, path, 12, 0);

								name = "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_BB, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_BB, ref_subset_atoms_BB, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_BB, ref_subset_atoms_BB, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_BB, ref_subset_atoms_BB, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "BACKBONE_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_BB, ref_subset_atoms_BB, 1);
					}
			}
	}

	private static void align_Hierarchical_Heavy_Atom_Subset()
	{
		System.out.println("\n\tProcessing the HIERARCHICAL HEAVY ATOM Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_hierarchical_HA, subset_reference_PDB_coordinates_HHA);

		aligned_subset_PDB_coordinates_hierarchical_HA = tf_coords.get_SS_Transformed_coords();
		aligned_subset_REF_PDB_coordinates_hierarchical_HA = tf_coords.get_Aligned_Reference_Coordinates();

		String out = Input_Parameters.OUT_DIR + HHA + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_HHA_SS + "_" + number_of_atoms_hierarchical_HA_SS + "_Hierarchical_Heavy_Atom_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_hierarchical_HA, path, 9, 3);
			}

		aligned_atomic_RMSFs_Hierarchical_HA = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_HHA, aligned_atomic_RMSFs_Hierarchical_HA);
		name = "ss_" + number_of_residues_HHA_SS + "_" + number_of_atoms_hierarchical_HA_SS + "_Hierarchical_HA_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_Hierarchical_HA, ref_subset_atoms_HHA);
		RMSD_Plot.create_RMSF_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_Hierarchical_HA);
		RMSD_Plot.create_RMSD_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "HIERARCHICAL_HEAVY_ATOM_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "HIERARCHICAL_HEAVY_ATOM_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_Hierarchical_HA, path, 12);
				name = "HIERARCHICAL_HEAVY_ATOM_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_HHA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_HHA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_HHA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_HHA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_HHA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_HHA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_HHA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_HHA, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(HHA, ref_subset_atoms_HHA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(HHA);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_hierarchical_Outliers_REMOVED_HA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_hierarchical_Outliers_SELECTED_HA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HHA = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_HHA = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HHA, path, 12, 0);

								name = "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HHA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HHA, ref_subset_atoms_HHA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_HHA = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_HHA = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_HHA, path, 12, 0);

								name = "HIERARCHICAL_ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_HHA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_HHA, ref_subset_atoms_HHA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "HIERARCHICAL_HEAVY_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_HHA, ref_subset_atoms_HHA, 1);
					}
			}
	}

	private static void align_All_Atom_Subset()
	{
		System.out.println("\n\tProcessing the ALL ATOM Subset: \n");
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(subset_PDB_coordinates_AA, subset_reference_PDB_coordinates_AA);

		aligned_subset_REF_PDB_coordinates_AA = tf_coords.get_Aligned_Reference_Coordinates();
		aligned_subset_PDB_coordinates_AA = tf_coords.get_SS_Transformed_coords();

		String out = Input_Parameters.OUT_DIR + cAA + File.separatorChar;
		create_directory(out);

		if (outputAlignedCoords)
			{
				name = "ss_" + number_of_residues_AA_SS + "_" + number_of_atoms_AA_SS + "_All_Atom_Aligned_Coordinates.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(aligned_subset_PDB_coordinates_AA, path, 9, 3);
			}

		aligned_atomic_RMSFs_AA = tf_coords.get_SS_RMSF();
		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_subset_atoms_AA, aligned_atomic_RMSFs_AA);

		name = "ss_" + number_of_residues_AA_SS + "_" + number_of_atoms_AA_SS + "_AA_RMSF_edited.pdb";
		path = out + name;
		PDB_IO.Write_PDB(path, edited_atoms);

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		RMSD_Plot.create_RMSF_Line_Chart(out, "ALL_ATOM_Subset_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs_AA, ref_subset_atoms_AA);
		RMSD_Plot.create_RMSF_XY_Chart(out, "ALL_ATOM_Subset_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs_AA);
		RMSD_Plot.create_RMSD_XY_Chart(out, "ALL_ATOM_Subset_Conformation_RMSDs", "Conformation Number", conf_rmsds);

		Matrix stats = tf_coords.get_SS_coordinate_STATS();

		if (verbose)
			{
				name = "ALL_ATOM_Subset_Coordinate_Stats.txt";
				path = out + name;
				Matrix_IO.write_Matrix(stats, path, 12, 6);
				name = "ALL_ATOM_Subset_Atomic_RMSFs.txt";
				path = out + name;
				List_IO.write_Double_List(aligned_atomic_RMSFs_AA, path, 12);
				name = "ALL_ATOM_Subset_Conformation_RMSDs.txt";
				path = out + name;
				List_IO.write_Double_List(conf_rmsds, path, 12);
			}

		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Variable_Means", "Atom Number", "Mean", stats, ref_subset_atoms_AA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Variable_Variances", "Atom Number", "Variance", stats, ref_subset_atoms_AA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Variable_Skews", "Atom Number", "Skew", stats, ref_subset_atoms_AA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_subset_atoms_AA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Variable_Means_XY", "Atom Index", "Mean", stats, ref_subset_atoms_AA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_subset_atoms_AA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_subset_atoms_AA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_subset_atoms_AA, 4);

		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(AA, ref_subset_atoms_AA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(AA);
			}

		if (Input_Parameters.doOutlierProcessing)
			{
				subset_PDB_coordinates_Outliers_REMOVED_AA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
				subset_PDB_coordinates_Outliers_SELECTED_AA = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_AA = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
						adjustments_per_variable_SELECT_AA = tf_coords.getAdjustments_per_variable_SELECT_MAD();
						if (verbose)
							{
								name = "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_AA, path, 12, 0);

								name = "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_AA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_MAD_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);
					}

				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						adjustments_per_variable_REMOVE_AA = tf_coords.getAdjustments_per_variable_REMOVE_Z();
						adjustments_per_variable_SELECT_AA = tf_coords.getAdjustments_per_variable_SELECT_Z();
						if (verbose)
							{
								name = "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_AA, path, 12, 0);

								name = "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z.txt";
								path = out + name;
								Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_AA, path, 12, 0);
							}
						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_REMOVE_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_REMOVE_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_Line_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z", "Atom Number", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);

						STATS_Plot.create_Variables_Stat_XY_Chart(out, "ALL_ATOM_Subset_Adjustments_per_Variable_SELECT_Outliers_Z_XY", "Atom Index", "Counts",
								adjustments_per_variable_SELECT_AA, ref_subset_atoms_AA, 1);
					}
			}
	}

	private static void process_Distance_Pairs_Subset()
	{
		Outlier_Processing op = new Outlier_Processing(distances);
		Matrix stats = op.getCoordinateStats();
		int number_of_pairs = atom_list_dp1.size();

		String out = Input_Parameters.OUT_DIR + DP + File.separatorChar;
		create_directory(out);

		if (Input_Parameters.verbose)
			{
				String name = "Distance_Pair_4_Moment_Stats.txt.bz2";
				path = out + name;
				Matrix_IO.write_BZ2_Matrix(stats, path, 12, 6);
			}

		MODES_Plot.createLineChart1SeriesDist(out, "Distance_Pair_Means", "Atom Index", "Mean", stats, atom_list_dp_original1, atom_list_dp_original2, 1);
		MODES_Plot.createLineChart1SeriesDist(out, "Distance_Pair_Variances", "Atom Index", "Variance", stats, atom_list_dp_original1, atom_list_dp_original2, 2);
		MODES_Plot.createLineChart1SeriesDist(out, "Distance_Pair_Skews", "Atom Index", "Skew", stats, atom_list_dp_original1, atom_list_dp_original2, 3);
		MODES_Plot.createLineChart1SeriesDist(out, "Distance_Pair_Kurtosis", "Atom Index", "Kurtosis", stats, atom_list_dp_original1, atom_list_dp_original2, 4);

		if (Input_Parameters.doOutlierProcessing)
			{
				if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						op.remove_outliers_MAD();
						distances_Outliers_REMOVED = op.remove_outliers_MAD();
						Matrix counts = op.getCounts_remove_MAD();

						if (Input_Parameters.verbose) Matrix_IO.write_BZ2_Matrix(counts, out + "ss_" + number_of_pairs + "_Adjusted_Outliers_Per_Variable_MAD.txt.bz2", 6, 0);
						MODES_Plot.createLineChart1SeriesDist(out, "Adjustments to Outliers", "Atom Pair Index", "Counts", counts, atom_list_dp_original1, atom_list_dp_original2,
								1);

						op.select_outliers_MAD();
						distances_Outliers_SELECTED = op.select_outliers_MAD();
						counts = op.getCounts_select_MAD();

						if (Input_Parameters.verbose) Matrix_IO.write_BZ2_Matrix(counts, out + "ss_" + number_of_pairs + "_Adjusted_Inliers_Per_Variable_MAD.txt.bz2", 6, 0);
						Plot_Line_Chart.createChart_One_Series(out, "Adjustments to Inliers", "Atom Pair Index", "Counts", counts);

					}
				if (Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						op.remove_outliers_Z();
						distances_Outliers_REMOVED = op.remove_outliers_Z();
						Matrix counts = op.getCounts_remove_Z();

						if (Input_Parameters.verbose) Matrix_IO.write_BZ2_Matrix(counts, out + "ss_" + number_of_pairs + "_Adjusted_Outliers_Per_Variable_Z.txt.bz2", 6, 0);
						MODES_Plot.createLineChart1SeriesDist(out, "Adjustments to Outliers", "Atom Pair Index", "Counts", counts, atom_list_dp_original1, atom_list_dp_original2,
								1);

						op.select_outliers_Z();
						distances_Outliers_SELECTED = op.select_outliers_Z();
						counts = op.getCounts_select_Z();

						if (Input_Parameters.verbose) Matrix_IO.write_BZ2_Matrix(counts, out + "ss_" + number_of_pairs + "_Adjusted_Inliers_Per_Variable_Z.txt.bz2", 6, 0);
						Plot_Line_Chart.createChart_One_Series(out, "Adjustments to Inliers", "Atom Pair Index", "Counts", counts);
					}
			}
		atomic_distance_means = op.get_means();
		atomic_distance_std_devs = op.get_std_deviations();
		path = out + "ss_" + number_of_pairs + "_Distance_Pair_Means_and_Std_Devs.txt";
		if (Input_Parameters.verbose)
			{
				Write_DP_Stats.write_stats(path, atomic_distance_means, atomic_distance_std_devs, chain_idents1, chain_idents2, atom_list_dp_original1, atom_list_dp_original2);
			}
	}

	private static void do_All_Atom_Residue_Pairs_Analysis()
	{
		System.out.println("\n\tProcessing the RESIDUE PAIRS Subset: \n");
		JEDi_Do_Residue resPairAA = new JEDi_Do_Residue(residue_list_pairs_orig, subset_PDB_coordinates_pairs, subset_reference_PDB_coordinates_pairs, Residue_Atoms_List_pairs);

		resPairAA.setDoFES(false);
		resPairAA.setDoKPCA(false);
		resPairAA.setDo_CORR(false);
		resPairAA.setDo_PCORR(false);
		resPairAA.setDoReduce(false);
		resPairAA.setDo_ESSENTIAL_VIZ(false);
		resPairAA.setDo_INDIVIDUAL_VIZ(false);
		resPairAA.set_out_dir_pca(Input_Parameters.OUT_DIR + Res_PAIRS + File.separatorChar);
		resPairAA.do_Residue_Pairs(residue_reference_coordinates_list, residue_coordinates_list);
		oc.write_Residue_Pair_Analysis_Log();
	}

	private static void do_Individual_All_Atom_Residue_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, String outDirVIZ, Matrix coords)
	{
		JEDi_Do_Residue LCresAA = new JEDi_Do_Residue(residue_list_local_orig, subset_PDB_coordinates_local, subset_reference_PDB_coordinates_local, coords,
				Residue_Atoms_List_Local);

		LCresAA.set_out_dir_pca(outDirPCA);
		LCresAA.set_out_dir_ssa(outDirSSA);
		if (Input_Parameters.doFES) LCresAA.set_out_dir_fes(outDirFES);
		if (Input_Parameters.doKPCA) LCresAA.set_out_dir_kpca(outDirKPCA);
		if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz) LCresAA.set_out_dir_viz(outDirVIZ);
		LCresAA.do_Individual_Residue_Analysis();
		oc.write_Individual_Residue_PCA_Log();
	}

	private static void do_Global_All_Atom_Residue_PCA(String outputDirectory, Matrix coords)
	{
		JEDi_Do_Residue GCresAA = new JEDi_Do_Residue(residue_list_HAA_orig, Input_Parameters.MODES_EIGEN_RESIDUE_AA, aligned_subset_PDB_coordinates_hierarchical_AA,
				aligned_subset_REF_PDB_coordinates_hierarchical_AA, coords, Residue_Atoms_List_AA, 0);

		GCresAA.setDoAllAtoms(true);
		GCresAA.setDoHeavyAtoms(false);
		GCresAA.do_Global_Residue_Analysis();

		Residue_Generalized_Coordinates_Global_AA_Cart = GCresAA.get_Residue_Generalized_Coordinates_Global();
		Residue_Eigenvectors_Global_AA = GCresAA.get_Residue_Eigenvectors();
		Residue_Delta_Vectors_Global_AA = GCresAA.getResidues_Centered_Data(); // For PCs, we need the mean centered data as the DVs...
		Residue_Centered_Coordinates_Global_AA = GCresAA.getResidues_Centered_Data();
	}

	private static void do_Global_Heavy_Atom_Residue_PCA(String outputDirectory, Matrix coords)
	{
		JEDi_Do_Residue GCresHA = new JEDi_Do_Residue(residue_list_HHA_orig, Input_Parameters.MODES_EIGEN_RESIDUE_HA, aligned_subset_PDB_coordinates_hierarchical_HA,
				aligned_subset_REF_PDB_coordinates_hierarchical_HA, coords, Residue_Atoms_List_HA, 0);

		GCresHA.setDoAllAtoms(false);
		GCresHA.setDoHeavyAtoms(true);
		GCresHA.do_Global_Residue_Analysis();

		Residue_Generalized_Coordinates_Global_HA_Cart = GCresHA.get_Residue_Generalized_Coordinates_Global();
		Residue_Eigenvectors_Global_HA = GCresHA.get_Residue_Eigenvectors();
		Residue_Delta_Vectors_Global_HA = GCresHA.getResidues_Centered_Data(); // For PCs, we need the mean centered data as the DVs...
		Residue_Centered_Coordinates_Global_HA = GCresHA.getResidues_Centered_Data();
	}

	private static void do_Hierarchical_All_Atom_PCA(String outDirPCA, String outDirFES, String outDirKPCA)
	{
		JEDi_Do_Hierarchical hss_AA_cart = new JEDi_Do_Hierarchical(Input_Parameters.MODES_EIGEN_RESIDUE_AA, Input_Parameters.MODES_HIERARCHICAL_AA,
				Residue_Generalized_Coordinates_Global_AA_Cart, Residue_Eigenvectors_Global_AA, Residue_Delta_Vectors_Global_AA);

		hss_AA_cart.setNumber_of_Atoms(number_of_atoms_hierarchical_AA_SS);
		hss_AA_cart.set_Out_Dir(outDirPCA);
		hss_AA_cart.do_Hierarchical_PCA();

		trace_HAA = hss_AA_cart.get_Trace();
		cond_HAA = hss_AA_cart.get_cond();
		det_HAA = hss_AA_cart.get_Det();
		rank_HAA = hss_AA_cart.get_Rank();
		shrinkage_HAA = hss_AA_cart.get_Shrinkage();

		top_eigenvalues_HAA = hss_AA_cart.getTop_eigenvalues();
		top_evectors_HAA = hss_AA_cart.getTop_evectors();
		square_pca_modes_HAA = hss_AA_cart.getConvoluted_square_pca_modes();
		Matrix wsm = hss_AA_cart.getConvoluted_weighted_square_pca_modes();

		projections_COV = hss_AA_cart.getProjections();
		normed_projections_COV = hss_AA_cart.getNormed_projections_rgc();
		weighted_normed_projections_COV = hss_AA_cart.getWeighted_normed_projections();
		weighted_projections_COV = hss_AA_cart.getWeighted_projections_rgc();

		Convoluted_Eigenvectors_AA_Cart = hss_AA_cart.getConvoluted_Eigenvectors();
		Convoluted_Square_Modes_AA_Cart = hss_AA_cart.getConvoluted_square_pca_modes();
		Convoluted_DVs_AA_Cart = hss_AA_cart.getConvoluted_DVs();

		convoluted_pca_mode_maxes_AA = hss_AA_cart.getConvoluted_pca_mode_max();
		convoluted_pca_mode_mins_AA = hss_AA_cart.getConvoluted_pca_mode_min();

		if (Input_Parameters.doFES) hss_AA_cart.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) hss_AA_cart.do_KPCA(outDirKPCA);

		PC_Plot.createChart4Series(outDirPCA, "Weighted_Projections_Hierarchical", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(outDirPCA, "Top_2_Hierarchical_Weighted_Square_PCA_Modes", wsm, ref_subset_atoms_HAA);
		MODES_Plot.create_XY_Chart_2Series(outDirPCA, "Top_2_Hierarchical_Weighted_Square_PCA_Modes_XY", wsm);
	}

	private static void do_Hierarchical_Heavy_Atom_PCA(String outDirPCA, String outDirFES, String outDirKPCA)
	{
		JEDi_Do_Hierarchical hss_HA_cart = new JEDi_Do_Hierarchical(Input_Parameters.MODES_EIGEN_RESIDUE_HA, Input_Parameters.MODES_HIERARCHICAL_HA,
				Residue_Generalized_Coordinates_Global_HA_Cart, Residue_Eigenvectors_Global_HA, Residue_Delta_Vectors_Global_HA);

		hss_HA_cart.setNumber_of_Atoms(number_of_atoms_hierarchical_HA_SS);
		hss_HA_cart.set_Out_Dir(outDirPCA);
		hss_HA_cart.do_Hierarchical_PCA();

		trace_HHA = hss_HA_cart.get_Trace();
		cond_HHA = hss_HA_cart.get_cond();
		det_HHA = hss_HA_cart.get_Det();
		rank_HHA = hss_HA_cart.get_Rank();
		shrinkage_HHA = hss_HA_cart.get_Shrinkage();

		top_eigenvalues_HHA = hss_HA_cart.getTop_eigenvalues();
		top_evectors_HHA = hss_HA_cart.getTop_evectors();
		square_pca_modes_HHA = hss_HA_cart.getConvoluted_square_pca_modes();
		Matrix wsm = hss_HA_cart.getConvoluted_weighted_square_pca_modes();

		projections_COV = hss_HA_cart.getProjections_rgc();
		normed_projections_COV = hss_HA_cart.getNormed_projections_rgc();
		weighted_normed_projections_COV = hss_HA_cart.getWeighted_normed_projections_rgc();
		weighted_projections_COV = hss_HA_cart.getWeighted_projections_rgc();

		Convoluted_Eigenvectors_HA_Cart = hss_HA_cart.getConvoluted_Eigenvectors();
		Convoluted_Square_Modes_HA_Cart = hss_HA_cart.getConvoluted_square_pca_modes();
		Convoluted_DVs_HA_Cart = hss_HA_cart.getConvoluted_DVs();

		convoluted_pca_mode_maxes_HA = hss_HA_cart.getConvoluted_pca_mode_max();
		convoluted_pca_mode_mins_HA = hss_HA_cart.getConvoluted_pca_mode_min();

		if (Input_Parameters.doFES) hss_HA_cart.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) hss_HA_cart.do_KPCA(outDirKPCA);

		PC_Plot.createChart4Series(outDirPCA, "Weighted_Projections_Hierarchical", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(outDirPCA, "Top_2_Hierarchical_Weighted_Square_PCA_Modes", wsm, ref_subset_atoms_HHA);
		MODES_Plot.create_XY_Chart_2Series(outDirPCA, "Top_2_Hierarchical_Weighted_Square_PCA_Modes_XY", wsm);
	}

	private static void do_Atom_List_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Cartesian cSS_AL = new JEDi_Do_Cartesian(Input_Parameters.MODES_ATOMS_LIST, aligned_subset_PDB_coordinates_AL, aligned_subset_REF_PDB_coordinates_AL, coords, cAL);
		cSS_AL.setNumber_of_residues(number_of_atoms_AL_SS);
		cSS_AL.set_Out_dir(outDirPCA);
		cSS_AL.do_Cartesian();

		trace_AL = cSS_AL.get_Trace_COV();
		cond_AL = cSS_AL.get_Cond_COV();
		det_AL = cSS_AL.get_Det_COV();
		rank_AL = cSS_AL.get_Rank_COV();
		shrinkage_AA = cSS_AL.get_Shrinkage();

		Matrix corr = cSS_AL.getCorr();
		Matrix pcorr = cSS_AL.getPcorr();

		Matrix msa = KMO_MSA.get_MSA(corr, pcorr);
		KMO_AL = KMO_MSA.getKMO();

		top_eigenvalues_COV_AL = cSS_AL.getTop_cartesian_eigenvalues_COV();
		top_evectors_COV_AL = cSS_AL.getTop_cartesian_evectors_COV();
		square_pca_modes_COV_AL = cSS_AL.getSquare_pca_modes_COV();
		Matrix wsm = cSS_AL.getWeighted_square_pca_modes_COV();

		pca_mode_maxes_COV_AL = cSS_AL.getPca_mode_max_COV();
		pca_mode_mins_COV_AL = cSS_AL.getPca_mode_min_COV();

		projections_COV = cSS_AL.getProjections_COV();
		normed_projections_COV = cSS_AL.getNormed_projections_COV();
		weighted_normed_projections_COV = cSS_AL.getWeighted_normed_projections_COV();
		weighted_projections_COV = cSS_AL.getWeighted_projections_COV();

		String out = outDirPCA + "COV" + File.separatorChar;

		if (verbose) Matrix_IO.write_BZ2_Matrix(msa, outDirPCA + "MSA_matrix.txt.bz2", 12, 6);
		STATS_Plot.create_Variables_Stat_Line_Chart(outDirPCA, "MSA Scores", "Atom Number", "MSA", msa, ref_subset_atoms_AL, 1);
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV", wsm, ref_subset_atoms_AL);
		MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV_XY", wsm);

		if (Input_Parameters.doCORR)
			{
				top_eigenvalues_CORR_AL = cSS_AL.getTop_cartesian_eigenvalues_CORR();
				top_evectors_CORR_AL = cSS_AL.getTop_cartesian_evectors_CORR();
				square_pca_modes_CORR_AL = cSS_AL.getSquare_pca_modes_CORR();
				wsm = cSS_AL.getWeighted_square_pca_modes_CORR();
				pca_mode_maxes_CORR_AL = cSS_AL.getPca_mode_max_CORR();
				pca_mode_mins_CORR_AL = cSS_AL.getPca_mode_min_CORR();
				Matrix weighted_projections_CORR = cSS_AL.getWeighted_projections_CORR();

				out = outDirPCA + "CORR" + File.separatorChar;

				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_CORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR", wsm, ref_subset_atoms_AL);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_XY", wsm);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_CORR_SPARSE_AL = cSS_AL.getTop_cartesian_eigenvalues_CORR_SPARSE();
						top_evectors_CORR_SPARSE_AL = cSS_AL.getTop_cartesian_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE_AL = cSS_AL.getSquare_pca_modes_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_AL = cSS_AL.getPca_mode_max_CORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_AL = cSS_AL.getPca_mode_min_CORR_SPARSE();
						wsm = cSS_AL.getWeighted_square_pca_modes_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_AL = cSS_AL.getPca_mode_max_CORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_AL = cSS_AL.getPca_mode_min_CORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = cSS_AL.getWeighted_projections_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", wsm, ref_subset_atoms_AL);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE_XY", wsm);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_eigenvalues_PCORR_AL = cSS_AL.getTop_cartesian_eigenvalues_PCORR();
				top_evectors_PCORR_AL = cSS_AL.getTop_cartesian_evectors_PCORR();
				square_pca_modes_PCORR_AL = cSS_AL.getSquare_pca_modes_PCORR();
				wsm = cSS_AL.getWeighted_square_pca_modes_PCORR();
				pca_mode_maxes_PCORR_AL = cSS_AL.getPca_mode_max_PCORR();
				pca_mode_mins_PCORR_AL = cSS_AL.getPca_mode_min_PCORR();
				Matrix weighted_projections_PCORR = cSS_AL.getWeighted_projections_PCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;

				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_PCORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR", wsm, ref_subset_atoms_AL);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_XY", wsm);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_PCORR_SPARSE_AL = cSS_AL.getTop_cartesian_eigenvalues_PCORR_SPARSE();
						top_evectors_PCORR_SPARSE_AL = cSS_AL.getTop_cartesian_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE_AL = cSS_AL.getSquare_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_AL = cSS_AL.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_AL = cSS_AL.getPca_mode_min_PCORR_SPARSE();
						wsm = cSS_AL.getWeighted_square_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_AL = cSS_AL.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_AL = cSS_AL.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = cSS_AL.getWeighted_projections_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", wsm, ref_subset_atoms_AL);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE_XY", wsm);
					}
			}
		if (Input_Parameters.doFES) cSS_AL.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) cSS_AL.do_KPCA(outDirKPCA);
		do_AL_SSA(outDirSSA);
		// oc.write_PCA_Log(cAL, number_of_atoms_AL_SS, number_of_atoms_AL_SS, Input_Parameters.MODES_ATOMS_LIST, rank_AL, trace_AL, cond_AL, det_AL, KMO_AL, shrinkage_AL);
	}

	private static void do_All_Atom_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Cartesian cSS_AA = new JEDi_Do_Cartesian(Input_Parameters.MODES_ALL_ATOM, aligned_subset_PDB_coordinates_AA, aligned_subset_REF_PDB_coordinates_AA, coords, cAA);

		cSS_AA.setNumber_of_residues(number_of_residues_AA_SS);
		cSS_AA.set_Out_dir(outDirPCA);
		cSS_AA.do_Cartesian();

		trace_AA = cSS_AA.get_Trace_COV();
		cond_AA = cSS_AA.get_Cond_COV();
		det_AA = cSS_AA.get_Det_COV();
		rank_AA = cSS_AA.get_Rank_COV();
		shrinkage_AA = cSS_AA.get_Shrinkage();

		Matrix corr_aa = cSS_AA.getCorr();
		Matrix pcorr_aa = cSS_AA.getPcorr();

		Matrix msa_aa = KMO_MSA.get_MSA(corr_aa, pcorr_aa);
		KMO_AA = KMO_MSA.getKMO();

		top_eigenvalues_COV_AA = cSS_AA.getTop_cartesian_eigenvalues_COV();
		top_evectors_COV_AA = cSS_AA.getTop_cartesian_evectors_COV();
		square_pca_modes_COV_AA = cSS_AA.getSquare_pca_modes_COV();
		Matrix wsm_aa = cSS_AA.getWeighted_square_pca_modes_COV();

		pca_mode_maxes_COV_AA = cSS_AA.getPca_mode_max_COV();
		pca_mode_mins_COV_AA = cSS_AA.getPca_mode_min_COV();

		projections_COV = cSS_AA.getProjections_COV();
		normed_projections_COV = cSS_AA.getNormed_projections_COV();
		weighted_normed_projections_COV = cSS_AA.getWeighted_normed_projections_COV();
		weighted_projections_COV = cSS_AA.getWeighted_projections_COV();

		String out = outDirPCA + "COV" + File.separatorChar;

		if (verbose) Matrix_IO.write_BZ2_Matrix(msa_aa, outDirPCA + "MSA_matrix.txt.bz2", 12, 6);
		STATS_Plot.create_Variables_Stat_Line_Chart(outDirPCA, "MSA Scores", "Atom Number", "MSA", msa_aa, ref_subset_atoms_AA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(outDirPCA, "MSA Scores_XY", "Atom Index", "MSA", msa_aa, ref_subset_atoms_AA, 1);
		MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV", wsm_aa, ref_subset_atoms_AA);
		MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV_XY", wsm_aa);
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_COV);

		if (Input_Parameters.doCORR)
			{
				top_eigenvalues_CORR_AA = cSS_AA.getTop_cartesian_eigenvalues_CORR();
				top_evectors_CORR_AA = cSS_AA.getTop_cartesian_evectors_CORR();
				square_pca_modes_CORR_AA = cSS_AA.getSquare_pca_modes_CORR();
				wsm_aa = cSS_AA.getWeighted_square_pca_modes_CORR();
				pca_mode_maxes_CORR_AA = cSS_AA.getPca_mode_max_CORR();
				pca_mode_mins_CORR_AA = cSS_AA.getPca_mode_min_CORR();
				Matrix weighted_projections_CORR = cSS_AA.getWeighted_projections_CORR();

				out = outDirPCA + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_CORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR", wsm_aa, ref_subset_atoms_AA);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_XY", wsm_aa);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_CORR_SPARSE_AA = cSS_AA.getTop_cartesian_eigenvalues_CORR_SPARSE();
						top_evectors_CORR_SPARSE_AA = cSS_AA.getTop_cartesian_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE_AA = cSS_AA.getSquare_pca_modes_CORR_SPARSE();
						wsm_aa = cSS_AA.getWeighted_square_pca_modes_CORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = cSS_AA.getWeighted_projections_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_AA = cSS_AA.getPca_mode_max_CORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_AA = cSS_AA.getPca_mode_min_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", wsm_aa, ref_subset_atoms_AA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE_XY", wsm_aa);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_eigenvalues_PCORR_AA = cSS_AA.getTop_cartesian_eigenvalues_PCORR();
				top_evectors_PCORR_AA = cSS_AA.getTop_cartesian_evectors_PCORR();
				square_pca_modes_PCORR_AA = cSS_AA.getSquare_pca_modes_PCORR();
				wsm_aa = cSS_AA.getWeighted_square_pca_modes_PCORR();
				pca_mode_maxes_PCORR_AA = cSS_AA.getPca_mode_max_PCORR();
				pca_mode_mins_PCORR_AA = cSS_AA.getPca_mode_min_PCORR();
				Matrix weighted_projections_PCORR = cSS_AA.getWeighted_projections_PCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_PCORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR", wsm_aa, ref_subset_atoms_AA);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_XY", wsm_aa);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_PCORR_SPARSE_AA = cSS_AA.getTop_cartesian_eigenvalues_PCORR_SPARSE();
						top_evectors_PCORR_SPARSE_AA = cSS_AA.getTop_cartesian_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE_AA = cSS_AA.getSquare_pca_modes_PCORR_SPARSE();
						wsm_aa = cSS_AA.getWeighted_square_pca_modes_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = cSS_AA.getWeighted_projections_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_AA = cSS_AA.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_AA = cSS_AA.getPca_mode_min_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", wsm_aa, ref_subset_atoms_AA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE_XY", wsm_aa);
					}
			}
		if (Input_Parameters.doFES) cSS_AA.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) cSS_AA.do_KPCA(outDirKPCA);
		do_AA_SSA(outDirSSA);
	}

	private static void do_Backbone_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Cartesian cSS_BB = new JEDi_Do_Cartesian(Input_Parameters.MODES_BACKBONE, aligned_subset_PDB_coordinates_BB, aligned_subset_REF_PDB_coordinates_BB, coords, cBB);

		cSS_BB.setNumber_of_residues(number_of_residues_BB_SS);
		cSS_BB.set_Out_dir(outDirPCA);
		cSS_BB.do_Cartesian();

		trace_BB = cSS_BB.get_Trace_COV();
		cond_BB = cSS_BB.get_Cond_COV();
		det_BB = cSS_BB.get_Det_COV();
		rank_BB = cSS_BB.get_Rank_COV();
		shrinkage_BB = cSS_BB.get_Shrinkage();

		Matrix corr = cSS_BB.getCorr();
		Matrix pcorr = cSS_BB.getPcorr();

		Matrix msa = KMO_MSA.get_MSA(corr, pcorr);
		KMO_BB = KMO_MSA.getKMO();

		top_eigenvalues_COV_BB = cSS_BB.getTop_cartesian_eigenvalues_COV();
		top_evectors_COV_BB = cSS_BB.getTop_cartesian_evectors_COV();
		square_pca_modes_COV_BB = cSS_BB.getSquare_pca_modes_COV();
		Matrix wsm = cSS_BB.getWeighted_square_pca_modes_COV();
		pca_mode_maxes_COV_BB = cSS_BB.getPca_mode_max_COV();
		pca_mode_mins_COV_BB = cSS_BB.getPca_mode_min_COV();

		projections_COV = cSS_BB.getProjections_COV();
		normed_projections_COV = cSS_BB.getNormed_projections_COV();
		weighted_normed_projections_COV = cSS_BB.getWeighted_normed_projections_COV();
		weighted_projections_COV = cSS_BB.getWeighted_projections_COV();

		String out = outDirPCA + "COV" + File.separatorChar;

		if (verbose) Matrix_IO.write_BZ2_Matrix(msa, outDirPCA + "MSA_matrix.txt.bz2", 12, 6);
		STATS_Plot.create_Variables_Stat_Line_Chart(outDirPCA, "MSA Scores", "Atom Number", "MSA", msa, ref_subset_atoms_BB, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(outDirPCA, "MSA Scores_XY", "Atom Index", "MSA", msa, ref_subset_atoms_BB, 1);
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV", wsm, ref_subset_atoms_BB);
		MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV_XY", wsm);

		if (Input_Parameters.doCORR)
			{
				top_eigenvalues_CORR_BB = cSS_BB.getTop_cartesian_eigenvalues_CORR();
				top_evectors_CORR_BB = cSS_BB.getTop_cartesian_evectors_CORR();
				square_pca_modes_CORR_BB = cSS_BB.getSquare_pca_modes_CORR();
				wsm = cSS_BB.getWeighted_square_pca_modes_CORR();
				pca_mode_maxes_CORR_BB = cSS_BB.getPca_mode_max_CORR();
				pca_mode_mins_CORR_BB = cSS_BB.getPca_mode_min_CORR();
				Matrix weighted_projections_CORR = cSS_BB.getWeighted_projections_CORR();

				out = outDirPCA + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_CORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR", wsm, ref_subset_atoms_BB);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_XY", wsm);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_CORR_SPARSE_BB = cSS_BB.getTop_cartesian_eigenvalues_CORR_SPARSE();
						top_evectors_CORR_SPARSE_BB = cSS_BB.getTop_cartesian_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE_BB = cSS_BB.getSquare_pca_modes_CORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = cSS_BB.getWeighted_projections_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_BB = cSS_BB.getPca_mode_max_CORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_BB = cSS_BB.getPca_mode_min_CORR_SPARSE();
						wsm = cSS_BB.getWeighted_square_pca_modes_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", wsm, ref_subset_atoms_BB);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE_XY", wsm);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_eigenvalues_PCORR_BB = cSS_BB.getTop_cartesian_eigenvalues_PCORR();
				top_evectors_PCORR_BB = cSS_BB.getTop_cartesian_evectors_PCORR();
				square_pca_modes_PCORR_BB = cSS_BB.getSquare_pca_modes_PCORR();
				wsm = cSS_BB.getWeighted_square_pca_modes_PCORR();
				pca_mode_maxes_PCORR_BB = cSS_BB.getPca_mode_max_PCORR();
				pca_mode_mins_PCORR_BB = cSS_BB.getPca_mode_min_PCORR();
				Matrix weighted_projections_PCORR = cSS_BB.getWeighted_projections_PCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_PCORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR", wsm, ref_subset_atoms_BB);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_XY", wsm);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_PCORR_SPARSE_BB = cSS_BB.getTop_cartesian_eigenvalues_PCORR_SPARSE();
						top_evectors_PCORR_SPARSE_BB = cSS_BB.getTop_cartesian_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE_BB = cSS_BB.getSquare_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_BB = cSS_BB.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_BB = cSS_BB.getPca_mode_min_PCORR_SPARSE();
						wsm = cSS_BB.getWeighted_square_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_BB = cSS_BB.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_BB = cSS_BB.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = cSS_BB.getWeighted_projections_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", wsm, ref_subset_atoms_BB);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE_XY", wsm);
					}
			}
		if (Input_Parameters.doFES) cSS_BB.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) cSS_BB.do_KPCA(outDirKPCA);
		do_BB_SSA(outDirSSA);
	}

	private static void do_Heavy_Atom_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Cartesian cSS_HA = new JEDi_Do_Cartesian(Input_Parameters.MODES_HEAVY_ATOM, aligned_subset_PDB_coordinates_HA, aligned_subset_REF_PDB_coordinates_HA, coords, cHA);

		cSS_HA.setNumber_of_residues(number_of_residues_HA_SS);
		cSS_HA.set_Out_dir(outDirPCA);
		cSS_HA.do_Cartesian();

		trace_HA = cSS_HA.get_Trace_COV();
		cond_HA = cSS_HA.get_Cond_COV();
		det_HA = cSS_HA.get_Det_COV();
		rank_HA = cSS_HA.get_Rank_COV();
		shrinkage_HA = cSS_HA.get_Shrinkage();

		Matrix corr_ha = cSS_HA.getCorr();
		Matrix pcorr_ha = cSS_HA.getPcorr();

		Matrix msa_ha = KMO_MSA.get_MSA(corr_ha, pcorr_ha);
		KMO_HA = KMO_MSA.getKMO();

		top_eigenvalues_COV_HA = cSS_HA.getTop_cartesian_eigenvalues_COV();
		top_evectors_COV_HA = cSS_HA.getTop_cartesian_evectors_COV();
		square_pca_modes_COV_HA = cSS_HA.getSquare_pca_modes_COV();
		Matrix wsm_ha = cSS_HA.getWeighted_square_pca_modes_COV();
		pca_mode_maxes_COV_HA = cSS_HA.getPca_mode_max_COV();
		pca_mode_mins_COV_HA = cSS_HA.getPca_mode_min_COV();

		projections_COV = cSS_HA.getProjections_COV();
		normed_projections_COV = cSS_HA.getNormed_projections_COV();
		weighted_normed_projections_COV = cSS_HA.getWeighted_normed_projections_COV();
		weighted_projections_COV = cSS_HA.getWeighted_projections_COV();

		String out = outDirPCA + "COV" + File.separatorChar;

		if (verbose) Matrix_IO.write_BZ2_Matrix(msa_ha, outDirPCA + "MSA_matrix.txt.bz2", 12, 6);
		STATS_Plot.create_Variables_Stat_Line_Chart(outDirPCA, "MSA Scores", "Atom Number", "MSA", msa_ha, ref_subset_atoms_HA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(outDirPCA, "MSA Scores_XY", "Atom Index", "MSA", msa_ha, ref_subset_atoms_HA, 1);
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV", wsm_ha, ref_subset_atoms_HA);
		MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV_XY", wsm_ha);

		if (Input_Parameters.doCORR)
			{
				top_eigenvalues_CORR_HA = cSS_HA.getTop_cartesian_eigenvalues_CORR();
				top_evectors_CORR_HA = cSS_HA.getTop_cartesian_evectors_CORR();
				square_pca_modes_CORR_HA = cSS_HA.getSquare_pca_modes_CORR();
				wsm_ha = cSS_HA.getWeighted_square_pca_modes_CORR();
				pca_mode_maxes_CORR_HA = cSS_HA.getPca_mode_max_CORR();
				pca_mode_mins_CORR_HA = cSS_HA.getPca_mode_min_CORR();
				Matrix weighted_projections_CORR = cSS_HA.getWeighted_projections_CORR();

				out = outDirPCA + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_CORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR", wsm_ha, ref_subset_atoms_HA);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_XY", wsm_ha);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_CORR_SPARSE_HA = cSS_HA.getTop_cartesian_eigenvalues_CORR_SPARSE();
						top_evectors_CORR_SPARSE_HA = cSS_HA.getTop_cartesian_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE_HA = cSS_HA.getSquare_pca_modes_CORR_SPARSE();
						wsm_ha = cSS_HA.getWeighted_square_pca_modes_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_HA = cSS_HA.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_HA = cSS_HA.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = cSS_HA.getWeighted_projections_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", wsm_ha, ref_subset_atoms_HA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE_XY", wsm_ha);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_eigenvalues_PCORR_HA = cSS_HA.getTop_cartesian_eigenvalues_PCORR();
				top_evectors_PCORR_HA = cSS_HA.getTop_cartesian_evectors_PCORR();
				square_pca_modes_PCORR_HA = cSS_HA.getSquare_pca_modes_PCORR();
				wsm_ha = cSS_HA.getWeighted_square_pca_modes_PCORR();
				pca_mode_maxes_PCORR_HA = cSS_HA.getPca_mode_max_PCORR();
				pca_mode_mins_PCORR_HA = cSS_HA.getPca_mode_min_PCORR();
				Matrix weighted_projections_PCORR = cSS_HA.getWeighted_projections_PCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_PCORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR", wsm_ha, ref_subset_atoms_HA);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_XY", wsm_ha);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_PCORR_SPARSE_HA = cSS_HA.getTop_cartesian_eigenvalues_PCORR_SPARSE();
						top_evectors_PCORR_SPARSE_HA = cSS_HA.getTop_cartesian_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE_HA = cSS_HA.getSquare_pca_modes_PCORR_SPARSE();
						wsm_ha = cSS_HA.getWeighted_square_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_HA = cSS_HA.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_HA = cSS_HA.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = cSS_HA.getWeighted_projections_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", wsm_ha, ref_subset_atoms_HA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE_XY", wsm_ha);
					}
			}
		if (Input_Parameters.doFES) cSS_HA.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) cSS_HA.do_KPCA(outDirKPCA);
		do_HA_SSA(outDirSSA);
	}

	private static void do_Alpha_Carbon_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Cartesian cSS_CA = new JEDi_Do_Cartesian(Input_Parameters.MODES_ALPHA_CARBON, aligned_subset_PDB_coordinates_CA, aligned_subset_REF_PDB_coordinates_CA, coords,
				cAC);

		cSS_CA.setNumber_of_residues(number_of_residues_CA_SS);
		cSS_CA.set_Out_dir(outDirPCA);
		cSS_CA.do_Cartesian();

		trace_CA = cSS_CA.get_Trace_COV();
		cond_CA = cSS_CA.get_Cond_COV();
		det_CA = cSS_CA.get_Det_COV();
		rank_CA = cSS_CA.get_Rank_COV();
		shrinkage_CA = cSS_CA.get_Shrinkage();

		Matrix corr = cSS_CA.getCorr();
		Matrix pcorr = cSS_CA.getPcorr();

		Matrix msa = KMO_MSA.get_MSA(corr, pcorr);
		KMO_CA = KMO_MSA.getKMO();

		top_eigenvalues_COV_CA = cSS_CA.getTop_cartesian_eigenvalues_COV();
		top_evectors_COV_CA = cSS_CA.getTop_cartesian_evectors_COV();
		square_pca_modes_COV_CA = cSS_CA.getSquare_pca_modes_COV();
		Matrix wsm = cSS_CA.getWeighted_square_pca_modes_COV();

		pca_mode_maxes_COV_CA = cSS_CA.getPca_mode_max_COV();
		pca_mode_mins_COV_CA = cSS_CA.getPca_mode_min_COV();

		projections_COV = cSS_CA.getProjections_COV();
		normed_projections_COV = cSS_CA.getNormed_projections_COV();
		weighted_normed_projections_COV = cSS_CA.getWeighted_normed_projections_COV();
		weighted_projections_COV = cSS_CA.getWeighted_projections_COV();

		String out = outDirPCA + "COV" + File.separatorChar;

		if (verbose) Matrix_IO.write_BZ2_Matrix(msa, outDirPCA + "MSA_matrix.txt.bz2", 12, 6);
		STATS_Plot.create_Variables_Stat_Line_Chart(outDirPCA, "MSA Scores", "Atom Number", "MSA", msa, ref_subset_atoms_CA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(outDirPCA, "MSA Scores_XY", "Atom Index", "MSA", msa, ref_subset_atoms_CA, 1);
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_COV);
		MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV", wsm, ref_subset_atoms_CA);
		MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_COV_XY", wsm);

		if (Input_Parameters.doCORR)
			{
				top_eigenvalues_CORR_CA = cSS_CA.getTop_cartesian_eigenvalues_CORR();
				top_evectors_CORR_CA = cSS_CA.getTop_cartesian_evectors_CORR();
				square_pca_modes_CORR_CA = cSS_CA.getSquare_pca_modes_CORR();
				wsm = cSS_CA.getWeighted_square_pca_modes_CORR();
				pca_mode_maxes_CORR_CA = cSS_CA.getPca_mode_max_CORR();
				pca_mode_mins_CORR_CA = cSS_CA.getPca_mode_min_CORR();
				Matrix weighted_projections_CORR = cSS_CA.getWeighted_projections_CORR();

				out = outDirPCA + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_CORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR", wsm, ref_subset_atoms_CA);
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_XY", wsm);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_CORR_SPARSE_CA = cSS_CA.getTop_cartesian_eigenvalues_CORR_SPARSE();
						top_evectors_CORR_SPARSE_CA = cSS_CA.getTop_cartesian_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE_CA = cSS_CA.getSquare_pca_modes_CORR_SPARSE();
						wsm = cSS_CA.getWeighted_square_pca_modes_CORR_SPARSE();
						pca_mode_maxes_CORR_SPARSE_CA = cSS_CA.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_CORR_SPARSE_CA = cSS_CA.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = cSS_CA.getWeighted_projections_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE", wsm, ref_subset_atoms_CA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_CORR_SPARSE_XY", wsm);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_eigenvalues_PCORR_CA = cSS_CA.getTop_cartesian_eigenvalues_PCORR();
				top_evectors_PCORR_CA = cSS_CA.getTop_cartesian_evectors_PCORR();
				square_pca_modes_PCORR_CA = cSS_CA.getSquare_pca_modes_PCORR();
				wsm = cSS_CA.getWeighted_square_pca_modes_PCORR();
				pca_mode_maxes_PCORR_CA = cSS_CA.getPca_mode_max_PCORR();
				pca_mode_mins_PCORR_CA = cSS_CA.getPca_mode_min_PCORR();

				Matrix weighted_projections_PCORR = cSS_CA.getWeighted_projections_PCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;

				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_PCORR);
				MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR", wsm, ref_subset_atoms_CA);

				if (Input_Parameters.doSPARSIFY)
					{
						top_eigenvalues_PCORR_SPARSE_CA = cSS_CA.getTop_cartesian_eigenvalues_PCORR_SPARSE();
						top_evectors_PCORR_SPARSE_CA = cSS_CA.getTop_cartesian_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE_CA = cSS_CA.getSquare_pca_modes_PCORR_SPARSE();
						wsm = cSS_CA.getWeighted_square_pca_modes_PCORR_SPARSE();
						pca_mode_maxes_PCORR_SPARSE_CA = cSS_CA.getPca_mode_max_PCORR_SPARSE();
						pca_mode_mins_PCORR_SPARSE_CA = cSS_CA.getPca_mode_min_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = cSS_CA.getWeighted_projections_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.create_Line_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE", wsm, ref_subset_atoms_CA);
						MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_SPARSE_XY", wsm);
					}
				MODES_Plot.create_XY_Chart_2Series(out, "Top_2_Weighted_Square_PCA_Modes_PCORR_XY", wsm);
			}
		if (Input_Parameters.doFES) cSS_CA.do_FES(outDirFES);
		if (Input_Parameters.doKPCA) cSS_CA.do_KPCA(outDirKPCA);
		do_CA_SSA(outDirSSA);
	}

	private static void do_Distance_Pair_PCA(String outDirPCA, String outDirSSA, String outDirFES, String outDirKPCA, Matrix coords)
	{
		JEDi_Do_Dist_Pairs dpSS = new JEDi_Do_Dist_Pairs(coords, reference_distances, outDirPCA, outDirSSA, outDirFES, outDirKPCA);

		dpSS.do_Dist();

		trace_d = dpSS.getTrace_dist_COV();
		cond_d = dpSS.getCond_cov();
		det_d = dpSS.getDet_cov();
		rank_d = dpSS.getRank_cov();
		shrinkage_d = dpSS.getShrinkage();
		KMO_DP = dpSS.getKMO();

		top_distance_evectors_COV = dpSS.getTop_distance_evectors_COV();
		projections_dist_COV = dpSS.getProjections_dist_COV();
		normed_projections_dist_COV = dpSS.getNormed_projections_dist_COV();
		weighted_projections_dist_COV = dpSS.getWeighted_projections_dist_COV();
		weighted_normed_projections_dist_COV = dpSS.getWeighted_normed_projections_dist_COV();
		weighted_square_pca_modes_COV_dist = dpSS.getWeightedSquarePCAmodesCOV();

		Matrix msa = dpSS.getMSA();
		String name = "MSA_Scores_Matrix.txt.bz2";
		if (verbose) Matrix_IO.write_BZ2_Matrix(msa, outDirPCA + name, 12, 6);
		MODES_Plot.createLineChart1SeriesDist(outDirPCA, "MSA Scores", "Distance Pairs", "MSA", msa, atom_list_dp_original1, atom_list_dp_original2, 1);

		String out = outDirPCA + "COV" + File.separatorChar;
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_dist_COV);

		String title = "Top_2_Weighted_Square_PCA_Modes_COV";
		MODES_Plot.createLineChart2SeriesDist(out, title, weighted_square_pca_modes_COV_dist, ref_subset_atoms_DP1, ref_subset_atoms_DP2);

		if (Input_Parameters.doCORR)
			{
				top_distance_evectors_CORR = dpSS.getTop_distance_evectors_CORR();
				projections_dist_CORR = dpSS.getProjections_dist_CORR();
				normed_projections_dist_CORR = dpSS.getNormed_projections_dist_CORR();
				weighted_projections_dist_CORR = dpSS.getWeighted_projections_dist_CORR();
				weighted_normed_projections_dist_CORR = dpSS.getWeighted_normed_projections_dist_CORR();
				weighted_square_pca_modes_CORR_dist = dpSS.getWeightedSquarePCAmodesCORR();

				out = outDirPCA + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_dist_CORR);
				title = "Top_2_Weighted_Square_PCA_Modes_CORR_Line";
				MODES_Plot.createLineChart2SeriesDist(out, title, weighted_square_pca_modes_CORR_dist, ref_subset_atoms_DP1, ref_subset_atoms_DP2);

				if (Input_Parameters.doSPARSIFY)
					{
						Matrix wsm = dpSS.getWeighted_Square_PCA_modes_CORR_SPARSE();
						Matrix weighted_projections_CORR_SPARSE = dpSS.getWeighted_projections_dist_CORR_SPARSE();

						out = outDirPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar;

						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_CORR_SPARSE);
						MODES_Plot.createLineChart2SeriesDist(out, title, wsm, ref_subset_atoms_DP1, ref_subset_atoms_DP2);
					}
			}
		if (Input_Parameters.doPCORR)
			{
				top_distance_evectors_PCORR = dpSS.getTop_distance_evectors_PCORR();
				projections_dist_PCORR = dpSS.getProjections_dist_PCORR();
				normed_projections_dist_PCORR = dpSS.getNormed_projections_dist_PCORR();
				weighted_projections_dist_PCORR = dpSS.getWeighted_projections_dist_PCORR();
				weighted_normed_projections_dist_PCORR = dpSS.getWeighted_normed_projections_dist_PCORR();
				weighted_square_pca_modes_PCORR_dist = dpSS.getWeightedSquarePCAmodesPCORR();

				out = outDirPCA + "PCORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_dist_PCORR);
				title = "Top_2_Weighted_Square_PCA_Modes_PCORR_Line";
				MODES_Plot.createLineChart2SeriesDist(out, title, weighted_square_pca_modes_PCORR_dist, ref_subset_atoms_DP1, ref_subset_atoms_DP2);

				if (Input_Parameters.doSPARSIFY)
					{
						Matrix wsm = dpSS.getWeighted_Square_PCA_modes_PCORR_SPARSE();
						Matrix weighted_projections_PCORR_SPARSE = dpSS.getWeighted_projections_dist_PCORR_SPARSE();

						out = outDirPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_PCORR_SPARSE);
						MODES_Plot.createLineChart2SeriesDist(out, title, wsm, ref_subset_atoms_DP1, ref_subset_atoms_DP2);
					}
			}
	}

	public static void do_PreProcess()
	{
		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_coordinates_AA, original_reference_PDB_coordinates_AA);
		tf_coords.get_Aligned_Reference_Coordinates();
		Matrix aligned_coords = tf_coords.get_SS_Transformed_coords();
		tf_coords.get_SS_coordinate_STATS();

		outputAlignedCoords = Input_Parameters.doOutputCoordinates;
		if (outputAlignedCoords)
			{
				name = "Aligned_Original_PDB_Coordinates_" + number_of_residues_REF + "_" + number_of_atoms_REF + ".txt.bz2";
				path = Input_Parameters.OUT_DIR + name;
				Matrix_IO.write_BZ2_Matrix(aligned_coords, path, 9, 3);
			}

		List<Double> conf_rmsds = tf_coords.get_SS_Conformation_RMSDs();
		List<Double> aligned_atomic_RMSFs = tf_coords.get_SS_RMSF();

		List<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSFs(ref_atoms_AA, aligned_atomic_RMSFs);
		String ref_PDB_file_name = "ss_" + number_of_residues_REF + "_Atomic_RMSF_edited.pdb";
		path = Input_Parameters.OUT_DIR + ref_PDB_file_name;
		PDB_IO.Write_PDB(path, edited_atoms);

		Matrix stats = tf_coords.getStats();

		verbose = Input_Parameters.verbose;
		if (verbose)
			{
				name = "Original_PDB_Coordinate_Stats.txt.bz2";
				path = Input_Parameters.OUT_DIR + name;
				Matrix_IO.write_BZ2_Matrix(stats, path, 12, 6);
			}

		RMSD_Plot.create_RMSD_XY_Chart(Input_Parameters.OUT_DIR, "PreProcess_Conformation_RMSDs", "Conformation Number", conf_rmsds);
		RMSD_Plot.create_RMSF_Line_Chart(Input_Parameters.OUT_DIR, "PreProcess_Atomic_RMSFs", "Atom Number", aligned_atomic_RMSFs, ref_atoms_AA);
		RMSD_Plot.create_RMSD_XY_Chart(Input_Parameters.OUT_DIR, "PreProcess_Atomic_RMSFs_XY", "Atom Index", aligned_atomic_RMSFs);

		STATS_Plot.create_Variables_Stat_Line_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Means", "Atom Number", "Mean", stats, ref_atoms_AA, 1);
		STATS_Plot.create_Variables_Stat_Line_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Variances", "Atom Number", "Variance", stats, ref_atoms_AA, 2);
		STATS_Plot.create_Variables_Stat_Line_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Skews", "Atom Number", "Skew", stats, ref_atoms_AA, 3);
		STATS_Plot.create_Variables_Stat_Line_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Kurtosis", "Atom Number", "Kurtosis", stats, ref_atoms_AA, 4);

		STATS_Plot.create_Variables_Stat_XY_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Means_XY", "Atom Index", "Mean", stats, ref_atoms_AA, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Variances_XY", "Atom Index", "Variance", stats, ref_atoms_AA, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Skews_XY", "Atom Index", "Skew", stats, ref_atoms_AA, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(Input_Parameters.OUT_DIR, "Original_PDB_Coordinate_Variable_Kurtosis_XY", "Atom Index", "Kurtosis", stats, ref_atoms_AA, 4);

		doStatThresholds = Input_Parameters.do_StatThresholds;
		if (doStatThresholds)
			{
				Select_Variables_by_Statistical_Threshold sv = new Select_Variables_by_Statistical_Threshold(PP, ref_atoms_AA, stats);
				sv.doThresholding();
				oc.write_Stat_Thresholding_Log(PP);
			}
		oc.write_PreProcess_Log(number_of_residues_REF, number_of_atoms_REF, number_conformations);
	}

	private static void do_Mode_Visualization_AL(String type)
	{
		JEDi_Get_PCA_Mode_Vizualization cov = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AL, top_eigenvalues_COV_AL, top_evectors_COV_AL, square_pca_modes_COV_AL,
				pca_mode_maxes_COV_AL, pca_mode_mins_COV_AL, Input_Parameters.MODES_VIZ, type, Q);

		cov.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + Q + File.separatorChar);

		if (Input_Parameters.doModeViz) cov.get_Mode_Visualizations_All_Atom();
		cov.get_Essential_Visualization_All_Atom();

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization corr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AL, top_eigenvalues_CORR_AL, top_evectors_CORR_AL,
						square_pca_modes_CORR_AL, pca_mode_maxes_CORR_AL, pca_mode_mins_CORR_AL, Input_Parameters.MODES_VIZ, type, R);

				corr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar);

				if (Input_Parameters.doModeViz) corr.get_Mode_Visualizations_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AL, top_eigenvalues_CORR_SPARSE_AL,
								top_evectors_CORR_SPARSE_AL, square_pca_modes_CORR_SPARSE_AL, pca_mode_maxes_CORR_SPARSE_AL, pca_mode_mins_CORR_SPARSE_AL,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}

		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization pcorr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AL, top_eigenvalues_PCORR_AL, top_evectors_PCORR_AL,
						square_pca_modes_PCORR_AL, pca_mode_maxes_PCORR_AL, pca_mode_mins_PCORR_AL, Input_Parameters.MODES_VIZ, type, P);

				pcorr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar);

				if (Input_Parameters.doModeViz) pcorr.get_Mode_Visualizations_All_Atom();
				pcorr.get_Essential_Visualization_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AL, top_eigenvalues_PCORR_SPARSE_AL,
								top_evectors_PCORR_SPARSE_AL, square_pca_modes_PCORR_SPARSE_AL, pca_mode_maxes_PCORR_SPARSE_AL, pca_mode_mins_PCORR_SPARSE_AL,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}
	}

	private static void do_Mode_Visualization_AA(String type)
	{
		JEDi_Get_PCA_Mode_Vizualization cov = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AA, top_eigenvalues_COV_AA, top_evectors_COV_AA, square_pca_modes_COV_AA,
				pca_mode_maxes_COV_AA, pca_mode_mins_COV_AA, Input_Parameters.MODES_VIZ, type, Q);

		cov.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + Q + File.separatorChar);

		if (Input_Parameters.doModeViz) cov.get_Mode_Visualizations_All_Atom();
		cov.get_Essential_Visualization_All_Atom();

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization corr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AA, top_eigenvalues_CORR_AA, top_evectors_CORR_AA,
						square_pca_modes_CORR_AA, pca_mode_maxes_CORR_AA, pca_mode_mins_CORR_AA, Input_Parameters.MODES_VIZ, type, R);

				corr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar);

				if (Input_Parameters.doModeViz) corr.get_Mode_Visualizations_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AA, top_eigenvalues_CORR_SPARSE_AA,
								top_evectors_CORR_SPARSE_AA, square_pca_modes_CORR_SPARSE_AA, pca_mode_maxes_CORR_SPARSE_AA, pca_mode_mins_CORR_SPARSE_AA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}

		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization pcorr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AA, top_eigenvalues_PCORR_AA, top_evectors_PCORR_AA,
						square_pca_modes_PCORR_AA, pca_mode_maxes_PCORR_AA, pca_mode_mins_PCORR_AA, Input_Parameters.MODES_VIZ, type, P);

				pcorr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar);

				if (Input_Parameters.doModeViz) pcorr.get_Mode_Visualizations_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_AA, top_eigenvalues_PCORR_SPARSE_AA,
								top_evectors_PCORR_SPARSE_AA, square_pca_modes_PCORR_SPARSE_AA, pca_mode_maxes_PCORR_SPARSE_AA, pca_mode_mins_PCORR_SPARSE_AA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
				pcorr.get_Essential_Visualization_All_Atom();
			}
	}

	private static void do_Mode_Visualization_BB(String type)
	{
		JEDi_Get_PCA_Mode_Vizualization cov = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_BB, top_eigenvalues_COV_BB, top_evectors_COV_BB, square_pca_modes_COV_BB,
				pca_mode_maxes_COV_BB, pca_mode_mins_COV_BB, Input_Parameters.MODES_VIZ, type, Q);

		cov.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + Q + File.separatorChar);

		if (Input_Parameters.doModeViz) cov.get_Mode_Visualizations_All_Atom();
		cov.get_Essential_Visualization_All_Atom();

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization corr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_BB, top_eigenvalues_CORR_BB, top_evectors_CORR_BB,
						square_pca_modes_CORR_BB, pca_mode_maxes_CORR_BB, pca_mode_mins_CORR_BB, Input_Parameters.MODES_VIZ, type, R);

				corr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar);

				if (Input_Parameters.doModeViz) corr.get_Mode_Visualizations_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_BB, top_eigenvalues_CORR_SPARSE_BB,
								top_evectors_CORR_SPARSE_BB, square_pca_modes_CORR_SPARSE_BB, pca_mode_maxes_CORR_SPARSE_BB, pca_mode_mins_CORR_SPARSE_BB,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}
		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization pcorr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_BB, top_eigenvalues_PCORR_BB, top_evectors_PCORR_BB,
						square_pca_modes_PCORR_BB, pca_mode_maxes_PCORR_BB, pca_mode_mins_PCORR_BB, Input_Parameters.MODES_VIZ, type, P);

				pcorr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar);

				if (Input_Parameters.doModeViz) pcorr.get_Mode_Visualizations_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_BB, top_eigenvalues_PCORR_SPARSE_BB,
								top_evectors_PCORR_SPARSE_BB, square_pca_modes_PCORR_SPARSE_BB, pca_mode_maxes_PCORR_SPARSE_BB, pca_mode_mins_PCORR_SPARSE_BB,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
				pcorr.get_Essential_Visualization_All_Atom();
			}
	}

	private static void do_Mode_Visualization_CA(String type)
	{
		/* Check for OXT, OC2 in amended atom list and remove if there. Need to check each chain... */
		int index = 0;
		boolean oxt = false, oc2 = false;
		List<Integer> removeList = new ArrayList<Integer>(10);
		for (Atom a : amended_subset_atoms_CA_BB_VIZ)
			{
				if (a.symbol.equals("OXT"))
					{
						oxt = true;
						index = amended_subset_atoms_CA_BB_VIZ.indexOf(a);
						removeList.add(index);
						// System.out.println("\t\tOXT detected... Removing from the amended set of atoms...");
					}
				if (a.symbol.equals("OC2"))
					{
						oc2 = true;
						index = amended_subset_atoms_CA_BB_VIZ.indexOf(a);
						removeList.add(index);
						// System.out.println("\t\tOC2 detected... Removing from the amended set of atoms...");
					}
			}
		if (oxt || oc2)
			{
				int size = removeList.size();
				int last = size - 1;
				for (int i = 0; i < size; i++)
					{
						int x = removeList.get(last - i);
						amended_subset_atoms_CA_BB_VIZ.remove(x);
					}
			}

		JEDi_Get_PCA_Mode_Vizualization cov = new JEDi_Get_PCA_Mode_Vizualization(amended_subset_atoms_CA_BB_VIZ, top_eigenvalues_COV_CA, top_evectors_COV_CA,
				square_pca_modes_COV_CA, pca_mode_maxes_COV_CA, pca_mode_mins_COV_CA, Input_Parameters.MODES_VIZ, type, Q);

		cov.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + Q + File.separatorChar);

		if (Input_Parameters.doModeViz) cov.get_Mode_Visualizations_BackBone();
		cov.get_Essential_Visualization_BackBone();

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization corr = new JEDi_Get_PCA_Mode_Vizualization(amended_subset_atoms_CA_BB_VIZ, top_eigenvalues_CORR_CA, top_evectors_CORR_CA,
						square_pca_modes_CORR_CA, pca_mode_maxes_CORR_CA, pca_mode_mins_CORR_CA, Input_Parameters.MODES_VIZ, type, R);

				corr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar);

				if (Input_Parameters.doModeViz) corr.get_Mode_Visualizations_BackBone();
				corr.get_Essential_Visualization_BackBone();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_CA, top_eigenvalues_CORR_SPARSE_CA,
								top_evectors_CORR_SPARSE_CA, square_pca_modes_CORR_SPARSE_CA, pca_mode_maxes_CORR_SPARSE_CA, pca_mode_mins_CORR_SPARSE_CA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}
		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization pcorr = new JEDi_Get_PCA_Mode_Vizualization(amended_subset_atoms_CA_BB_VIZ, top_eigenvalues_PCORR_CA, top_evectors_PCORR_CA,
						square_pca_modes_PCORR_CA, pca_mode_maxes_PCORR_CA, pca_mode_mins_PCORR_CA, Input_Parameters.MODES_VIZ, type, P);

				pcorr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar);

				if (Input_Parameters.doModeViz) pcorr.get_Mode_Visualizations_BackBone();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_CA, top_eigenvalues_PCORR_SPARSE_CA,
								top_evectors_PCORR_SPARSE_CA, square_pca_modes_PCORR_SPARSE_CA, pca_mode_maxes_PCORR_SPARSE_CA, pca_mode_mins_PCORR_SPARSE_CA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
				pcorr.get_Essential_Visualization_BackBone();
			}
	}

	private static void do_Mode_Visualization_HA(String type)
	{
		JEDi_Get_PCA_Mode_Vizualization cov = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HA, top_eigenvalues_COV_HA, top_evectors_COV_HA, square_pca_modes_COV_HA,
				pca_mode_maxes_COV_HA, pca_mode_mins_COV_HA, Input_Parameters.MODES_VIZ, type, Q);

		cov.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + Q + File.separatorChar);

		if (Input_Parameters.doModeViz) cov.get_Mode_Visualizations_All_Atom();
		cov.get_Essential_Visualization_All_Atom();

		if (Input_Parameters.doCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization corr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HA, top_eigenvalues_CORR_HA, top_evectors_CORR_HA,
						square_pca_modes_CORR_HA, pca_mode_maxes_CORR_HA, pca_mode_mins_CORR_HA, Input_Parameters.MODES_VIZ, type, R);

				corr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar);

				if (Input_Parameters.doModeViz) corr.get_Mode_Visualizations_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HA, top_eigenvalues_CORR_SPARSE_HA,
								top_evectors_CORR_SPARSE_HA, square_pca_modes_CORR_SPARSE_HA, pca_mode_maxes_CORR_SPARSE_HA, pca_mode_mins_CORR_SPARSE_HA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + R + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
			}

		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_PCA_Mode_Vizualization pcorr = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HA, top_eigenvalues_PCORR_HA, top_evectors_PCORR_HA,
						square_pca_modes_PCORR_HA, pca_mode_maxes_PCORR_HA, pca_mode_mins_PCORR_HA, Input_Parameters.MODES_VIZ, type, P);

				pcorr.set_Output_Directory(Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar);

				if (Input_Parameters.doModeViz) pcorr.get_Mode_Visualizations_All_Atom();

				if (Input_Parameters.doSPARSIFY)
					{
						JEDi_Get_PCA_Mode_Vizualization corrS = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HA, top_eigenvalues_PCORR_SPARSE_HA,
								top_evectors_PCORR_SPARSE_HA, square_pca_modes_PCORR_SPARSE_HA, pca_mode_maxes_PCORR_SPARSE_HA, pca_mode_mins_PCORR_SPARSE_HA,
								Input_Parameters.MODES_VIZ, type, R);

						corrS.set_Output_Directory(
								Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + type + File.separatorChar + P + File.separatorChar + S + File.separatorChar);

						if (Input_Parameters.doModeViz) corrS.get_Mode_Visualizations_All_Atom();
						corrS.get_Essential_Visualization_All_Atom();
					}
				pcorr.get_Essential_Visualization_All_Atom();
			}
	}

	private static void do_Mode_Visualization_HAA_Cart()
	{
		JEDi_Get_PCA_Mode_Vizualization covAA = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HAA, top_eigenvalues_HAA, Convoluted_Eigenvectors_AA_Cart,
				Convoluted_Square_Modes_AA_Cart, convoluted_pca_mode_maxes_AA, convoluted_pca_mode_mins_AA, Input_Parameters.MODES_VIZ, HAA, Q);

		covAA.set_Output_Directory(
				Input_Parameters.DIRECTORY + "JEDi_RESULTS_" + Input_Parameters.DESCRIPTION + File.separatorChar + "VIZ" + File.separatorChar + HAA + File.separatorChar);

		if (Input_Parameters.doModeViz) covAA.get_Mode_Visualizations_All_Atom();
		covAA.get_Essential_Visualization_All_Atom();
	}

	private static void do_Mode_Visualization_HHA_Cart()
	{
		JEDi_Get_PCA_Mode_Vizualization covHA = new JEDi_Get_PCA_Mode_Vizualization(ref_subset_atoms_HHA, top_eigenvalues_HHA, Convoluted_Eigenvectors_HA_Cart,
				Convoluted_Square_Modes_HA_Cart, convoluted_pca_mode_maxes_HA, convoluted_pca_mode_mins_HA, Input_Parameters.MODES_VIZ, HHA, Q);

		covHA.set_Output_Directory(
				Input_Parameters.DIRECTORY + "JEDi_RESULTS_" + Input_Parameters.DESCRIPTION + File.separatorChar + "VIZ" + File.separatorChar + HHA + File.separatorChar);

		if (Input_Parameters.doModeViz) covHA.get_Mode_Visualizations_All_Atom();
		covHA.get_Essential_Visualization_All_Atom();
	}

	private static void do_Comparative_AA_SSA(String outSSA, Matrix evects1, Matrix evects2)
	{
		if (Input_Parameters.do_cartesian_AA
				& (number_of_atoms_hierarchical_AA_SS == number_of_atoms_AA_SS && Input_Parameters.MODES_ALL_ATOM == Input_Parameters.MODES_HIERARCHICAL_AA))
			{
				if (Convoluted_Eigenvectors_AA_Cart.getRowDimension() == top_evectors_COV_AA.getRowDimension()
						& Convoluted_Eigenvectors_AA_Cart.getColumnDimension() == top_evectors_COV_AA.getColumnDimension())
					{
						doInterSSA_AA = true;
						JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(AA, "Hierarchical_vs_Direct_Cartesian", evects1, evects2);
						ssa.setOut_dir(outSSA);
						ssa.get_SSA();
						ssa.get_FSSA_Iterated();
					}
			}
	}

	private static void do_Comparative_HA_SSA(String outSSA, Matrix evects1, Matrix evects2)
	{
		if (Input_Parameters.do_cartesian_HA
				& (number_of_atoms_hierarchical_HA_SS == number_of_atoms_HA_SS & Input_Parameters.MODES_HEAVY_ATOM == Input_Parameters.MODES_HIERARCHICAL_HA))
			{

				if (Convoluted_Eigenvectors_HA_Cart.getRowDimension() == top_evectors_COV_HA.getRowDimension()
						& Convoluted_Eigenvectors_HA_Cart.getColumnDimension() == top_evectors_COV_HA.getColumnDimension())
					{
						doInterSSA_HA = true;
						JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(HA, "Hierarchical_vs_Direct_Cartesian", evects1, evects2);
						ssa.setOut_dir(outSSA);
						ssa.get_SSA();
						ssa.get_FSSA_Iterated();
					}
			}
	}

	private static void do_AA_SSA(String outSSA)
	{
		JEDi_Get_Subspace_Analysis ssa;
		if (Input_Parameters.doCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAA, "COV_vs_CORR", top_evectors_CORR_AA, top_evectors_COV_AA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAA, "COV_vs_PCORR", top_evectors_PCORR_AA, top_evectors_COV_AA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAA, "CORR_vs_PCORR", top_evectors_PCORR_AA, top_evectors_CORR_AA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAA, "CORR_vs_CORR_SPARSE", top_evectors_CORR_AA, top_evectors_CORR_SPARSE_AA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAA, "PCORR_vs_PCORR_SPARSE", top_evectors_PCORR_AA, top_evectors_PCORR_SPARSE_AA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	private static void do_AL_SSA(String outSSA)
	{
		JEDi_Get_Subspace_Analysis ssa;
		if (Input_Parameters.doCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAL, "COV_vs_CORR", top_evectors_CORR_AL, top_evectors_COV_AL);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAL, "COV_vs_PCORR", top_evectors_PCORR_AL, top_evectors_COV_AL);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAL, "CORR_vs_PCORR", top_evectors_PCORR_AL, top_evectors_CORR_AL);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();

			}
		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAL, "CORR_vs_CORR_SPARSE", top_evectors_CORR_AL, top_evectors_CORR_SPARSE_AL);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAL, "PCORR_vs_PCORR_SPARSE", top_evectors_PCORR_AL, top_evectors_PCORR_SPARSE_AL);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	private static void do_HA_SSA(String outSSA)
	{
		JEDi_Get_Subspace_Analysis ssa;
		if (Input_Parameters.doCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cHA, "COV_vs_CORR", top_evectors_CORR_HA, top_evectors_COV_HA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cHA, "COV_vs_PCORR", top_evectors_PCORR_HA, top_evectors_COV_HA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cHA, "CORR_vs_PCORR", top_evectors_PCORR_HA, top_evectors_CORR_HA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();

			}
		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cHA, "CORR_vs_CORR_SPARSE", top_evectors_CORR_HA, top_evectors_CORR_SPARSE_HA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cHA, "PCORR_vs_PCORR_SPARSE", top_evectors_PCORR_HA, top_evectors_PCORR_SPARSE_HA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	private static void do_BB_SSA(String outSSA)
	{
		JEDi_Get_Subspace_Analysis ssa;
		if (Input_Parameters.doCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cBB, "COV_vs_CORR", top_evectors_CORR_BB, top_evectors_COV_BB);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cBB, "COV_vs_PCORR", top_evectors_PCORR_BB, top_evectors_COV_BB);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cBB, "CORR_vs_PCORR", top_evectors_PCORR_BB, top_evectors_CORR_BB);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}

		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cBB, "CORR_vs_CORR_SPARSE", top_evectors_CORR_BB, top_evectors_CORR_SPARSE_BB);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cBB, "PCORR_vs_PCORR_SPARSE", top_evectors_PCORR_BB, top_evectors_PCORR_SPARSE_BB);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	private static void do_CA_SSA(String outSSA)
	{
		JEDi_Get_Subspace_Analysis ssa;
		if (Input_Parameters.doCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAC, "COV_vs_CORR", top_evectors_CORR_CA, top_evectors_COV_CA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAC, "COV_vs_PCORR", top_evectors_PCORR_CA, top_evectors_COV_CA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAC, "CORR_vs_PCORR", top_evectors_PCORR_CA, top_evectors_CORR_CA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}

		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAC, "CORR_vs_CORR_SPARSE", top_evectors_CORR_CA, top_evectors_CORR_SPARSE_CA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				ssa = new JEDi_Get_Subspace_Analysis(cAC, "PCORR_vs_PCORR_SPARSE", top_evectors_PCORR_CA, top_evectors_PCORR_SPARSE_CA);
				ssa.setOut_dir(outSSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	// ******************************************** MAIN METHOD ************************************************************* //

	public static void main(String[] args)
	{
		jedi = new JEDi_Driver_MT();
		jediStart = System.nanoTime();
		System.out.println("Running JEDi Driver - Multi Threaded: \n");
		try
			{
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_file;
				if (args.length == 0)
					{
						System.out.println("\tUsing the default input file: 'JEDi_Parameters.txt'\n");
					}
				if (args.length == 1)
					{
						in_path = args[0];
						System.out.println("\tUsing a specified input file\n");
					}
				System.out.println("\tCurrent User Directory (calling the JVM) = " + WD + "\n");
				System.out.println("\tFull Path to Input File = " + in_path + "\n");
				boolean check = new File(in_path).exists();
				if (!check)
					{
						System.err.println("The Input File can not be found: " + in_path);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				input_reader = new BufferedReader(new FileReader(in_path));
			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the input file: " + input_file);
				System.err.println("Program terminating.\n");
				e.printStackTrace();
				System.exit(0);
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		System.out.println("Reading the input file: \n");
		read_input_file();
		System.out.println("Processing Parameters from Input File... \n");
		ip = new Input_Parameters(parameters);
		oc = new Output_Control();
		oc.initialize_JED_Log();
		oc.write_Parameter_Log(lines);
		System.out.println("Processing PDB Reference File... \n");
		startTime = System.nanoTime();
		process_Reference_PDB();
		endTime = System.nanoTime();
		totalTime = endTime - startTime;
		System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds) \n");
		System.out.println("----------------------------------------------------------------------------------------------\n");
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.doPREPROCESS)
			{
				System.out.println("PROCESSING PDB files to create matrix of coordinates for JEDi Analysis... \n");
				System.out.println("\tReading PDB files... \n");
				startTime = System.nanoTime();
				read_PDB_Files();
				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
				System.out.println("Computing preprocess statistics... \n");
				startTime = System.nanoTime();
				do_PreProcess();
				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
				System.out.println("----------------------------------------------------------------------------------------------\n");
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (!Input_Parameters.doPREPROCESS)
			{
				System.out.println("Reading ALL ATOM Input Coordinate File... \n");
				read_All_Atom_Coordinate_File();

				if (Input_Parameters.doDownSample)
					{
						System.out.println("Down Sampling the frames in the input file... \n");
						do_Down_Sample();
					}

				if (Input_Parameters.doFrameSelect)
					{
						System.out.println("Selecting frames... \n");
						do_Frame_Select();
					}

				System.out.println("Getting Subsets and Aligning where required... \n");
				get_Subsets();
				align_Subsets();
				System.out.println("----------------------------------------------------------------------------------------------\n");
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_residue_pairs)
			{
				future0 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						do_All_Atom_Residue_Pairs_Analysis();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed RESIDUE PAIRS Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_residue_individual)
			{
				future1 = service.submit(new Runnable()
				{

					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA, outVIZ;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + Res_AA + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + Res_AA + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + Res_AA + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + Res_AA + File.separatorChar + "INLIERS" + File.separatorChar;
								outVIZ = Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + Res_AA + File.separatorChar + "INLIERS" + File.separatorChar;

								do_Individual_All_Atom_Residue_PCA(outPCA, outSSA, outFES, outKPCA, outVIZ, aligned_subset_PDB_coordinates_Outliers_REMOVED_local);

								outPCA = Input_Parameters.OUT_DIR + Res_AA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + Res_AA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + Res_AA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + Res_AA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outVIZ = Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + Res_AA + File.separatorChar + "OUTLIERS" + File.separatorChar;

								do_Individual_All_Atom_Residue_PCA(outPCA, outSSA, outFES, outKPCA, outVIZ, aligned_subset_PDB_coordinates_Outliers_SELECTED_local);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + Res_AA + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + Res_AA + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + Res_AA + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + Res_AA + File.separatorChar;
								outVIZ = Input_Parameters.OUT_DIR + "VIZ" + File.separatorChar + Res_AA + File.separatorChar;
								do_Individual_All_Atom_Residue_PCA(outPCA, outSSA, outFES, outKPCA, outVIZ, aligned_subset_PDB_coordinates_local);
							}
						oc.write_Individual_Residue_PCA_Log();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed INDIVIDUAL RESIDUE PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_hierarchical_AA)
			{
				future2 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outFES, outKPCA;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + HAA + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HAA + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HAA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Global_All_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_Outliers_REMOVED_AA);
								do_Hierarchical_All_Atom_PCA(outPCA, outFES, outKPCA);
								HAA_Evects_OR = Convoluted_Eigenvectors_AA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HAA_Cart();
									}

								oc.write_Hierarchical_PCA_Log(HAA + " INLIERS. ", number_of_residues_HAA_SS, number_of_atoms_hierarchical_AA_SS,
										Input_Parameters.MODES_EIGEN_RESIDUE_AA, Input_Parameters.MODES_HIERARCHICAL_AA, rank_HAA, trace_HAA, cond_HAA, det_HAA, shrinkage_HAA);

								outPCA = Input_Parameters.OUT_DIR + HAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Global_All_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_Outliers_SELECTED_AA);
								do_Hierarchical_All_Atom_PCA(outPCA, outFES, outKPCA);
								HAA_Evects_OS = Convoluted_Eigenvectors_AA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HAA_Cart();
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(HAA, "Inliers_vs_Outliers", HAA_Evects_OR, HAA_Evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + HAA + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_Hierarchical_PCA_Log(HAA + " OUTLIERS. ", number_of_residues_HAA_SS, number_of_atoms_hierarchical_AA_SS,
										Input_Parameters.MODES_EIGEN_RESIDUE_AA, Input_Parameters.MODES_HIERARCHICAL_AA, rank_HAA, trace_HAA, cond_HAA, det_HAA, shrinkage_HAA);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + HAA + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HAA + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HAA + File.separatorChar;
								do_Global_All_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_AA);
								do_Hierarchical_All_Atom_PCA(outPCA, outFES, outKPCA);
								HAA_Evects = Convoluted_Eigenvectors_AA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HAA_Cart();
									}
								oc.write_Hierarchical_PCA_Log(HAA, number_of_residues_HAA_SS, number_of_atoms_hierarchical_AA_SS, Input_Parameters.MODES_EIGEN_RESIDUE_AA,
										Input_Parameters.MODES_HIERARCHICAL_AA, rank_HAA, trace_HAA, cond_HAA, det_HAA, shrinkage_HAA);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed HIERARCHICAL ALL ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_hierarchical_HA)

			{
				future3 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outFES, outKPCA;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + HHA + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HHA + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HHA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Global_Heavy_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_Outliers_REMOVED_HA);
								do_Hierarchical_Heavy_Atom_PCA(outPCA, outFES, outKPCA);
								HHA_Evects_OR = Convoluted_Eigenvectors_HA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HHA_Cart();
									}

								oc.write_Hierarchical_PCA_Log(HHA + " INLIERS", number_of_residues_HHA_SS, number_of_atoms_hierarchical_HA_SS,
										Input_Parameters.MODES_EIGEN_RESIDUE_HA, Input_Parameters.MODES_HIERARCHICAL_HA, rank_HHA, trace_HHA, cond_HHA, det_HHA, shrinkage_HHA);

								outPCA = Input_Parameters.OUT_DIR + HHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Global_Heavy_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_Outliers_SELECTED_HA);
								do_Hierarchical_Heavy_Atom_PCA(outPCA, outFES, outKPCA);
								HHA_Evects_OS = Convoluted_Eigenvectors_HA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HHA_Cart();
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(HHA, "Inliers_vs_Outliers", HHA_Evects_OR, HHA_Evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + HHA + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_Hierarchical_PCA_Log(HHA + " OUTLIERS", number_of_residues_HHA_SS, number_of_atoms_hierarchical_HA_SS,
										Input_Parameters.MODES_EIGEN_RESIDUE_HA,

										Input_Parameters.MODES_HIERARCHICAL_HA, rank_HHA, trace_HHA, cond_HHA, det_HHA, shrinkage_HHA);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + HHA + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + HHA + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + HHA + File.separatorChar;
								do_Global_Heavy_Atom_Residue_PCA(outPCA, subset_PDB_coordinates_hierarchical_HA);
								do_Hierarchical_Heavy_Atom_PCA(outPCA, outFES, outKPCA);
								HHA_Evects = Convoluted_Eigenvectors_HA_Cart;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HHA_Cart();
									}
								oc.write_Hierarchical_PCA_Log(HHA, number_of_residues_HHA_SS, number_of_atoms_hierarchical_HA_SS, Input_Parameters.MODES_EIGEN_RESIDUE_HA,
										Input_Parameters.MODES_HIERARCHICAL_HA, rank_HHA, trace_HHA, cond_HHA, det_HHA, shrinkage_HHA);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed HIERARCHICAL HEAVY ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_cartesian_AA)
			{
				future4 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + cAA + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAA + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_All_Atom_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_REMOVED_AA);
								AA_Evects_OR = top_evectors_COV_AA;

								oc.write_PCA_Log(cAA + " INLIERS", number_of_residues_AA_SS, number_of_atoms_AA_SS, Input_Parameters.MODES_ALL_ATOM, rank_AA, trace_AA, cond_AA,
										det_AA, KMO_AA, shrinkage_AA);

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AA(cAA + File.separatorChar + "INLIERS");
									}

								outPCA = Input_Parameters.OUT_DIR + cAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_All_Atom_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_SELECTED_AA);
								AA_Evects_OS = top_evectors_COV_AA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AA(cAA + File.separatorChar + "OUTLIERS");
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cAA, "Inliers_vs_Outliers", AA_Evects_OR, AA_Evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_PCA_Log(cAA + " OUTLIERS", number_of_residues_AA_SS, number_of_atoms_AA_SS, Input_Parameters.MODES_ALL_ATOM, rank_AA, trace_AA, cond_AA,
										det_AA, KMO_AA, shrinkage_AA);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + cAA + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAA + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAA + File.separatorChar;
								do_All_Atom_PCA(outPCA, outSSA, outFES, outKPCA, aligned_subset_PDB_coordinates_AA);

								AA_Evects = top_evectors_COV_AA;
								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AA(cAA);
									}
								oc.write_PCA_Log(cAA, number_of_residues_AA_SS, number_of_atoms_AA_SS, Input_Parameters.MODES_ALL_ATOM, rank_AA, trace_AA, cond_AA, det_AA, KMO_AA,
										shrinkage_AA);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed ALL ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_cartesian_HA)
			{
				future5 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + cHA + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cHA + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cHA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Heavy_Atom_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_REMOVED_HA);
								HA_Evects_OR = top_evectors_COV_HA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HA(cHA + File.separatorChar + "INLIERS");
									}

								oc.write_PCA_Log(cHA + " INLIERS", number_of_residues_HA_SS, number_of_atoms_HA_SS, Input_Parameters.MODES_HEAVY_ATOM, rank_HA, trace_HA, cond_HA,
										det_HA, KMO_HA, shrinkage_HA);

								outPCA = Input_Parameters.OUT_DIR + cHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Heavy_Atom_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_SELECTED_HA);
								HA_Evects_OS = top_evectors_COV_HA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HA(cHA + File.separatorChar + "OUTLIERS");
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cHA, "Inliers_vs_Outliers", HA_Evects_OR, HA_Evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_PCA_Log(cHA + " OUTLIERS", number_of_residues_HA_SS, number_of_atoms_HA_SS, Input_Parameters.MODES_HEAVY_ATOM, rank_HA, trace_HA, cond_HA,
										det_HA, KMO_HA, shrinkage_HA);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + cHA + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cHA + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cHA + File.separatorChar;
								do_Heavy_Atom_PCA(outPCA, outSSA, outFES, outKPCA, aligned_subset_PDB_coordinates_HA);
								HA_Evects = top_evectors_COV_HA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_HA(cHA);
									}
							}
						oc.write_PCA_Log(cHA, number_of_residues_HA_SS, number_of_atoms_HA_SS, Input_Parameters.MODES_HEAVY_ATOM, rank_HA, trace_HA, cond_HA, det_HA, KMO_HA,
								shrinkage_HA);
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed HEAVY ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_hierarchical_AA && Input_Parameters.do_cartesian_AA)
			{
				try
					{
						if (future2 != null) future2.get();
						if (future4 != null) future4.get();
					}
				catch (Exception e)
					{
						e.printStackTrace();
					}
				String outSSA;
				System.out.println("Performing ALL ATOM Comparative Subspace Analysis: Hierarchical versus Direct: \n");
				try
					{
						if (Input_Parameters.doOutlierProcessing)
							{
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Comparative_AA_SSA(outSSA, HAA_Evects_OR, AA_Evects_OR);

								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Comparative_AA_SSA(outSSA, HAA_Evects_OS, AA_Evects_OS);
							}
						else
							{
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar;
								do_Comparative_AA_SSA(outSSA, HAA_Evects, AA_Evects);
							}
					}
				catch (Exception e)
					{
						e.printStackTrace();
					}
				System.out.println("----------------------------------------------------------------------------------------------\n");
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_hierarchical_HA && Input_Parameters.do_cartesian_HA)
			{
				try
					{
						if (future3 != null) future3.get();
						if (future5 != null) future5.get();
					}
				catch (Exception e)
					{
						e.printStackTrace();
					}
				String outSSA;
				System.out.println("Performing HEAVY ATOM Comparative Subspace Analysis: Hierarchical versus Direct \n");
				try
					{
						if (Input_Parameters.doOutlierProcessing)
							{
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Comparative_HA_SSA(outSSA, HHA_Evects_OR, HA_Evects_OR);

								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cHA + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Comparative_HA_SSA(outSSA, HHA_Evects_OS, HA_Evects_OS);
							}
						else
							{

								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAA + File.separatorChar;
								do_Comparative_HA_SSA(outSSA, HHA_Evects, HA_Evects);
							}
					}
				catch (Exception e)
					{
						e.printStackTrace();
					}
				System.out.println("----------------------------------------------------------------------------------------------\n");
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_cartesian_BB)
			{
				future6 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						Matrix evects_OR, evects_OS;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + cBB + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cBB + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cBB + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cBB + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Backbone_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_REMOVED_BB);
								evects_OR = top_evectors_COV_BB;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_BB(cBB + File.separatorChar + "INLIERS");
									}

								oc.write_PCA_Log(cBB + " INLIERS", number_of_residues_BB_SS, number_of_atoms_BB_SS, Input_Parameters.MODES_BACKBONE, rank_BB, trace_BB, cond_BB,
										det_BB, KMO_BB, shrinkage_BB);

								outPCA = Input_Parameters.OUT_DIR + cBB + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cBB + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cBB + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cBB + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Backbone_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_SELECTED_BB);
								evects_OS = top_evectors_COV_BB;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_BB(cBB + File.separatorChar + "OUTLIERS");
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cBB, "Inliers_vs_Outliers", evects_OR, evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cBB + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_PCA_Log(cBB + " OUTLIERS", number_of_residues_BB_SS, number_of_atoms_BB_SS, Input_Parameters.MODES_BACKBONE, rank_BB, trace_BB, cond_BB,
										det_BB, KMO_BB, shrinkage_BB);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + cBB + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cBB + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cBB + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cBB + File.separatorChar;
								do_Backbone_PCA(outPCA, outSSA, outFES, outKPCA, aligned_subset_PDB_coordinates_BB);

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_BB(cBB);
									}
								oc.write_PCA_Log(cBB, number_of_residues_BB_SS, number_of_atoms_BB_SS, Input_Parameters.MODES_BACKBONE, rank_BB, trace_BB, cond_BB, det_BB, KMO_BB,
										shrinkage_BB);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed BACKBONE ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_cartesian_CA)
			{
				future7 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						Matrix evects_OR, evects_OS;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + cAC + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAC + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAC + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAC + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Alpha_Carbon_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_REMOVED_CA);
								evects_OR = top_evectors_COV_CA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_CA(cAC + File.separatorChar + "INLIERS");
									}

								oc.write_PCA_Log(cAC + " INLIERS", number_of_residues_CA_SS, number_of_residues_CA_SS, Input_Parameters.MODES_ALPHA_CARBON, rank_CA, trace_CA,
										cond_CA, det_CA, KMO_CA, shrinkage_CA);

								outPCA = Input_Parameters.OUT_DIR + cAC + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAC + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAC + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAC + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Alpha_Carbon_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_SELECTED_CA);
								evects_OS = top_evectors_COV_CA;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_CA(cAC + File.separatorChar + "OUTLIERS");
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cAC, "Inliers_vs_Outliers", evects_OR, evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAC + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_PCA_Log(cAC + " OUTLIERS", number_of_residues_CA_SS, number_of_residues_CA_SS, Input_Parameters.MODES_ALPHA_CARBON, rank_CA, trace_CA,
										cond_CA, det_CA, KMO_CA, shrinkage_CA);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + cAC + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAC + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAC + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAC + File.separatorChar;
								do_Alpha_Carbon_PCA(outPCA, outSSA, outFES, outKPCA, aligned_subset_PDB_coordinates_CA);

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_CA(cAC);
									}
								oc.write_PCA_Log(cAC, number_of_residues_CA_SS, number_of_residues_CA_SS, Input_Parameters.MODES_ALPHA_CARBON, rank_CA, trace_CA, cond_CA, det_CA,
										KMO_CA, shrinkage_CA);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed ALPHA CARBON ATOM PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_atom_list)
			{
				future8 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						Matrix evects_OR, evects_OS;
						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + cAL + File.separatorChar + "INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAL + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAL + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAL + File.separatorChar + "INLIERS" + File.separatorChar;
								do_Atom_List_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_REMOVED_AL);
								evects_OR = top_evectors_COV_AL;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AL(cAL + File.separatorChar + "INLIERS");
									}

								oc.write_PCA_Log(cAL + " INLIERS", number_of_atoms_AL_SS, number_of_atoms_AL_SS, Input_Parameters.MODES_ATOMS_LIST, rank_AL, trace_AL, cond_AL,
										det_AL, KMO_AL, shrinkage_AL);

								outPCA = Input_Parameters.OUT_DIR + cAL + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAL + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAL + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAL + File.separatorChar + "OUTLIERS" + File.separatorChar;
								do_Atom_List_PCA(outPCA, outSSA, outFES, outKPCA, subset_PDB_coordinates_Outliers_SELECTED_AL);
								evects_OS = top_evectors_COV_AL;

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AL(cAL + File.separatorChar + "OUTLIERS");
									}

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(cAL, "Inliers_vs_Outliers", evects_OR, evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAL + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();

								oc.write_PCA_Log(cAL + " OUTLIERS", number_of_atoms_AL_SS, number_of_atoms_AL_SS, Input_Parameters.MODES_ATOMS_LIST, rank_AL, trace_AL, cond_AL,
										det_AL, KMO_AL, shrinkage_AL);
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + cAL + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + cAL + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + cAL + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + cAL + File.separatorChar;
								do_Atom_List_PCA(outPCA, outSSA, outFES, outKPCA, aligned_subset_PDB_coordinates_AL);

								if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz)
									{
										do_Mode_Visualization_AL(cAL);
									}

								oc.write_PCA_Log(cAL, number_of_atoms_AL_SS, number_of_atoms_AL_SS, Input_Parameters.MODES_ATOMS_LIST, rank_AL, trace_AL, cond_AL, det_AL, KMO_AL,
										shrinkage_AL);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed ATOM LIST PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		if (Input_Parameters.do_dist_pairs)
			{
				future9 = service.submit(new Runnable()
				{
					@Override
					public void run()
					{
						startTime = System.nanoTime();
						String outPCA, outSSA, outFES, outKPCA;
						Matrix evects_OR, evects_OS;

						if (Input_Parameters.doOutlierProcessing)
							{
								outPCA = Input_Parameters.OUT_DIR + DP + File.separatorChar + " INLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + DP + File.separatorChar + "INLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + DP + File.separatorChar + "INLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + DP + File.separatorChar + "INLIERS" + File.separatorChar;

								do_Distance_Pair_PCA(outPCA, outSSA, outFES, outKPCA, distances_Outliers_REMOVED);
								evects_OR = top_distance_evectors_COV;
								oc.write_DP_PCA_Log("INLIERS", number_of_atom_pairs_DP_SS, rank_d, trace_d, cond_d, det_d, KMO_DP, shrinkage_d, chain_idents1,
										atom_list_dp_original1, chain_idents2, atom_list_dp_original2, atomic_distance_means, atomic_distance_std_devs);

								outPCA = Input_Parameters.OUT_DIR + DP + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + DP + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + DP + File.separatorChar + "OUTLIERS" + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + DP + File.separatorChar + "OUTLIERS" + File.separatorChar;

								do_Distance_Pair_PCA(outPCA, outSSA, outFES, outKPCA, distances_Outliers_SELECTED);
								evects_OS = top_distance_evectors_COV;
								oc.write_DP_PCA_Log(" OUTLIERS", number_of_atom_pairs_DP_SS, rank_d, trace_d, cond_d, det_d, KMO_DP, shrinkage_d, chain_idents1,
										atom_list_dp_original1, chain_idents2, atom_list_dp_original2, atomic_distance_means, atomic_distance_std_devs);

								JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis(DP, "Inliers_vs_Outliers", evects_OR, evects_OS);
								ssa.setOut_dir(Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + DP + File.separatorChar);
								ssa.get_SSA();
								ssa.get_FSSA_Iterated();
							}
						else
							{
								outPCA = Input_Parameters.OUT_DIR + DP + File.separatorChar;
								outSSA = Input_Parameters.OUT_DIR + "SSA" + File.separatorChar + DP + File.separatorChar;
								outFES = Input_Parameters.OUT_DIR + "FES" + File.separatorChar + DP + File.separatorChar;
								outKPCA = Input_Parameters.OUT_DIR + "KPCA" + File.separatorChar + DP + File.separatorChar;

								do_Distance_Pair_PCA(outPCA, outSSA, outFES, outKPCA, distances);
								oc.write_DP_PCA_Log(" No Outlier Processing", number_of_atom_pairs_DP_SS, rank_d, trace_d, cond_d, det_d, KMO_DP, shrinkage_d, chain_idents1,
										atom_list_dp_original1, chain_idents2, atom_list_dp_original2, atomic_distance_means, atomic_distance_std_devs);
							}
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Performed DISTANCE PAIR PCA Analysis... \n");
						System.out.println("Done. (" + nf3.format(totalTime / 1.000E9) + " seconds)\n");
						System.out.println("----------------------------------------------------------------------------------------------\n");
					}
				});
			}
		/* ----------------------------------------------------------------------------------------------------------- */
		try
			{
				if (future0 != null) future0.get();
				if (future1 != null) future1.get();
				if (future2 != null) future2.get();
				if (future3 != null) future3.get();
				if (future4 != null) future4.get();
				if (future5 != null) future5.get();
				if (future6 != null) future6.get();
				if (future7 != null) future7.get();
				if (future8 != null) future8.get();
				if (future9 != null) future9.get();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}

		if (!Input_Parameters.doPREPROCESS) oc.write_DVP_Log();
		if (Input_Parameters.doFES) oc.write_FES_Log();
		if (Input_Parameters.doKPCA) oc.write_KPCA_Log();
		if (Input_Parameters.doCORR || Input_Parameters.doPCORR || doInterSSA_AA || doInterSSA_HA) oc.write_SSA_Log(doInterSSA_AA, doInterSSA_HA);
		if (Input_Parameters.doEssentialViz || Input_Parameters.doModeViz) oc.write_VIZ_Log();

		date = DateUtils.now();
		jediEnd = System.nanoTime();
		totalTime = jediEnd - jediStart;

		System.out.println("\n\nJEDi Analysis Completed: " + date);
		System.out.println("Total Compute Time = " + nf3.format(totalTime / (60.000E9)) + " minutes.");
		System.out.println("----------------------------------------------------------------------------------------------\n");

		oc.close_JEDi_Log(date, totalTime);
		service.shutdown();
	}
}
