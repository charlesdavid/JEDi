package jedi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

/**
 * JED class JEDi_Driver: Driver program for running JEDi. Input file is "JEDi_Parameters.txt" Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Driver
	{

		static boolean doREAD, doMULTI, do_cartesian_AA, do_cartesian_HA, do_cartesian_CA, do_residue, do_hierarchical, do_dist_pairs, doVIZ, do_no_pca, exist, success;
		static int NUMBER_OF_MODES_RESIDUE, NUMBER_OF_MODES_ALL_ATOM_CARTESIAN, NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN, NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN,
				NUMBER_OF_MODES_HIERARCHCAL, NUMBER_OF_MODES_DISTANCE_PAIRS, NUMBER_OF_MODES_VIZ, number_of_residues_REF, number_of_residues_AA_Cart_SS,
				number_of_residues_HA_Cart_SS, number_of_residues_CA_Cart_SS, number_of_residues_Hierarchical_SS, number_Of_Input_Lines, number_of_atoms_REF,
				number_of_heavy_atoms_REF, number_of_atoms_AA_Cart_SS, number_of_atoms_HA_Cart_SS, number_of_atom_pairs, number_conformations, ROWS, ROWS_CA, ROWS_HA, COLS;
		static double Z_SCORE_CUTOFF, VIZ_MODE_AMPLITUDE, trace_COV, trace_d_cov, trace_CORR, trace_d_corr, trace_d_pcorr, trace_PCORR, cond_cov, cond_d_cov, cond_corr,
				cond_d_corr, cond_pcorr, cond_d_pcorr, det_cov, det_d_cov, det_corr, det_d_corr, det_pcorr, det_d_pcorr, rank_cov, rank_d_cov, rank_corr, rank_d_corr, rank_pcorr,
				rank_d_pcorr;
		static String key, value, path, name, line, date, type, model, DIRECTORY, DESCRIPTION, REFERENCE_PDB, RESIDUE_LIST_ALL_ATOM, RESIDUE_LIST_HEAVY_ATOM,
				RESIDUE_LIST_ALPHA_CARBON, RESIDUE_LIST_HIERARCHICAL, ATOM_PAIRS_LIST, ORIGINAL_COORDS_ALL_ATOM, ORIGINAL_COORDS_HEAVY_ATOM, ORIGINAL_COORDS_ALPHA_CARBONS, OUT_DIR;
		static final String input_file = "JEDi_Parameters.txt", delim = "=", Q = "COV", R = "CORR", P = "PCORR", AA = "All_Atom_cPCA", AC = "Alpha_Carbon_cPCA",
				HA = "Heavy_Atom_cPCA", Res = "Residue_cPCA", RGC = "Hierarchical_PCA", DP = "dpPCA";
		static double[] pca_mode_min_COV, pca_mode_max_COV, pca_mode_min_CORR, pca_mode_max_CORR, pca_mode_min_PCORR, pca_mode_max_PCORR, residue_distance_means,
				residue_distance_std_devs;
		static List<Integer> residue_list_cartesian_AA, residue_list_cartesian_AA_orig, residue_list_cartesian_HA, residue_list_cartesian_HA_orig, residue_list_cartesian_CA,
				residue_list_cartesian_CA_orig, residue_list_hierarchical, residue_list_hierarchical_orig, atom_list_dp1, atom_list_dp2, atom_list_dp_original1,
				atom_list_dp_original2, atoms_read, numbers_Of_Atoms_in_Residues, numbers_Of_Heavy_Atoms_in_Residues, residues_read;
		static List<Double> original_conformation_rmsds, transformed_conformation_rmsds, transformed_residue_rmsds, top_cartesian_eigenvalues_COV, top_distance_eigenvalues_COV,
				top_cartesian_eigenvalues_CORR, top_distance_eigenvalues_CORR, top_cartesian_eigenvalues_PCORR, top_distance_eigenvalues_PCORR;
		static List<String> lines, pdb_file_names, chain_idents, chain_idents1, chain_idents2, chain_ids_read;
		static List<Matrix> Residue_Coordinates, Residue_Coordinates_HA, Reference_Residue_Coordinates, Residue_Delta_Vectors, Residue_Centered_Coordinates,
				Residue_Eigenvectors_COV;
		static List<List<Double>> Residue_Eigenvalues_COV, Residue_Eigenvalues_CORR, Residue_Eigenvalues_PCORR;
		static Matrix original_reference_PDB_coordinates_AA, original_reference_PDB_coordinates_CA, original_reference_PDB_coordinates_HA, original_PDB_coordinates_AA,
				original_PDB_coordinates_CA, original_PDB_coordinates_HA, subset_PDB_coordinates_cartesian_AA, subset_reference_PDB_coordinates_cartesian_CA,
				subset_reference_PDB_coordinates_cartesian_AA, subset_PDB_coordinates_cartesian_CA, subset_reference_PDB_coordinates_cartesian_HA,
				subset_PDB_coordinates_cartesian_HA, subset_reference_PDB_coordinates_hierarchical, subset_PDB_coordinates_hierarchical, tranasformed_reference_PDB_coordinates;
		static Matrix weighted_pca_modes_COV, square_pca_modes_COV, weighted_square_pca_modes_COV, transformed_PDB_coordinates_cartesian, distance_matrix,
				top_cartesian_evectors_COV, top_distance_evectors_COV, weighted_pca_modes_CORR, square_pca_modes_CORR, square_pca_modes_PCORR, weighted_square_pca_modes_CORR,
				top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, top_distance_evectors_CORR, top_distance_evectors_PCORR, pca_modes_CORR, pca_modes_COV;
		static Matrix projections_COV, normed_projections_COV, projections_CORR, normed_projections_CORR, weighted_normed_projections_COV, weighted_normed_projections_CORR,
				weighted_projections_CORR, weighted_projections_COV, projections_PCORR, normed_projections_PCORR;
		static Matrix projections_dist_COV, normed_projections_dist_COV, projections_dist_CORR, normed_projections_dist_CORR, weighted_normed_projections_dist_COV,
				weighted_normed_projections_dist_CORR, weighted_projections_dist_CORR, weighted_projections_dist_COV, projections_dist_PCORR, normed_projections_dist_PCORR;

		static Matrix Residue_Generalized_Coordinates_COV, RGC_Eigenvectors, RGC_DVPs, G, G_Square_Modes, U, hDVPs;
		static Vector<Atom> ref_atoms_AA, ref_atoms_CA, ref_atoms_HA, ref_subset_atoms_AA, ref_subset_atoms_CA, ref_subset_atoms_HA, ref_subset_atoms_RES_HIER;
		static File log;
		static BufferedReader input_reader;
		static BufferedWriter log_writer;
		static DateUtils now;
		static StringTokenizer sToken;
		static NumberFormat nf3, nff0, nf6, df;
		static RoundingMode rm;
		static long startTime, endTime, totalTime;
		static Hashtable<String, String> parameters;
		static JEDi_PDB_Processing refPDB;
		static JEDi_Driver jedi = new JEDi_Driver();

		// ************************* CONSTRUCTOR ***************************** //

		public JEDi_Driver()
			{
				rm = RoundingMode.HALF_UP;
				df = new DecimalFormat("0.###E0");
				df.setRoundingMode(rm);

				nff0 = NumberFormat.getInstance();
				nff0.setRoundingMode(rm);
				nff0.setMaximumFractionDigits(0);
				nff0.setMinimumFractionDigits(0);

				nf3 = NumberFormat.getInstance();
				nf3.setMaximumFractionDigits(3);
				nf3.setMinimumFractionDigits(3);
				nf3.setRoundingMode(rm);

				nf6 = NumberFormat.getInstance();
				nf6.setMaximumFractionDigits(6);
				nf6.setMinimumFractionDigits(6);
				nf6.setRoundingMode(rm);
			}

		// ***************************************************** PRIVATE METHODS ********************************************************************** //

		private static void read_input_file()
			{
				parameters = new Hashtable<>();
				number_Of_Input_Lines = 0;
				lines = new ArrayList<>();
				System.out.println("Reading the input file: " + input_file);
				System.out.println("Below are the parameters to be processed:");
				System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
				try
					{
						while ((line = input_reader.readLine()) != null && line.length() >= 1)
							{
								lines.add(line);
								System.out.println(line);
								sToken = new StringTokenizer(line);
								key = sToken.nextToken(delim);
								value = sToken.nextToken(delim);
								parameters.put(key, value);
								number_Of_Input_Lines++;
							}

						input_reader.close();
						System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
						System.out.println("The number of lines of parameters in the input file is: " + number_Of_Input_Lines + "\n");
					} catch (IOException e)
					{
						System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void assign_parameters()
			{
				if (parameters.get("doREAD").equals("true"))
					doREAD = true;
				if (parameters.get("doMULTI").equals("true"))
					doMULTI = true;
				DIRECTORY = parameters.get("DIRECTORY");
				DESCRIPTION = parameters.get("DESCRIPTION");
				REFERENCE_PDB = parameters.get("REFERENCE_PDB");
				Z_SCORE_CUTOFF = Double.valueOf(parameters.get("Z_SCORE_CUTOFF"));
				NUMBER_OF_MODES_HIERARCHCAL = Integer.valueOf(parameters.get("NUMBER_OF_MODES_HIERARCHCAL"));
				NUMBER_OF_MODES_ALL_ATOM_CARTESIAN = Integer.valueOf(parameters.get("NUMBER_OF_MODES_ALL_ATOM_CARTESIAN"));
				NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN = Integer.valueOf(parameters.get("NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN"));
				NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN = Integer.valueOf(parameters.get("NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN"));
				NUMBER_OF_MODES_DISTANCE_PAIRS = Integer.valueOf(parameters.get("NUMBER_OF_MODES_DISTANCE_PAIRS"));
				NUMBER_OF_MODES_VIZ = Integer.valueOf(parameters.get("NUMBER_OF_MODES_VIZ"));
				VIZ_MODE_AMPLITUDE = Double.valueOf(parameters.get("VIZ_MODE_AMPLITUDE"));
				RESIDUE_LIST_ALL_ATOM = parameters.get("RESIDUE_LIST_ALL_ATOM");
				RESIDUE_LIST_HEAVY_ATOM = parameters.get("RESIDUE_LIST_HEAVY_ATOM");
				RESIDUE_LIST_ALPHA_CARBON = parameters.get("RESIDUE_LIST_ALPHA_CARBON");
				RESIDUE_LIST_HIERARCHICAL = parameters.get("RESIDUE_LIST_HIERARCHICAL");
				ATOM_PAIRS_LIST = parameters.get("ATOM_PAIRS_LIST");
				ORIGINAL_COORDS_ALL_ATOM = parameters.get("ORIGINAL_COORDS_ALL_ATOM");
				ORIGINAL_COORDS_HEAVY_ATOM = parameters.get("ORIGINAL_COORDS_HEAVY_ATOM");
				ORIGINAL_COORDS_ALPHA_CARBONS = parameters.get("ORIGINAL_COORDS_ALPHA_CARBONS");

				NUMBER_OF_MODES_RESIDUE = NUMBER_OF_MODES_HIERARCHCAL;
				/* **************************************************************************************************************************************** */
				if (!(DIRECTORY.endsWith(File.separator)))
					{
						System.err.println("Expected the directory to end with " + File.separator + ", but got: " + DIRECTORY);
						System.err.println("Attempting to fix...");
						DIRECTORY = DIRECTORY + File.separator;
					}
				boolean dir = new File(DIRECTORY).isDirectory();
				if (!dir)
					{
						System.err.println("The entered directory is not a proper directory: " + DIRECTORY);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				/* **************************************************************************************************************************************** */
				boolean check = new File(DIRECTORY + REFERENCE_PDB).exists();
				if (!check)
					{
						System.err.println("The PDB Reference File can not be found: " + DIRECTORY + REFERENCE_PDB);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				OUT_DIR = DIRECTORY + "JEDi_RESULTS_" + DESCRIPTION + File.separatorChar;
				exist = new File(OUT_DIR).exists();
				if (!exist)
					{
						success = (new File(OUT_DIR)).mkdirs();
						if (!success)
							{
								System.err.println("Failed to create the output directory: " + OUT_DIR);
								System.exit(0);
							}
					}
				/* **************************************************************************************************************************************** */
				if (doREAD)
					{
						do_no_pca = true;
						do_cartesian_AA = false;
						do_cartesian_HA = false;
						do_cartesian_CA = false;
						do_residue = false;
						do_hierarchical = false;
						do_dist_pairs = false;
						doVIZ = false;
					}
				/* **************************************************************************************************************************************** */
				if (NUMBER_OF_MODES_ALL_ATOM_CARTESIAN > 0)
					do_cartesian_AA = true;
				if (NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN > 0)
					do_cartesian_HA = true;
				if (NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN > 0)
					do_cartesian_CA = true;
				if (NUMBER_OF_MODES_RESIDUE > 0)
					do_residue = true;
				if (NUMBER_OF_MODES_HIERARCHCAL > 0)
					do_hierarchical = true;
				if (NUMBER_OF_MODES_DISTANCE_PAIRS > 0)
					do_dist_pairs = true;
				if (NUMBER_OF_MODES_VIZ > 0)
					doVIZ = true;

				/* **************************************************************************************************************************************** */
				if (do_cartesian_AA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_ALL_ATOM).exists();
						if (!check)
							{
								System.err.println("The Residue List File can not be found: " + DIRECTORY + RESIDUE_LIST_ALL_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* **************************************************************************************************************************************** */
				if (do_cartesian_HA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_HEAVY_ATOM).exists();
						if (!check)
							{
								System.err.println("The Residue List File can not be found: " + DIRECTORY + RESIDUE_LIST_HEAVY_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* **************************************************************************************************************************************** */
				if (do_cartesian_CA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_ALPHA_CARBON).exists();
						if (!check)
							{
								System.err.println("The Residue List File can not be found: " + DIRECTORY + RESIDUE_LIST_ALPHA_CARBON);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* **************************************************************************************************************************************** */
				if (do_hierarchical)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_HIERARCHICAL).exists();
						if (!check)
							{
								System.err.println("The Residue List File can not be found: " + DIRECTORY + RESIDUE_LIST_HIERARCHICAL);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* **************************************************************************************************************************************** */
				if (do_dist_pairs)
					{
						check = new File(DIRECTORY + ATOM_PAIRS_LIST).exists();
						if (!check)
							{
								System.err.println("The Atom Pairs List File can not be found: " + DIRECTORY + ATOM_PAIRS_LIST);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* ************************************************************************************************************************************ */
				if (!doREAD)
					{
						check = new File(DIRECTORY + ORIGINAL_COORDS_ALL_ATOM).exists();
						if (!check)
							{
								System.err.println("The All-Atom Coordinates Matrix File can not be found: " + DIRECTORY + ORIGINAL_COORDS_ALL_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
						check = new File(DIRECTORY + ORIGINAL_COORDS_HEAVY_ATOM).exists();
						if (!check)
							{
								System.err.println("The NCO Coordinates Matrix File can not be found: " + DIRECTORY + ORIGINAL_COORDS_HEAVY_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
						check = new File(DIRECTORY + ORIGINAL_COORDS_ALPHA_CARBONS).exists();
						if (!check)
							{
								System.err.println("The CA Coordinates Matrix File can not be found: " + DIRECTORY + ORIGINAL_COORDS_ALPHA_CARBONS);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				/* ************************************************************************************************************************************* */
				System.gc();
			}

		private static void initialize_JED_Log()
			{
				try
					{
						log = new File(OUT_DIR + "JEDi_LOG.txt");
						log_writer = new BufferedWriter(new FileWriter(log));
						log_writer.write("JEDi: Java Essential Dynamics Inspector" + "\n");
						log_writer.write("Job Parameters:" + "\n");
						for (String line : lines)
							log_writer.write("\t" + line + "\n");
						log_writer.flush();
					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void process_Reference_PDB()
			{
				refPDB = new JEDi_PDB_Processing(DIRECTORY, REFERENCE_PDB);
				refPDB.read_Reference_PDB();

				ref_atoms_AA = refPDB.get_ref_Atoms();
				ref_atoms_CA = refPDB.get_ref_Atoms_CA();
				ref_atoms_HA = refPDB.get_ref_Atoms_HA();

				number_of_atoms_REF = refPDB.getNumber_of_atoms();
				number_of_residues_REF = refPDB.getNumber_of_residues();
				number_of_heavy_atoms_REF = ref_atoms_HA.size();

				original_reference_PDB_coordinates_AA = refPDB.getReference_PDB_coordinates();
				original_reference_PDB_coordinates_CA = refPDB.getReference_PDB_coordinates_CA();
				original_reference_PDB_coordinates_HA = refPDB.getReference_PDB_coordinates_HA();

				Reference_Residue_Coordinates = refPDB.getReference_Residue_Coords_List();

				atoms_read = refPDB.getAtoms_read();
				residues_read = refPDB.getResidues_read();
				numbers_Of_Atoms_in_Residues = refPDB.getNumber_of_Atoms_in_Residues();
				numbers_Of_Heavy_Atoms_in_Residues = refPDB.getNumber_of_Heavy_Atoms_in_Residues();

				ROWS = original_reference_PDB_coordinates_AA.getRowDimension();
				ROWS_CA = original_reference_PDB_coordinates_CA.getRowDimension();
				ROWS_HA = original_reference_PDB_coordinates_HA.getRowDimension();

				name = "original_Reference_PDB_coordinates_AA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_reference_PDB_coordinates_AA, path, 9, 3);

				name = "original_Reference_PDB_coordinates_CA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_reference_PDB_coordinates_CA, path, 9, 3);

				name = "original_Reference_PDB_coordinates_HA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_reference_PDB_coordinates_HA, path, 9, 3);

				name = "numbers_Of_Atoms_in_Residues.txt";
				path = OUT_DIR + name;
				List_IO.write_Integer_List(numbers_Of_Atoms_in_Residues, path);

				name = "numbers_Of_Heavy_Atoms_in_Residues.txt";
				path = OUT_DIR + name;
				List_IO.write_Integer_List(numbers_Of_Heavy_Atoms_in_Residues, path);

				if (!doMULTI)
					{
						name = "All_PDB_Residues_JED.txt";
						path = OUT_DIR + name;
						List_IO.write_Integer_List(residues_read, path);
						write_reference_PDB_Log();
					}

				if (doMULTI)
					{
						chain_ids_read = refPDB.getChain_ids_read();
						name = "All_PDB_Residues_Multi_JED.txt";
						path = OUT_DIR + name;
						List_IO.write_String_Integer_List(chain_ids_read, residues_read, path);
						write_reference_PDB_Log();
					}
			}

		private static void write_reference_PDB_Log()
			{
				try
					{
						log_writer.write("\nReading the Reference PDB file: " + REFERENCE_PDB + "\n");
						log_writer.write("\nReading all atom coordinates from the reference PDB file." + "\n");
						log_writer.write("The number of residues found in the Reference PDB file = " + number_of_residues_REF + "\n");
						log_writer.write("The number of atoms found in the Reference PDB file = " + number_of_atoms_REF + "\n");
						log_writer.write("The number of Heavy atoms found in the Reference PDB file = " + number_of_heavy_atoms_REF + "\n");
						log_writer.write("Reference All-Atom (AA) coordinates matrix created: 'original_Reference_PDB_Coordinates_AA'" + "\n");
						log_writer.write("Reference Heavy-Atom (HA) coordinates matrix created: 'original_Reference_PDB_Coordinates_HA'" + "\n");
						log_writer.write("Reference Alpha-Carbon (CA) coordinates matrix created: 'original_Reference_PDB_Coordinates_CA'" + "\n");
						if (!doMULTI)
							log_writer.write("The file " + name + "' contains all residue numbers found in the PDB file, formatted for Single Chain Analysis Input." + "\n");
						if (doMULTI)
							log_writer.write(
									"The file " + name + "' contains all chainID-residue number pairs found in the PDB file, formatted for Multi Chain Analysis Input." + "\n");
						log_writer.write("Number of atoms per residue list created: 'numbers_Of_Atoms_in_Residues'" + "\n");
						log_writer.write("Number of heavy atoms per residue list created: 'numbers_Of_Heavy_Atoms_in_Residues'" + "\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JED_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void read_PDB_Files()
			{
				original_PDB_coordinates_AA = refPDB.get_Original_PDB_Coords();
				original_PDB_coordinates_HA = refPDB.getOriginal_PDB_coordinates_HA();
				original_PDB_coordinates_CA = refPDB.getOriginal_PDB_coordinates_CA();

				COLS = original_PDB_coordinates_AA.getColumnDimension();
				pdb_file_names = refPDB.get_PDB_file_names();
				number_conformations = COLS;

				name = "original_PDB_Coordinates_AA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_PDB_coordinates_AA, path, 9, 3);

				name = "original_PDB_Coordinates_HA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_PDB_coordinates_HA, path, 9, 3);

				name = "original_PDB_Coordinates_CA.txt";
				path = OUT_DIR + name;
				Matrix_IO.write_Matrix(original_PDB_coordinates_CA, path, 9, 3);

				name = "PDB_Read_Log.txt";
				path = OUT_DIR + name;
				List_IO.write_String_List(pdb_file_names, path);

				if (!doMULTI)
					{
						name = "All_PDB_Residues_JED.txt";
						path = OUT_DIR + name;
						List_IO.write_Integer_List(residues_read, path);
					}
				if (doMULTI)
					{
						chain_ids_read = refPDB.getChain_ids_read();
						name = "All_PDB_Residues_Multi_JED.txt";
						path = OUT_DIR + name;
						List_IO.write_String_Integer_List(chain_ids_read, residues_read, path);
					}

				write_Read_Log();
				System.gc();
			}

		private static void write_Read_Log()
			{
				try
					{
						log_writer.write("\nPerforming Preliminary Processing Run.\n");
						log_writer.write("\nReading all atom coordinates from all PDB files in the Working Directory" + "\n");
						log_writer.write("The number of PDB files read = " + number_conformations + "\n");
						log_writer.write("The number of residues found in the PDB files = " + number_of_residues_REF + "\n");
						log_writer.write("The number of atoms found in the PDB files = " + number_of_atoms_REF + "\n");
						log_writer.write("All Atom Coordinates Matrix created: 'original_PDB_Coordinates'" + "\n");
						log_writer.write("\tThe dimension of this matrix is: " + number_of_atoms_REF * 3 + "  by  " + number_conformations + "\n");
						log_writer.write("PDB reference structure is: " + REFERENCE_PDB + "\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void read_All_Atom_Coordinate_File()
			{

				JEDi_Get_Coordinates_from_Matrix mat_coords = new JEDi_Get_Coordinates_from_Matrix(DIRECTORY, ORIGINAL_COORDS_ALL_ATOM);
					{
						original_PDB_coordinates_AA = mat_coords.get_Original_PDB_coordinates();
						number_conformations = original_PDB_coordinates_AA.getColumnDimension();
					}

					{
						try
							{
								log_writer.write("\nThe all atom coordinates were obtained from coordinates matrix file: " + ORIGINAL_COORDS_ALL_ATOM + "\n");
								log_writer.write("The dimension of the coordinates matrix is = " + ROWS + " by " + number_conformations + "\n");
								log_writer.write("Total number of atoms in matrix = " + (ROWS / 3) + "\n");
								log_writer.write("Total number of conformations in matrix = " + number_conformations + "\n");
								log_writer.write("Transformed PDB coordinates obtained by quaternion alignment to the reference structure." + "\n");
								log_writer.write("PDB reference structure is: " + REFERENCE_PDB + "\n");
								log_writer.flush();

							} catch (IOException e)
							{
								System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
								System.err.println("Program terminating.\n");
								e.printStackTrace();
								System.exit(0);
							}
						System.gc();
					}
			}

		private static void read_HA_Coordinate_File()
			{

				JEDi_Get_Coordinates_from_Matrix mat_coords = new JEDi_Get_Coordinates_from_Matrix(DIRECTORY, ORIGINAL_COORDS_HEAVY_ATOM);
					{
						original_PDB_coordinates_HA = mat_coords.get_Original_PDB_coordinates();
						number_conformations = original_PDB_coordinates_HA.getColumnDimension();
					}

					{
						try
							{
								log_writer.write("\nThe Heavy Atom (HA) coordinates were obtained from coordinates matrix file: " + ORIGINAL_COORDS_HEAVY_ATOM + "\n");
								log_writer.write("The dimension of the coordinates matrix is = " + ROWS_HA + " by " + number_conformations + "\n");
								log_writer.write("Total number of heavy atoms in matrix = " + (ROWS_HA / 3) + "\n");
								log_writer.write("Total number of conformations in matrix = " + number_conformations + "\n");
								log_writer.write("Transformed PDB coordinates obtained by quaternion alignment to the reference structure." + "\n");
								log_writer.write("PDB reference structure is: " + REFERENCE_PDB + "\n");
								log_writer.flush();

							} catch (IOException e)
							{
								System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
								System.err.println("Program terminating.\n");
								e.printStackTrace();
								System.exit(0);
							}
						System.gc();
					}
			}

		private static void read_CA_Coordinate_File()
			{

				JEDi_Get_Coordinates_from_Matrix mat_coords = new JEDi_Get_Coordinates_from_Matrix(DIRECTORY, ORIGINAL_COORDS_ALPHA_CARBONS);
					{
						original_PDB_coordinates_CA = mat_coords.get_Original_PDB_coordinates();
						number_conformations = original_PDB_coordinates_CA.getColumnDimension();
					}

					{
						try
							{
								log_writer.write("\nThe alpha-carbon atom coordinates were obtained from coordinates matrix file: " + ORIGINAL_COORDS_ALPHA_CARBONS + "\n");
								log_writer.write("The dimension of the coordinates matrix is = " + ROWS_CA + " by " + number_conformations + "\n");
								log_writer.write("Total number of atoms in the matrix = " + number_of_residues_REF + "\n");
								log_writer.write("Total number of conformations in matrix = " + number_conformations + "\n");
								log_writer.write("Transformed PDB coordinates obtained by quaternion alignment to the reference structure." + "\n");
								log_writer.write("PDB reference structure is: " + REFERENCE_PDB + "\n");
								log_writer.flush();

							} catch (IOException e)
							{
								System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
								System.err.println("Program terminating.\n");
								e.printStackTrace();
								System.exit(0);
							}
						System.gc();
					}
			}

		private static void get_Subsets()
			{
				try
					{
						if (do_cartesian_AA)
							{
								residue_list_cartesian_AA = refPDB.read_residue_list_Single(RESIDUE_LIST_ALL_ATOM);
								residue_list_cartesian_AA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_AA = refPDB.get_Reference_Subset_Single(ref_atoms_AA, residue_list_cartesian_AA_orig);
								number_of_atoms_AA_Cart_SS = ref_subset_atoms_AA.size();
								number_of_residues_AA_Cart_SS = residue_list_cartesian_AA.size();
								subset_reference_PDB_coordinates_cartesian_AA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_AA);
								subset_PDB_coordinates_cartesian_AA = refPDB.get_subset_Cartesian_Coords_All_Atom(Residue_Coordinates, residue_list_cartesian_AA,
										numbers_Of_Atoms_in_Residues);

								log_writer.write("\nThe Cartesian All-Atom Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_ALL_ATOM + "\n");
								log_writer.write("\tThe number of residues in the All-Atom Subset = " + number_of_residues_AA_Cart_SS + "\n");
								log_writer.write("\tThe number of atoms in the All-Atom Subset = " + number_of_atoms_AA_Cart_SS + "\n");
								log_writer.write(
										"The dimension of the All-Atom Reference Subset coordinates matrix is = " + subset_reference_PDB_coordinates_cartesian_AA.getRowDimension()
												+ " by " + subset_reference_PDB_coordinates_cartesian_AA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the All-Atom Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_AA.getRowDimension() + " by "
										+ subset_PDB_coordinates_cartesian_AA.getColumnDimension() + "\n");
								log_writer.write("\tThe number of atoms per residue in the AA subset are: \n");
								for (int i = 0; i < residue_list_cartesian_AA_orig.size(); i++)
									{
										int res = residue_list_cartesian_AA_orig.get(i);
										int num = numbers_Of_Atoms_in_Residues.get(res - 1);
										log_writer.write("\tResidue: " + res + "\t" + "Number of atoms: " + num + "\n");
									}
								log_writer.flush();
							}

						if (do_cartesian_HA)
							{
								residue_list_cartesian_HA = refPDB.read_residue_list_Single(RESIDUE_LIST_HEAVY_ATOM);
								residue_list_cartesian_HA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_HA = refPDB.get_Reference_Subset_Single(ref_atoms_HA, residue_list_cartesian_HA_orig);
								number_of_atoms_HA_Cart_SS = ref_subset_atoms_HA.size();
								number_of_residues_HA_Cart_SS = residue_list_cartesian_HA.size();
								subset_reference_PDB_coordinates_cartesian_HA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_HA);
								subset_PDB_coordinates_cartesian_HA = refPDB.get_subset_Cartesian_Coords_HA(Residue_Coordinates_HA, residue_list_cartesian_HA,
										numbers_Of_Heavy_Atoms_in_Residues);

								log_writer.write("\nThe Cartesian Heavy-Atom Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_HEAVY_ATOM + "\n");
								log_writer.write("\tThe number of residues in the Heavy Atom (HA) Subset = " + number_of_residues_HA_Cart_SS + "\n");
								log_writer.write("\tThe number of atoms in the Heavy Atom (HA) Subset = " + number_of_atoms_HA_Cart_SS + "\n");
								log_writer
										.write("The dimension of the HA Reference Subset coordinates matrix is = " + subset_reference_PDB_coordinates_cartesian_HA.getRowDimension()
												+ " by " + subset_reference_PDB_coordinates_cartesian_HA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the HA Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_HA.getRowDimension() + " by "
										+ subset_PDB_coordinates_cartesian_HA.getColumnDimension() + "\n");
								log_writer.write("\tThe number of heavy-atoms per residue in the HA subset are: \n");
								for (int i = 0; i < residue_list_cartesian_HA_orig.size(); i++)
									{
										int res = residue_list_cartesian_HA_orig.get(i);
										int num_heavy = numbers_Of_Heavy_Atoms_in_Residues.get(res - 1);
										log_writer.write("\tResidue: " + res + "\t" + "  Number of heavy atoms: " + num_heavy + "\n");
									}
								log_writer.flush();
							}

						if (do_cartesian_CA)
							{
								residue_list_cartesian_CA = refPDB.read_residue_list_Single(RESIDUE_LIST_ALPHA_CARBON);
								residue_list_cartesian_CA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_CA = refPDB.get_Reference_Subset_Single(ref_atoms_CA, residue_list_cartesian_CA_orig);
								number_of_residues_CA_Cart_SS = ref_subset_atoms_CA.size();
								subset_reference_PDB_coordinates_cartesian_CA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_CA);
								subset_PDB_coordinates_cartesian_CA = refPDB.get_subset_Cartesian_Coords_CA(original_PDB_coordinates_CA, residue_list_cartesian_CA);

								log_writer.write("\nThe Cartesian Alpha-Carbon Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_ALPHA_CARBON + "\n");
								log_writer.write("\tThe number of residues in the Alpha-Carbon Subset = " + number_of_residues_CA_Cart_SS + "\n");
								log_writer.write("The dimension of the Alpha-Carbon Reference Subset coordinates matrix is = "
										+ subset_reference_PDB_coordinates_cartesian_CA.getRowDimension() + " by "
										+ subset_reference_PDB_coordinates_cartesian_CA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the Alpha-Carbon Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_CA.getRowDimension()
										+ " by " + subset_PDB_coordinates_cartesian_CA.getColumnDimension() + "\n");
								log_writer.flush();
							}

						if (do_hierarchical)
							{
								residue_list_hierarchical = refPDB.read_residue_list_Single(RESIDUE_LIST_HIERARCHICAL);
								residue_list_hierarchical_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_RES_HIER = refPDB.get_Reference_Subset_Single(ref_atoms_AA, residue_list_hierarchical_orig);
								number_of_atoms_AA_Cart_SS = ref_subset_atoms_RES_HIER.size();
								number_of_residues_Hierarchical_SS = residue_list_hierarchical.size();
								subset_reference_PDB_coordinates_hierarchical = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_RES_HIER);
								subset_PDB_coordinates_hierarchical = refPDB.get_subset_Cartesian_Coords_All_Atom(Residue_Coordinates, residue_list_hierarchical,
										numbers_Of_Atoms_in_Residues);

								log_writer.write("\nThe Residue-Hierarchical Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_HIERARCHICAL + "\n");
								log_writer.write("\tThe number of residues in the Residue-Hierarchical Subset = " + number_of_residues_Hierarchical_SS + "\n");
								log_writer.write("The dimension of the Residue-Hierarchical Reference Subset coordinates matrix is = "
										+ subset_reference_PDB_coordinates_hierarchical.getRowDimension() + " by "
										+ subset_reference_PDB_coordinates_hierarchical.getColumnDimension() + "\n");
								log_writer.write("The dimension of the Residue-Hierarchical Subset coordinates matrix is = " + subset_PDB_coordinates_hierarchical.getRowDimension()
										+ " by " + subset_PDB_coordinates_hierarchical.getColumnDimension() + "\n");
								log_writer.flush();
							}

						if (do_dist_pairs)
							{
								refPDB.read_atom_pairs_Single(ATOM_PAIRS_LIST);
								atom_list_dp1 = refPDB.getAtom_list1();
								atom_list_dp2 = refPDB.getAtom_list2();
								atom_list_dp_original1 = refPDB.getAtom_list1_original();
								atom_list_dp_original2 = refPDB.getAtom_list2_original();
								number_of_atom_pairs = atom_list_dp1.size();
								log_writer.write("\nThe Distance-Pair Subset coordinates were obtained using this residue list: " + ATOM_PAIRS_LIST + "\n");
								log_writer.write("\tThe number of atoms in the Distance-Pair Subset = " + number_of_atom_pairs + "\n");
								log_writer.flush();
							}
					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JED_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void get_Subsets_Multi()
			{
				try
					{
						if (do_cartesian_AA)
							{
								residue_list_cartesian_AA = refPDB.read_residue_list_Single(RESIDUE_LIST_ALL_ATOM);
								residue_list_cartesian_AA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_AA = refPDB.get_Reference_Subset_Single(ref_atoms_AA, residue_list_cartesian_AA_orig);
								number_of_atoms_AA_Cart_SS = ref_subset_atoms_AA.size();
								number_of_residues_AA_Cart_SS = residue_list_cartesian_AA.size();
								subset_reference_PDB_coordinates_cartesian_AA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_RES_HIER);
								subset_PDB_coordinates_cartesian_AA = refPDB.get_subset_Cartesian_Coords_All_Atom(Residue_Coordinates, residue_list_cartesian_AA,
										numbers_Of_Atoms_in_Residues);
								log_writer.write("\nThe Cartesian All-Atom Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_ALL_ATOM + "\n");
								log_writer.write("\tThe number of residues in the All-Atom Subset = " + number_of_residues_AA_Cart_SS + "\n");
								log_writer.write("\tThe number of atoms in the All-Atom Subset = " + number_of_atoms_AA_Cart_SS + "\n");
								log_writer.write(
										"The dimension of the All-Atom Reference Subset coordinates matrix is = " + subset_reference_PDB_coordinates_cartesian_AA.getRowDimension()
												+ " by " + subset_reference_PDB_coordinates_cartesian_AA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the All-Atom Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_AA.getRowDimension() + " by "
										+ subset_PDB_coordinates_cartesian_AA.getColumnDimension() + "\n");
								log_writer.write("\tThe number of atoms per residue in the AA subset are: \n");
								for (int i = 0; i < residue_list_cartesian_AA_orig.size(); i++)
									{
										int res = residue_list_cartesian_AA_orig.get(i);
										int num = numbers_Of_Atoms_in_Residues.get(res - 1);
										log_writer.write("\tResidue: " + res + "\t" + "Number of atoms: " + num + "\n");
									}
								log_writer.flush();
							}

						if (do_cartesian_HA)
							{
								residue_list_cartesian_HA = refPDB.read_residue_list_Single(RESIDUE_LIST_HEAVY_ATOM);
								residue_list_cartesian_HA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_HA = refPDB.get_Reference_Subset_Single(ref_atoms_HA, residue_list_cartesian_HA_orig);
								number_of_atoms_HA_Cart_SS = ref_subset_atoms_HA.size();
								number_of_residues_HA_Cart_SS = residue_list_cartesian_HA.size();
								subset_reference_PDB_coordinates_cartesian_HA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_HA);
								subset_PDB_coordinates_cartesian_HA = refPDB.get_subset_Cartesian_Coords_HA(Residue_Coordinates_HA, residue_list_cartesian_HA,
										numbers_Of_Heavy_Atoms_in_Residues);
								log_writer.write("\nThe Cartesian Heavy-Atom Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_HEAVY_ATOM + "\n");
								log_writer.write("\tThe number of residues in the Heavy Atom (HA) Subset = " + number_of_residues_HA_Cart_SS + "\n");
								log_writer.write("\tThe number of atoms in the Heavy Atom (HA) Subset = " + number_of_atoms_HA_Cart_SS + "\n");
								log_writer
										.write("The dimension of the HA Reference Subset coordinates matrix is = " + subset_reference_PDB_coordinates_cartesian_HA.getRowDimension()
												+ " by " + subset_reference_PDB_coordinates_cartesian_HA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the HA Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_HA.getRowDimension() + " by "
										+ subset_PDB_coordinates_cartesian_HA.getColumnDimension() + "\n");
								log_writer.write("\tThe number of heavy-atoms per residue in the HA subset are: \n");
								for (int i = 0; i < residue_list_cartesian_HA_orig.size(); i++)
									{
										int res = residue_list_cartesian_HA_orig.get(i);
										int num_heavy = numbers_Of_Heavy_Atoms_in_Residues.get(res - 1);
										log_writer.write("\tResidue: " + res + "\t" + "  Number of heavy atoms: " + num_heavy + "\n");
									}
								log_writer.flush();
							}

						if (do_cartesian_CA)
							{
								residue_list_cartesian_CA = refPDB.read_residue_list_Single(RESIDUE_LIST_ALPHA_CARBON);
								residue_list_cartesian_CA_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_CA = refPDB.get_Reference_Subset_Single(ref_atoms_CA, residue_list_cartesian_CA_orig);
								number_of_residues_CA_Cart_SS = ref_subset_atoms_CA.size();
								subset_reference_PDB_coordinates_cartesian_CA = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_CA);
								subset_PDB_coordinates_cartesian_CA = refPDB.get_subset_Cartesian_Coords_CA(original_PDB_coordinates_CA, residue_list_cartesian_AA);
								log_writer.write("\nThe Cartesian Alpha-Carbon Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_ALPHA_CARBON + "\n");
								log_writer.write("\tThe number of residues in the Alpha-Carbon Subset = " + number_of_residues_CA_Cart_SS + "\n");
								log_writer.write("The dimension of the Alpha-Carbon Reference Subset coordinates matrix is = "
										+ subset_reference_PDB_coordinates_cartesian_CA.getRowDimension() + " by "
										+ subset_reference_PDB_coordinates_cartesian_CA.getColumnDimension() + "\n");
								log_writer.write("The dimension of the Alpha-Carbon Subset coordinates matrix is = " + subset_PDB_coordinates_cartesian_CA.getRowDimension()
										+ " by " + subset_PDB_coordinates_cartesian_CA.getColumnDimension() + "\n");
								log_writer.flush();
							}

						if (do_hierarchical)
							{
								residue_list_hierarchical = refPDB.read_residue_list_Single(RESIDUE_LIST_HIERARCHICAL);
								residue_list_hierarchical_orig = refPDB.getResidue_list_original();
								ref_subset_atoms_RES_HIER = refPDB.get_Reference_Subset_Single(ref_atoms_AA, residue_list_hierarchical_orig);
								number_of_atoms_AA_Cart_SS = ref_subset_atoms_RES_HIER.size();
								number_of_residues_Hierarchical_SS = residue_list_hierarchical.size();
								subset_reference_PDB_coordinates_hierarchical = refPDB.get_Cartesian_PDB_Coords(ref_subset_atoms_RES_HIER);
								subset_PDB_coordinates_hierarchical = refPDB.get_subset_Cartesian_Coords_All_Atom(Residue_Coordinates, residue_list_hierarchical,
										numbers_Of_Atoms_in_Residues);
								log_writer.write("\nThe Residue-Hierarchical Subset coordinates were obtained using this residue list: " + RESIDUE_LIST_HIERARCHICAL + "\n");
								log_writer.write("\tThe number of residues in the Residue-Hierarchical Subset = " + number_of_residues_Hierarchical_SS + "\n");
								log_writer.write("The dimension of the Residue-Hierarchical Reference Subset coordinates matrix is = "
										+ subset_reference_PDB_coordinates_hierarchical.getRowDimension() + " by "
										+ subset_reference_PDB_coordinates_hierarchical.getColumnDimension() + "\n");
								log_writer.write("The dimension of the Residue-Hierarchical Subset coordinates matrix is = " + subset_PDB_coordinates_hierarchical.getRowDimension()
										+ " by " + subset_PDB_coordinates_hierarchical.getColumnDimension() + "\n");
								log_writer.flush();
							}

						if (do_dist_pairs)
							{
								refPDB.read_atom_pairs_Multi(ATOM_PAIRS_LIST);
								atom_list_dp1 = refPDB.getAtom_list1();
								atom_list_dp2 = refPDB.getAtom_list2();
								atom_list_dp_original1 = refPDB.getAtom_list1_original();
								atom_list_dp_original2 = refPDB.getAtom_list2_original();
								chain_idents1 = refPDB.getChain_ids1();
								chain_idents2 = refPDB.getChain_ids2();
								number_of_atom_pairs = atom_list_dp1.size();
							}
					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JED_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void do_Residue_cPCA()
			{

				JEDi_Do_Residue_PCA cRES = new JEDi_Do_Residue_PCA(DIRECTORY, DESCRIPTION, residue_list_hierarchical, numbers_Of_Atoms_in_Residues, NUMBER_OF_MODES_RESIDUE,
						subset_PDB_coordinates_hierarchical, subset_reference_PDB_coordinates_hierarchical);

				cRES.set_z_cutoff(Z_SCORE_CUTOFF);
				// cRES.do_Cartesian_Residue_Local();
				cRES.do_Cartesian_Residue_Global();

				Residue_Generalized_Coordinates_COV = cRES.get_Residue_Generalized_Coordinates_COV();
				Matrix_IO.write_Matrix(Residue_Generalized_Coordinates_COV, OUT_DIR + "Residue_Generalized_Coordinates_Matrix.txt", 12, 9);

				Residue_Eigenvalues_COV = cRES.get_Residue_Eigenvalues_COV();
				Residue_Eigenvectors_COV = cRES.get_Residue_Eigenvectors_COV();
				Residue_Delta_Vectors = cRES.get_Residue_Delta_Vectors();
				Residue_Centered_Coordinates = cRES.getResidues_Centered_Data();

				path = OUT_DIR + "Residue_cPCA" + File.separatorChar;
				int i = 0;
				for (Matrix m : Residue_Eigenvectors_COV)
					{
						int res = residue_list_hierarchical_orig.get(i);
						String res_index = String.format("%03d", res);

						Matrix rdv = Residue_Delta_Vectors.get(i);
						Matrix rcc = Residue_Centered_Coordinates.get(i);
						List<Double> evals = Residue_Eigenvalues_COV.get(i);

						List_IO.write_Double_List(evals, path + "Residue_" + res_index + "_Top_" + NUMBER_OF_MODES_RESIDUE + "_Eigenvalues_COV.txt", 6);
						Matrix_IO.write_Matrix(rdv, path + "Residue_" + res_index + "_Delta_Vectors.txt", 12, 9);
						Matrix_IO.write_Matrix(m, path + "Residue_" + res_index + "_Top_" + NUMBER_OF_MODES_RESIDUE + "_Eigenvectors_COV.txt", 12, 6);
						Matrix_IO.write_Matrix(rcc, path + "Residue_" + res_index + "_Centered_Coordinates.txt", 12, 3);

						i++;
					}
				System.gc();
				write_Residue_cPCA_Log();
			}

		private static void write_Residue_cPCA_Log()
			{
				try
					{
						log_writer.write("\nPERFORMING Residue cPCA, Computing Top " + NUMBER_OF_MODES_RESIDUE + " modes." + "\n\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The coordinate values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("PCs Calculated Using: Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp" + "\n");
						log_writer.write("The COV PCs are the generalized coordinates for use in the RGC hierarchical PCA method" + "\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi_LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void do_Hierarchical_PCA()
			{

				JEDi_Do_Hierarchical_PCA cSS = new JEDi_Do_Hierarchical_PCA(DIRECTORY, DESCRIPTION, NUMBER_OF_MODES_HIERARCHCAL, Residue_Generalized_Coordinates_COV,
						Residue_Eigenvectors_COV, Residue_Delta_Vectors, NUMBER_OF_MODES_RESIDUE);

				cSS.do_Hierarchical_PCA();

				trace_COV = cSS.get_Trace();
				cond_cov = cSS.get_cond();
				det_cov = cSS.get_Det();
				rank_cov = cSS.get_Rank();
				top_cartesian_eigenvalues_COV = cSS.getTop_cartesian_eigenvalues();
				top_cartesian_evectors_COV = cSS.getTop_cartesian_evectors();
				square_pca_modes_COV = cSS.getSquare_pca_modes();
				pca_mode_max_COV = cSS.getPca_mode_max();
				pca_mode_min_COV = cSS.getPca_mode_min();

				projections_COV = cSS.getProjections();
				normed_projections_COV = cSS.getNormed_projections();
				weighted_normed_projections_COV = cSS.getWeighted_normed_projections();
				weighted_projections_COV = cSS.getWeighted_projections();

				G = cSS.getBig_G();
				G_Square_Modes = cSS.getG_Square_Modes();
				U = cSS.getBig_DV();
				hDVPs = cSS.getDVPs();

				System.gc();

				try
					{
						log_writer.write("\nPERFORMING Hierarchical PCA, Computing Top " + NUMBER_OF_MODES_HIERARCHCAL + " modes." + "\n\n");
						log_writer.write("Residue list for Hierarchical subset:  " + RESIDUE_LIST_HIERARCHICAL + "\n");
						log_writer.write("Number of residues in Cartesian subset: " + number_of_residues_Hierarchical_SS + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The coordinate values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_COV) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_cov) + "\n");
						log_writer.write("The DVPs from the COV model were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by sqrt of eigenvalue), and weighted normed dp" + "\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi_LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void do_All_Atom_cPCA()
			{

				JEDi_Do_Cartesian cSS = new JEDi_Do_Cartesian(DIRECTORY, DESCRIPTION, NUMBER_OF_MODES_ALL_ATOM_CARTESIAN, subset_PDB_coordinates_cartesian_AA,
						subset_reference_PDB_coordinates_cartesian_AA, AA);

				cSS.set_z_cutoff(Z_SCORE_CUTOFF);
				cSS.set_Out_dir(OUT_DIR);

				cSS.do_Cartesian();

				trace_COV = cSS.get_Trace_COV();
				trace_CORR = cSS.get_Trace_CORR();
				trace_PCORR = cSS.get_Trace_PCORR();
				cond_cov = cSS.get_cond_COV();
				det_cov = cSS.get_Det_COV();
				rank_cov = cSS.get_Rank_COV();
				transformed_residue_rmsds = cSS.getTransformed_residue_rmsd_list();
				top_cartesian_eigenvalues_COV = cSS.getTop_cartesian_eigenvalues_COV();
				top_cartesian_eigenvalues_CORR = cSS.getTop_cartesian_eigenvalues_CORR();
				top_cartesian_eigenvalues_PCORR = cSS.getTop_cartesian_eigenvalues_PCORR();
				top_cartesian_evectors_COV = cSS.getTop_cartesian_evectors_COV();
				top_cartesian_evectors_CORR = cSS.getTop_cartesian_evectors_CORR();
				top_cartesian_evectors_PCORR = cSS.getTop_cartesian_evectors_PCORR();
				square_pca_modes_COV = cSS.getSquare_pca_modes_COV();
				square_pca_modes_CORR = cSS.getSquare_pca_modes_CORR();
				square_pca_modes_PCORR = cSS.getSquare_pca_modes_PCORR();
				pca_mode_max_COV = cSS.getPca_mode_max_COV();
				pca_mode_max_CORR = cSS.getPca_mode_max_CORR();
				pca_mode_max_PCORR = cSS.getPca_mode_max_PCORR();
				pca_mode_min_COV = cSS.getPca_mode_min_COV();
				pca_mode_min_CORR = cSS.getPca_mode_min_CORR();
				pca_mode_min_PCORR = cSS.getPca_mode_min_PCORR();

				projections_CORR = cSS.getProjections_CORR();
				projections_COV = cSS.getProjections_COV();
				projections_PCORR = cSS.getProjections_PCORR();

				normed_projections_CORR = cSS.getNormed_projections_CORR();
				normed_projections_COV = cSS.getNormed_projections_COV();
				normed_projections_PCORR = cSS.getNormed_projections_PCORR();

				weighted_normed_projections_CORR = cSS.getWeighted_normed_projections_CORR();
				weighted_normed_projections_COV = cSS.getWeighted_normed_projections_COV();

				weighted_projections_CORR = cSS.getWeighted_projections_CORR();
				weighted_projections_COV = cSS.getWeighted_projections_COV();

				System.gc();

				Vector<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSDs(ref_subset_atoms_AA, transformed_residue_rmsds);
				String ref_PDB_file_name = "ss_" + number_of_residues_AA_Cart_SS + "_RMSF_edited.pdb";
				path = OUT_DIR + ref_PDB_file_name;
				PDB_IO.Write_PDB(path, edited_atoms);

				try
					{
						log_writer.write("\nPERFORMING All-Atom cPCA, Computing Top " + NUMBER_OF_MODES_ALL_ATOM_CARTESIAN + " modes." + "\n\n");
						log_writer.write("Residue list used for the All-Atom Cartesian subset:  " + RESIDUE_LIST_ALL_ATOM + "\n");
						log_writer.write("Number of Residues in the All-Atom Cartesian subset: " + number_of_residues_AA_Cart_SS + "\n");
						log_writer.write("Number of Atoms in the All-Atom Cartesian subset: " + number_of_atoms_AA_Cart_SS + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The coordinate values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_COV) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_cov) + "\n");
						log_writer.write("Trace of the Correlation Matrix = " + nff0.format(trace_CORR) + "\n");
						log_writer.write("Trace of the Partial Correlation Matrix = " + nff0.format(trace_PCORR) + "\n");
						log_writer.write("PDB file with B-factors replaced by residue RMSDs: ss_" + number_of_residues_AA_Cart_SS + "_RMSF_edited.pdb" + "\n");
						log_writer.write("The DVPs from the 3 different models were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp" + "\n");
						log_writer.write("Subspace analysis was done comparing the top vector spaces from the 3 different models." + "\n");
						log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to "
								+ NUMBER_OF_MODES_ALL_ATOM_CARTESIAN + "\n");
						log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
				do_Cartesian_SSA(AA);
			}

		private static void do_Heavy_Atom_cPCA()
			{

				JEDi_Do_Cartesian cSS = new JEDi_Do_Cartesian(DIRECTORY, DESCRIPTION, NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN, subset_PDB_coordinates_cartesian_HA,
						subset_reference_PDB_coordinates_cartesian_HA, HA);

				cSS.set_z_cutoff(Z_SCORE_CUTOFF);
				cSS.set_Out_dir(OUT_DIR);

				cSS.do_Cartesian();

				trace_COV = cSS.get_Trace_COV();
				trace_CORR = cSS.get_Trace_CORR();
				trace_PCORR = cSS.get_Trace_PCORR();
				cond_cov = cSS.get_cond_COV();
				det_cov = cSS.get_Det_COV();
				rank_cov = cSS.get_Rank_COV();
				transformed_residue_rmsds = cSS.getTransformed_residue_rmsd_list();
				top_cartesian_eigenvalues_COV = cSS.getTop_cartesian_eigenvalues_COV();
				top_cartesian_eigenvalues_CORR = cSS.getTop_cartesian_eigenvalues_CORR();
				top_cartesian_eigenvalues_PCORR = cSS.getTop_cartesian_eigenvalues_PCORR();
				top_cartesian_evectors_COV = cSS.getTop_cartesian_evectors_COV();
				top_cartesian_evectors_CORR = cSS.getTop_cartesian_evectors_CORR();
				top_cartesian_evectors_PCORR = cSS.getTop_cartesian_evectors_PCORR();
				square_pca_modes_COV = cSS.getSquare_pca_modes_COV();
				square_pca_modes_CORR = cSS.getSquare_pca_modes_CORR();
				square_pca_modes_PCORR = cSS.getSquare_pca_modes_PCORR();
				pca_mode_max_COV = cSS.getPca_mode_max_COV();
				pca_mode_max_CORR = cSS.getPca_mode_max_CORR();
				pca_mode_max_PCORR = cSS.getPca_mode_max_PCORR();
				pca_mode_min_COV = cSS.getPca_mode_min_COV();
				pca_mode_min_CORR = cSS.getPca_mode_min_CORR();
				pca_mode_min_PCORR = cSS.getPca_mode_min_PCORR();

				projections_CORR = cSS.getProjections_CORR();
				projections_COV = cSS.getProjections_COV();
				projections_PCORR = cSS.getProjections_PCORR();

				normed_projections_CORR = cSS.getNormed_projections_CORR();
				normed_projections_COV = cSS.getNormed_projections_COV();
				normed_projections_PCORR = cSS.getNormed_projections_PCORR();

				weighted_normed_projections_CORR = cSS.getWeighted_normed_projections_CORR();
				weighted_normed_projections_COV = cSS.getWeighted_normed_projections_COV();

				weighted_projections_CORR = cSS.getWeighted_projections_CORR();
				weighted_projections_COV = cSS.getWeighted_projections_COV();

				System.gc();

				Vector<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSDs(ref_subset_atoms_HA, transformed_residue_rmsds);
				String ref_PDB_file_name = "ss_" + number_of_residues_AA_Cart_SS + "_HA_RMSF_edited.pdb";
				path = OUT_DIR + ref_PDB_file_name;
				PDB_IO.Write_PDB(path, edited_atoms);

				try
					{
						log_writer.write("\nPERFORMING Heavy-Atom cPCA, Computing Top " + NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN + " modes." + "\n\n");
						log_writer.write("Residue list used for the Heavy-Atom Cartesian subset:  " + RESIDUE_LIST_HEAVY_ATOM + "\n");
						log_writer.write("Number of Residues in the Heavy-Atom Cartesian subset: " + number_of_residues_HA_Cart_SS + "\n");
						log_writer.write("Number of Atoms in the Heavy-Atom Cartesian subset: " + number_of_atoms_HA_Cart_SS + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The coordinate values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_COV) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_cov) + "\n");
						log_writer.write("Trace of the Correlation Matrix = " + nff0.format(trace_CORR) + "\n");
						log_writer.write("Trace of the Partial Correlation Matrix = " + nff0.format(trace_PCORR) + "\n");
						log_writer.write("PDB file with B-factors replaced by residue RMSDs: ss_" + number_of_residues_HA_Cart_SS + "_RMSF_edited.pdb" + "\n");
						log_writer.write("The DVPs from the 3 different models were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp" + "\n");
						log_writer.write("Subspace analysis was done comparing the top vector spaces from the 3 different models." + "\n");
						log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to "
								+ NUMBER_OF_MODES_HEAVY_ATOM_CARTESIAN + "\n");
						log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				do_Cartesian_SSA(HA);
			}

		private static void do_Alpha_Carbon_cPCA()
			{

				JEDi_Do_Cartesian cSS = new JEDi_Do_Cartesian(DIRECTORY, DESCRIPTION, NUMBER_OF_MODES_ALL_ATOM_CARTESIAN, subset_PDB_coordinates_cartesian_CA,
						subset_reference_PDB_coordinates_cartesian_CA, AC);

				cSS.set_z_cutoff(Z_SCORE_CUTOFF);
				cSS.set_Out_dir(OUT_DIR);

				cSS.do_Cartesian();

				trace_COV = cSS.get_Trace_COV();
				trace_CORR = cSS.get_Trace_CORR();
				trace_PCORR = cSS.get_Trace_PCORR();
				cond_cov = cSS.get_cond_COV();
				det_cov = cSS.get_Det_COV();
				rank_cov = cSS.get_Rank_COV();
				transformed_residue_rmsds = cSS.getTransformed_residue_rmsd_list();
				top_cartesian_eigenvalues_COV = cSS.getTop_cartesian_eigenvalues_COV();
				top_cartesian_eigenvalues_CORR = cSS.getTop_cartesian_eigenvalues_CORR();
				top_cartesian_eigenvalues_PCORR = cSS.getTop_cartesian_eigenvalues_PCORR();
				top_cartesian_evectors_COV = cSS.getTop_cartesian_evectors_COV();
				top_cartesian_evectors_CORR = cSS.getTop_cartesian_evectors_CORR();
				top_cartesian_evectors_PCORR = cSS.getTop_cartesian_evectors_PCORR();
				square_pca_modes_COV = cSS.getSquare_pca_modes_COV();
				square_pca_modes_CORR = cSS.getSquare_pca_modes_CORR();
				square_pca_modes_PCORR = cSS.getSquare_pca_modes_PCORR();
				pca_mode_max_COV = cSS.getPca_mode_max_COV();
				pca_mode_max_CORR = cSS.getPca_mode_max_CORR();
				pca_mode_max_PCORR = cSS.getPca_mode_max_PCORR();
				pca_mode_min_COV = cSS.getPca_mode_min_COV();
				pca_mode_min_CORR = cSS.getPca_mode_min_CORR();
				pca_mode_min_PCORR = cSS.getPca_mode_min_PCORR();

				projections_CORR = cSS.getProjections_CORR();
				projections_COV = cSS.getProjections_COV();
				projections_PCORR = cSS.getProjections_PCORR();

				normed_projections_CORR = cSS.getNormed_projections_CORR();
				normed_projections_COV = cSS.getNormed_projections_COV();
				normed_projections_PCORR = cSS.getNormed_projections_PCORR();

				weighted_normed_projections_CORR = cSS.getWeighted_normed_projections_CORR();
				weighted_normed_projections_COV = cSS.getWeighted_normed_projections_COV();

				weighted_projections_CORR = cSS.getWeighted_projections_CORR();
				weighted_projections_COV = cSS.getWeighted_projections_COV();

				System.gc();

				Vector<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSDs(ref_subset_atoms_CA, transformed_residue_rmsds);
				String ref_PDB_file_name = "ss_" + number_of_residues_AA_Cart_SS + "_CA_RMSF_edited.pdb";
				path = OUT_DIR + ref_PDB_file_name;
				PDB_IO.Write_PDB(path, edited_atoms);

				try
					{
						log_writer.write("\nPERFORMING Alpha-Carbon cPCA, Computing Top " + NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN + " modes." + "\n\n");
						log_writer.write("Residue list used for the Alpha-Carbon Cartesian subset:  " + RESIDUE_LIST_ALPHA_CARBON + "\n");
						log_writer.write("Number of Residues in the Alpha-Carbon Cartesian subset: " + number_of_residues_CA_Cart_SS + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The coordinate values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_COV) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_cov) + "\n");
						log_writer.write("Trace of the Correlation Matrix = " + nff0.format(trace_CORR) + "\n");
						log_writer.write("Trace of the Partial Correlation Matrix = " + nff0.format(trace_PCORR) + "\n");
						log_writer.write("PDB file with B-factors replaced by residue RMSDs: ss_" + number_of_residues_CA_Cart_SS + "_CA_RMSF_edited.pdb" + "\n");
						log_writer.write("The DVPs from the 3 different models were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp" + "\n");
						log_writer.write("Subspace analysis was done comparing the top vector spaces from the 3 different models." + "\n");
						log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to "
								+ NUMBER_OF_MODES_ALPHA_CARBON_CARTESIAN + "\n");
						log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				do_Cartesian_SSA(AC);
			}

		private static void do_dpPCA()
			{
				JEDi_Do_Atom_Dist_Pairs dpSS = new JEDi_Do_Atom_Dist_Pairs(DIRECTORY, DESCRIPTION, NUMBER_OF_MODES_DISTANCE_PAIRS, original_PDB_coordinates_AA,
						original_reference_PDB_coordinates_AA, atom_list_dp1, atom_list_dp_original1, atom_list_dp2, atom_list_dp_original2, chain_idents1, chain_idents2);
				dpSS.set_z_cutoff(Z_SCORE_CUTOFF);

				if (!doMULTI)
					dpSS.do_Dist();

				if (doMULTI)
					dpSS.do_Dist_Multi();

				trace_d_cov = dpSS.getTrace_dist_COV();
				cond_d_cov = dpSS.getCond_cov();
				det_d_cov = dpSS.getDet_cov();
				rank_d_cov = dpSS.getRank_cov();
				trace_d_corr = dpSS.getTrace_dist_CORR();
				trace_d_pcorr = dpSS.getTrace_dist_PCORR();

				top_distance_evectors_COV = dpSS.getTop_distance_evectors_COV();
				top_distance_evectors_CORR = dpSS.getTop_distance_evectors_CORR();
				top_distance_evectors_PCORR = dpSS.getTop_distance_evectors_PCORR();

				projections_dist_COV = dpSS.getProjections_dist_COV();
				normed_projections_dist_COV = dpSS.getNormed_projections_dist_COV();
				weighted_projections_dist_COV = dpSS.getWeighted_projections_dist_COV();
				weighted_normed_projections_dist_COV = dpSS.getWeighted_normed_projections_dist_COV();

				projections_dist_CORR = dpSS.getProjections_dist_CORR();
				normed_projections_dist_CORR = dpSS.getNormed_projections_dist_CORR();
				weighted_projections_dist_CORR = dpSS.getWeighted_projections_dist_CORR();
				weighted_normed_projections_dist_CORR = dpSS.getWeighted_normed_projections_dist_CORR();

				projections_dist_PCORR = dpSS.getProjections_dist_PCORR();
				normed_projections_dist_PCORR = dpSS.getNormed_projections_dist_PCORR();

				residue_distance_means = dpSS.getResidue_distance_means();
				residue_distance_std_devs = dpSS.getResidue_distance_std_devs();

				System.gc();

				if (!doMULTI)
					write_dpPCA_Log_Single();
				if (doMULTI)
					write_dpPCA_Log_Multi();

				do_Distance_Pair_SSA(DP);
			}

		private static void write_dpPCA_Log_Multi()
			{
				try
					{
						log_writer.write("\nPERFORMING dpPCA, Computing Top " + NUMBER_OF_MODES_DISTANCE_PAIRS + " modes.\n");
						log_writer.write("\nAtom Pair list:  " + ATOM_PAIRS_LIST + "\n");
						log_writer.write("Number of atom pairs: " + number_of_atom_pairs + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The distance values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_d_cov) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_d_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_d_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_d_cov) + "\n");
						log_writer.write("Trace of the Correlation Matrix = " + nff0.format(trace_d_corr) + "\n");
						log_writer.write("Trace of the Partial Correlation Matrix = " + nff0.format(trace_d_pcorr) + "\n");
						log_writer.write("\nMEANs and STANDARD DEVIATIONs for the Atom Pair Distances: " + "\n\n");
						log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Atom1", "Atom2", "Mean", "Std_Dev"));
						for (int i = 0; i < number_of_atom_pairs; i++)
							{
								log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", chain_idents1.get(i) + atom_list_dp_original1.get(i),
										chain_idents2.get(i) + atom_list_dp_original2.get(i), nf3.format(residue_distance_means[i]), nf6.format(residue_distance_std_devs[i])));
							}
						log_writer.write("\nThe DVPs (PCs) from from the 3 different models were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp." + "\n");
						log_writer.write("Subspace analysis was done to compare the top vector spaces from the 3 different models." + "\n");
						log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to "
								+ NUMBER_OF_MODES_DISTANCE_PAIRS + "\n");
						log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JED_LOG file: " + OUT_DIR + "JED_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void write_dpPCA_Log_Single()
			{
				try
					{
						log_writer.write("\nPERFORMING dpPCA, Computing Top " + NUMBER_OF_MODES_DISTANCE_PAIRS + " modes.\n");
						log_writer.write("\nAtom Pair list:  " + ATOM_PAIRS_LIST + "\n");
						log_writer.write("Number of atom pairs: " + number_of_atom_pairs + "\n");
						if (Z_SCORE_CUTOFF > 0)
							log_writer.write("The distance values with Z-scores beyond |" + Z_SCORE_CUTOFF + "| were set to their mean value.\n\n");
						else
							log_writer.write("No coordinate outliers were adjusted.\n");
						log_writer.write("Trace of the Covariance Matrix = " + df.format(trace_d_cov) + "\n");
						log_writer.write("Condition Number of the Covariance Matrix = " + df.format(cond_d_cov) + "\n");
						log_writer.write("Determinant of the Covariance Matrix = " + df.format(det_d_cov) + "\n");
						log_writer.write("Rank of the Covariance Matrix = " + nff0.format(rank_d_cov) + "\n");
						log_writer.write("Trace of the Correlation Matrix = " + nff0.format(trace_d_corr) + "\n");
						log_writer.write("Trace of the Partial Correlation Matrix = " + nff0.format(trace_d_pcorr) + "\n");
						log_writer.write("\nMEANs and STANDARD DEVIATIONs for the Residue Pair Distances: " + "\n\n");
						log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", "Atom1", "Atom2", "Mean", "Std_Dev"));
						for (int i = 0; i < number_of_atom_pairs; i++)
							{
								log_writer.write(String.format("%-12s%-16s%-16s%-16s%n", atom_list_dp_original1.get(i), atom_list_dp_original2.get(i),
										nf3.format(residue_distance_means[i]), nf6.format(residue_distance_std_devs[i])));
							}
						log_writer.write("\nThe DVPs (PCs) from from the 3 different models were calculated using:" + "\n");
						log_writer.write("Standard dot product(dp), normed dp, weighted dp (by eigenvalue), and weighted normed dp." + "\n");
						log_writer.write("Subspace analysis was done to compare the top vector spaces from the 3 different models." + "\n");
						log_writer.write("Comparators include RMSIP and Principle Angles, for the essential subspace and iterated comparisons from dim 1 to "
								+ NUMBER_OF_MODES_DISTANCE_PAIRS + "\n");
						log_writer.write("Additional log files can be found in the /SSA directory tree." + "\n\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JED_LOG file: " + OUT_DIR + "JED_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		public static void doNoPCA()
			{
				JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_coordinates_CA, original_reference_PDB_coordinates_CA, DIRECTORY,
						DESCRIPTION);
					{
						tf_coords.set_z_cutoff(Z_SCORE_CUTOFF);
						tf_coords.set_Output_Directory(DIRECTORY + "JEDi_RESULTS_" + DESCRIPTION + File.separatorChar + "NoPCA" + File.separatorChar);
						tf_coords.get_Transformed_reference_coordinates();
						tf_coords.get_SS_Transformed_coords();
						tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
						tf_coords.get_SS_Conformation_RMSDs();
						tf_coords.get_conf_Z_scores();
						transformed_residue_rmsds = tf_coords.get_SS_RMSF();
					}

				Vector<Atom> edited_atoms = refPDB.edit_B_Factors_with_RMSDs(ref_atoms_CA, transformed_residue_rmsds);
				String ref_PDB_file_name = "ss_" + number_of_residues_REF + "_CA_RMSF_edited.pdb";
				path = OUT_DIR + ref_PDB_file_name;
				PDB_IO.Write_PDB(path, edited_atoms);

				tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_coordinates_AA, original_reference_PDB_coordinates_AA, DIRECTORY, DESCRIPTION);
					{
						tf_coords.set_z_cutoff(Z_SCORE_CUTOFF);
						tf_coords.set_Output_Directory(DIRECTORY + "JEDi_RESULTS_" + DESCRIPTION + File.separatorChar + "NoPCA" + File.separatorChar);
						tf_coords.get_Transformed_reference_coordinates();
						tf_coords.get_SS_Transformed_coords();
						tf_coords.get_SS_transformed_coordinates_adjusted_ROWS();
						tf_coords.get_SS_Conformation_RMSDs();
						tf_coords.get_conf_Z_scores();
						transformed_residue_rmsds = tf_coords.get_SS_RMSF();
					}

				edited_atoms = refPDB.edit_B_Factors_with_RMSDs(ref_atoms_AA, transformed_residue_rmsds);
				ref_PDB_file_name = "ss_" + number_of_residues_REF + "_RMSF_edited.pdb";
				path = OUT_DIR + ref_PDB_file_name;
				PDB_IO.Write_PDB(path, edited_atoms);

				try
					{
						log_writer.write("\nConformation RMSDs and Residue RMSDs (RMSF) were calculated." + "\n");
						log_writer.write("\nPDB files with B-factors replaced by residue RMSF were created: \n");
						log_writer.write("\tss_" + number_of_residues_REF + "_RMSF_edited.pdb" + "\n");
						log_writer.write("\tss_" + number_of_residues_REF + "_CA_RMSF_edited.pdb" + "\n");
						log_writer.write("\nVariable Z-scores were calculated.\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void do_Mode_Visualization_All_Atom()
			{

				JEDi_PCA_Mode_Vizualization cov = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_AA, top_cartesian_eigenvalues_COV,
						top_cartesian_evectors_COV, square_pca_modes_COV, pca_mode_max_COV, pca_mode_min_COV, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AA, Q);
				cov.get_Mode_Visualizations_SS_All_Atom();
				cov.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization corr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_AA, top_cartesian_eigenvalues_CORR,
						top_cartesian_evectors_CORR, square_pca_modes_CORR, pca_mode_max_CORR, pca_mode_min_CORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AA, R);
				corr.get_Mode_Visualizations_SS_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization pcorr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_AA, top_cartesian_eigenvalues_PCORR,
						top_cartesian_evectors_PCORR, square_pca_modes_PCORR, pca_mode_max_PCORR, pca_mode_min_PCORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AA, P);
				pcorr.get_Mode_Visualizations_SS_All_Atom();
				pcorr.get_Essential_Visualization_All_Atom();

				System.gc();
				write_VIZ_Log(AA);
			}

		private static void do_Mode_Visualization_CA()
			{

				JEDi_PCA_Mode_Vizualization cov = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_CA, top_cartesian_eigenvalues_COV,
						top_cartesian_evectors_COV, square_pca_modes_COV, pca_mode_max_COV, pca_mode_min_COV, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AC, Q);
				cov.get_Mode_Visualizations_SS_All_Atom();
				cov.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization corr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_CA, top_cartesian_eigenvalues_CORR,
						top_cartesian_evectors_CORR, square_pca_modes_CORR, pca_mode_max_CORR, pca_mode_min_CORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AC, R);
				corr.get_Mode_Visualizations_SS_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization pcorr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_CA, top_cartesian_eigenvalues_PCORR,
						top_cartesian_evectors_PCORR, square_pca_modes_PCORR, pca_mode_max_PCORR, pca_mode_min_PCORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, AC, P);
				pcorr.get_Mode_Visualizations_SS_All_Atom();
				pcorr.get_Essential_Visualization_All_Atom();

				System.gc();
				write_VIZ_Log(AC);
			}

		private static void do_Mode_Visualization_HA()
			{

				JEDi_PCA_Mode_Vizualization cov = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_HA, top_cartesian_eigenvalues_COV,
						top_cartesian_evectors_COV, square_pca_modes_COV, pca_mode_max_COV, pca_mode_min_COV, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, HA, Q);
				cov.get_Mode_Visualizations_SS_All_Atom();
				cov.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization corr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_HA, top_cartesian_eigenvalues_CORR,
						top_cartesian_evectors_CORR, square_pca_modes_CORR, pca_mode_max_CORR, pca_mode_min_CORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, HA, R);
				corr.get_Mode_Visualizations_SS_All_Atom();
				corr.get_Essential_Visualization_All_Atom();

				System.gc();

				JEDi_PCA_Mode_Vizualization pcorr = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_HA, top_cartesian_eigenvalues_PCORR,
						top_cartesian_evectors_PCORR, square_pca_modes_PCORR, pca_mode_max_PCORR, pca_mode_min_PCORR, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, HA, P);
				pcorr.get_Mode_Visualizations_SS_All_Atom();
				pcorr.get_Essential_Visualization_All_Atom();

				System.gc();
				write_VIZ_Log(HA);
			}

		private static void do_Mode_Visualization_Hierarchical()
			{

				JEDi_PCA_Mode_Vizualization cov = new JEDi_PCA_Mode_Vizualization(DIRECTORY, DESCRIPTION, ref_subset_atoms_RES_HIER, top_cartesian_eigenvalues_COV, G,
						G_Square_Modes, pca_mode_max_COV, pca_mode_min_COV, VIZ_MODE_AMPLITUDE, NUMBER_OF_MODES_VIZ, RGC, Q);
				cov.get_Mode_Visualizations_SS_All_Atom();
				cov.get_Essential_Visualization_All_Atom();

				System.gc();
				write_VIZ_Log(RGC);
			}

		private static void write_VIZ_Log(String type)
			{
				try
					{
						log_writer.write("Performing " + type + " Cartesian Mode Visualization on Top  " + NUMBER_OF_MODES_VIZ + "  cPCA modes.\n");
						log_writer.write("Sets of 20 structures were generated to animate each selected cPCA mode, for the available PCA models.\n");
						log_writer.write("\tNote that the Hierarchical method only uses the Covariance PCA model.\n");
						if (NUMBER_OF_MODES_VIZ >= 5)
							log_writer.write("A set of 20 structures was generated to animate the essential subspace composed of the top 5 modes for the available PCA models.\n");
						if (NUMBER_OF_MODES_VIZ < 5)
							log_writer.write("A set of 20 structures was generated to animate the essential subspace composed of the top " + NUMBER_OF_MODES_VIZ
									+ " modes, for the COV, CORR, and PCORR PCA models.\n");
						log_writer.write("For the individual modes, atoms of each residue were perturbed along the mode eigenvector using a sine function over one full period.\n");
						log_writer.write("For the essential motion, a super-position of the top modes was created spanning 2 full periods for the first mode.\n");
						log_writer.write("Pymol(TM) scripts were generated for each individual mode and the essential subspace, to play the mode structures as movies.\n");
						log_writer.write("MODE AMPLITUDE = " + nf3.format(VIZ_MODE_AMPLITUDE) + "\n");
						log_writer.flush();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi_LOG file: " + OUT_DIR + "JEDi_LOG.txt");
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}
				System.gc();
			}

		private static void do_Cartesian_SSA(String type_pca)
			{
				JEDi_Get_Subspace_Analysis jssa1 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "COV_vs_CORR", top_cartesian_evectors_CORR,
						top_cartesian_evectors_COV);
					{
						jssa1.get_SSA_JED();
						jssa1.get_FSSA_Iterated_JED();
						System.gc();
					}
				JEDi_Get_Subspace_Analysis jssa2 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "COV_vs_PCORR", top_cartesian_evectors_PCORR,
						top_cartesian_evectors_COV);
					{
						jssa2.get_SSA_JED();
						jssa2.get_FSSA_Iterated_JED();
						System.gc();
					}
				JEDi_Get_Subspace_Analysis jssa3 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "CORR_vs_PCORR", top_cartesian_evectors_PCORR,
						top_cartesian_evectors_CORR);
					{
						jssa3.get_SSA_JED();
						jssa3.get_FSSA_Iterated_JED();
						System.gc();
					}
			}

		private static void do_Distance_Pair_SSA(String type_pca)
			{
				JEDi_Get_Subspace_Analysis jssa1 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "COV_vs_CORR", top_distance_evectors_CORR,
						top_distance_evectors_COV);
					{
						jssa1.get_SSA_JED();
						jssa1.get_FSSA_Iterated_JED();
						System.gc();
					}
				JEDi_Get_Subspace_Analysis jssa2 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "COV_vs_PCORR", top_distance_evectors_PCORR,
						top_distance_evectors_COV);
					{
						jssa2.get_SSA_JED();
						jssa2.get_FSSA_Iterated_JED();
						System.gc();
					}
				JEDi_Get_Subspace_Analysis jssa3 = new JEDi_Get_Subspace_Analysis(DIRECTORY, DESCRIPTION, type_pca, "CORR_vs_PCORR", top_distance_evectors_PCORR,
						top_distance_evectors_CORR);
					{
						jssa3.get_SSA_JED();
						jssa3.get_FSSA_Iterated_JED();
						System.gc();
					}
			}

		private static void write_FES_Log()
			{
				if (do_cartesian_AA || do_cartesian_HA || do_cartesian_CA || do_hierarchical || do_dist_pairs)
					{
						try
							{
								log_writer.write("\nPerforming Free Energy calculations using two order parameters to create a Free-Energy Surface: \n");
								log_writer.write("\tThe first and second PC (weighted-projection). \n");
								log_writer.write("See the FES Logs in the FES subdirectory for more details.\n\n");
								log_writer.flush();

							} catch (IOException e)
							{
								System.err.println("Could not write the JEDi_LOG file: " + OUT_DIR + "JEDi_LOG.txt");
								System.err.println("Program terminating.\n");
								e.printStackTrace();
								System.exit(0);
							}
						System.gc();
					}
			}

		// ************************************** MAIN METHOD ************************************************************* //

		public static void main(String[] args)
			{
				jedi = new JEDi_Driver();
				System.out.println("Running JEDi Driver: ");
				System.out.println("Reading the input file: ");
				try
					{
						String WD = System.getProperty("user.dir");
						String in_path = WD + File.separator + input_file;
						if (args.length == 1)
							{
								in_path = args[0];
								System.out.println("User specified input file (must be the first program argument)");
								System.out.println("These are the specified program args:");
								for (int i = 0; i < args.length; i++)
									{
										System.out.println("Arg " + (i + 1) + " Value = " + args[i]);
									}
							}
						System.out.println("Working Directory = " + WD);
						System.out.println("Input File Path = " + in_path);
						boolean check = new File(in_path).exists();
						if (!check)
							{
								System.err.println("The Input File can not be found: " + in_path);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
						input_reader = new BufferedReader(new FileReader(in_path));

					} catch (FileNotFoundException e)
					{
						System.err.println("Could not find the input file: " + input_file);
						System.err.println("Program terminating.\n");
						e.printStackTrace();
						System.exit(0);
					}

				System.out.println("Assigning Paramaters from Input File... ");
				read_input_file();
				assign_parameters();
				initialize_JED_Log();
				System.out.println("Done.");

				System.out.println("Processing PDB Reference File... ");
				startTime = System.nanoTime();
				process_Reference_PDB();
				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");

				if (doREAD)
					{
						System.out.print("Reading PDB files... ");
						startTime = System.nanoTime();
						read_PDB_Files();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
					}

				if (!doREAD)
					{
						System.out.print("Reading Coordinate Files... ");
						startTime = System.nanoTime();
						read_All_Atom_Coordinate_File();
						read_HA_Coordinate_File();
						read_CA_Coordinate_File();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");

						System.out.print("Getting Residue Coordinate Matrices... ");
						startTime = System.nanoTime();
						Residue_Coordinates = refPDB.get_Residue_Coords(original_PDB_coordinates_AA, numbers_Of_Atoms_in_Residues);
						Residue_Coordinates_HA = refPDB.get_Residue_Coords_HA(original_PDB_coordinates_HA, numbers_Of_Heavy_Atoms_in_Residues);
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");

						System.out.print("Getting Required Residue / Atom subsets... ");
						startTime = System.nanoTime();
						get_Subsets();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
					}

				if (do_residue)
					{
						System.out.print("Performing Residue cPCA analysis... ");
						startTime = System.nanoTime();
						do_Residue_cPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
					}

				if (do_hierarchical)
					{
						System.out.print("Performing Hierarchical PCA analysis... ");
						startTime = System.nanoTime();
						do_Hierarchical_PCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
						if (doVIZ)
							{
								System.out.print("Performing Hierarchical PCA mode visualization... ");
								startTime = System.nanoTime();
								do_Mode_Visualization_Hierarchical();
								endTime = System.nanoTime();
								totalTime = endTime - startTime;
								System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
							}
					}

				if (do_cartesian_AA)
					{
						System.out.print("Performing All-Atom cPCA analysis... ");
						startTime = System.nanoTime();
						do_All_Atom_cPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
						if (doVIZ)
							{
								System.out.print("Performing All Atom cPCA mode visualization... ");
								startTime = System.nanoTime();
								do_Mode_Visualization_All_Atom();
								write_VIZ_Log(AA);
								endTime = System.nanoTime();
								totalTime = endTime - startTime;
								System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
							}
					}

				if (do_cartesian_HA)
					{
						System.out.print("Performing Heavy-Atom cPCA analysis... ");
						startTime = System.nanoTime();
						do_Heavy_Atom_cPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
						if (doVIZ)
							{
								System.out.print("Performing Heavy-Atom cPCA mode visualization... ");
								startTime = System.nanoTime();
								do_Mode_Visualization_HA();
								write_VIZ_Log(HA);
								endTime = System.nanoTime();
								totalTime = endTime - startTime;
								System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
							}
					}

				if (do_cartesian_CA)
					{
						System.out.print("Performing Alpha-Carbon cPCA analysis... ");
						startTime = System.nanoTime();
						do_Alpha_Carbon_cPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
						if (doVIZ)
							{
								System.out.print("Performing Alpha-Carbon cPCA mode visualization... ");
								startTime = System.nanoTime();
								do_Mode_Visualization_CA();
								write_VIZ_Log(AC);
								endTime = System.nanoTime();
								totalTime = endTime - startTime;
								System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
							}
					}

				if (do_dist_pairs)
					{
						System.out.print("Performing dpPCA analysis... ");
						startTime = System.nanoTime();
						do_dpPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
					}

				if (do_no_pca)
					{
						System.out.print("Pre-Processing PDB files for JEDi Analysis... ");
						startTime = System.nanoTime();
						doNoPCA();
						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done. (" + nf3.format(totalTime / 1000000000.0) + " seconds)");
					}

				write_FES_Log();
				try
					{
						date = DateUtils.now();
						log_writer.write("\nAnalysis completed: " + date);
						log_writer.close();

					} catch (IOException e)
					{
						System.err.println("Could not write the JEDi_LOG file");
						e.printStackTrace();
					}
				System.out.println("JEDi successfully completed: " + date);
			}
	}
