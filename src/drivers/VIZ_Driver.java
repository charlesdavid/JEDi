package drivers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.StringTokenizer;

import Jama.Matrix;
import support.Atom;
import support.Projector;
import supportIO.DateUtils;
import supportIO.FortranFormat;
import supportIO.Input_Parameters_VIZ;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportIO.PDB_File_Parser;
import supportIO.PDB_IO;

/**
 * JEDi class VIZ_Driver: Driver program to visualize Cartesian PCA modes.
 * 
 * Reads input file: "VIZ_Parameters.txt"
 * 
 * The expectation is that the PDB file corresponds EXACTLY to the eigenvector file.
 * 
 * Constructs sets of PDB files and Pymol(TM) Scripts to animate the individual PCA modes derived from the Cartesian subset.
 * 
 * Constructs a set of frames and a Pymol(TM) Script to animate the essential motion from a combined group of PCA modes.
 * 
 * NOTE1: This program will NOT animate modes from alpha carbon analysis only. The pdb MUST have at least back bone atoms.
 * 
 * NOTE2: Log coloring based on pca modes: SQRT(x^2 + y^2 + z^2) * SQRT(eigenvalue), so units are angstroms.
 * 
 * NOTE3: Affine increase of mode amplitude means that large proteins will need to be viewed as cartoon or ribbon, not sticks.
 * 
 * Copyright (C) 2012, 2020 Dr. Charles David
 * 
 * @author Dr. Charles David
 * 
 */
public class VIZ_Driver
{
	static String input_file = "VIZ_Parameters.txt", key, value, delim = "=", line, input_path, out_dir, PDB, eigenvectors, evals, OUT_INDIVIDUAL, OUT_ESSENTIAL, file_name_head,
			date, model, description;
	static int number_Of_Input_Lines, line_count, ROWS_Evectors, ROWS_Modes, COLS, number_of_atoms, number_of_modes;
	static final double delta_y = 99, pi = Math.PI;
	static List<Double> eigenvalues, pca_mode_maxs, pca_mode_mins;
	static Matrix top_evectors, top_square_pca_modes;
	static BufferedWriter output_file_writer;
	static List<Atom> atoms;
	static File log;
	static BufferedReader input_reader;
	static BufferedWriter log_writer;
	static DateUtils now;
	static StringTokenizer sToken;
	static NumberFormat nf;
	static RoundingMode rm;
	static long startTime, endTime, totalTime;
	static List<String> lines;
	static FortranFormat formatter;
	static PDB_File_Parser parser;
	static boolean do_individual, do_essential, exist, success;
	static Hashtable<String, String> parameters;
	static Input_Parameters_VIZ ip;

	private static void read_input_file()
	{
		number_Of_Input_Lines = 0;
		System.out.println("Reading file: " + input_file);
		lines = new ArrayList<String>();
		parameters = new Hashtable<>();
		try
			{
				while ((line = input_reader.readLine()) != null && line.length() >= 1)
					{
						if (!(line.startsWith("#")))
							{
								lines.add(line);
								sToken = new StringTokenizer(line);
								key = sToken.nextToken(delim);
								value = sToken.nextToken(delim);
								parameters.put(key, value);
								number_Of_Input_Lines++;
							}
					}
				input_reader.close();
				System.out.println("\tThe total number of KEY=VALUE Pairs in the input file is: " + number_Of_Input_Lines + "\n");
				ip = new Input_Parameters_VIZ(parameters);
			}
		catch (IOException e)
			{
				System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void write_Job_Log()
	{
		try
			{
				log = new File(out_dir + "VIZ_LOG.txt");
				log_writer = new BufferedWriter(new FileWriter(log));
				log_writer.write("JEDi Cartesian Mode Visualization Driver version 1.0" + "\n");
				log_writer.write("Description: " + Input_Parameters_VIZ.DESCRIPTION + "\n");
				log_writer.write("Number of Atoms in Subset: " + number_of_atoms + "\n");
				log_writer.write("PCA Model: " + model + "\n");
				log_writer.write("Reference PDB file: " + Input_Parameters_VIZ.REFERENCE_PDB + "\n");
				log_writer.write("Eigenvalues file: " + Input_Parameters_VIZ.EVALS + "\n");
				log_writer.write("Eigenvectors file: " + Input_Parameters_VIZ.EVECTS + "\n");
				log_writer.write("Output directory: " + Input_Parameters_VIZ.OUT_DIR + "\n");

				log_writer.write("First Mode = " + Input_Parameters_VIZ.MODE_START + "\n");
				log_writer.write("Last Mode = " + Input_Parameters_VIZ.MODE_END + "\n");
				log_writer.write("Number of Modes = " + number_of_modes + "\n");

				log_writer.write("Number of Frames = " + Input_Parameters_VIZ.numberOfFrames + "\n");
				log_writer.write("Number of Cycles = " + Input_Parameters_VIZ.numberOfCycles + "\n");
				log_writer.write("Number of Essential Frames = " + Input_Parameters_VIZ.numberOfFramesEssential + "\n");
				log_writer.write("Number of Essential Cycles = " + Input_Parameters_VIZ.numberOfCyclesEssential + "\n");
				log_writer.write("Mode amplitutde = " + Input_Parameters_VIZ.mode_amplitude + "\n");

				if (do_individual) log_writer.write("Performed Individual Motion Visualization. \n");
				if (do_essential) log_writer.write("Performed Essential Motion Visualization. \n");

				log_writer.write("\nMode Visualization Completed: " + date);
				log_writer.close();
			}

		catch (IOException e)
			{
				System.err.println("Could not write the VIZ_LOG file: " + out_dir + "VIZ_LOG.txt");
				System.err.println("Program terminating.\n");
				e.printStackTrace();
				System.exit(0);
			}
	}

	@SuppressWarnings("unchecked")
	private static void read_data_files()
	{
		PDB_IO pdbio = new PDB_IO(Input_Parameters_VIZ.REFERENCE_PDB);
		atoms = pdbio.parse_PDB_File();
		number_of_atoms = atoms.size();
		number_of_modes = Input_Parameters_VIZ.numberOfModes;
		top_evectors = Matrix_IO.read_Matrix(Input_Parameters_VIZ.EVECTS);
		eigenvalues = List_IO.read_List_From_File(Input_Parameters_VIZ.EVALS, "Double");

		ROWS_Evectors = top_evectors.getRowDimension();
		do_individual = Input_Parameters_VIZ.doINDIVIDUAL;
		do_essential = Input_Parameters_VIZ.doESSENTIAL;
	}

	private static void Output_Control()
	{
		model = Input_Parameters_VIZ.modelPCA;
		description = Input_Parameters_VIZ.DESCRIPTION;

		out_dir = Input_Parameters_VIZ.OUT_DIR + "JEDi_VIZ_Driver_" + description + "_ss_" + number_of_atoms + "_" + model + File.separatorChar;
		exist = new File(out_dir).exists();
		if (!exist)
			{
				success = (new File(out_dir)).mkdirs();
				if (!success)
					{
						System.err.println("Sorry, unable to create the output directory. Please check synatx and file permissions: " + out_dir);
						System.exit(0);
					}
			}

		OUT_INDIVIDUAL = out_dir + "Individual_Motion" + File.separatorChar;
		exist = new File(OUT_INDIVIDUAL).exists();
		if (!exist)
			{
				success = (new File(OUT_INDIVIDUAL)).mkdirs();
				if (!success)
					{
						System.err.println("Sorry, unable to create the output directory. Please check synatx and file permissions: " + OUT_INDIVIDUAL);
						System.exit(0);
					}
			}

		OUT_ESSENTIAL = out_dir + "Essential_Motion" + File.separatorChar;
		exist = new File(OUT_ESSENTIAL).exists();
		if (!exist)
			{
				success = (new File(OUT_ESSENTIAL)).mkdirs();
				if (!success)
					{
						System.err.println("Sorry, unable to create the output directory. Please check synatx and file permissions: " + OUT_ESSENTIAL);
						System.exit(0);
					}
			}
	}

	private static void construct_PCA_Modes()
	{
		ROWS_Modes = (ROWS_Evectors / 3);
		top_square_pca_modes = new Matrix(ROWS_Modes, number_of_modes);
		pca_mode_maxs = new ArrayList<Double>();
		pca_mode_mins = new ArrayList<Double>();

		for (int i = 0; i < number_of_modes; i++)
			{
				double max = 0;
				double min = 1;
				double value = Math.sqrt(eigenvalues.get(i));
				for (int j = 0; j < ROWS_Modes; j++)
					{
						double x = top_evectors.get(j, i);
						double y = top_evectors.get(j + ROWS_Modes, i);
						double z = top_evectors.get((j + 2 * ROWS_Modes), i);
						double sq = value * Math.sqrt(x * x + y * y + z * z);
						top_square_pca_modes.set(j, i, sq);
						if (sq > max) max = sq;
						if (sq < min) min = sq;
					}
				pca_mode_maxs.add(i, max);
				pca_mode_mins.add(i, min);
			}
	}

	private static void get_Individual_Mode_Visualizations_AA()
	{
		int start = Input_Parameters_VIZ.MODE_START - 1;
		for (int mode = start; mode < number_of_modes; mode++) // MODE LOOP
			{
				Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
				Matrix sq_mode = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, mode, mode);

				/* Establish the log coloring scheme using logs, Square Mode components as B-factors */
				double MODE_MIN = pca_mode_mins.get(mode);
				double MODE_MAX = pca_mode_maxs.get(mode);
				if (MODE_MIN < Input_Parameters_VIZ.FLOOR) MODE_MIN = Input_Parameters_VIZ.FLOOR;
				double LOG_MODE_MIN = Math.log10(MODE_MIN);
				double LOG_MODE_MAX = Math.log10(MODE_MAX);
				double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
				double slope = (delta_y / delta_x);

				/* Get the MAX eigenvector component. This controls the maximum deflection... */
				ArrayList<Double> evects = new ArrayList<Double>();
				for (int i = 0; i < ROWS_Evectors; i++)
					{
						double component = evector.get(i, 0);
						double compAbs = Math.abs(component);
						evects.add(compAbs);
					}
				Collections.sort(evects, Collections.reverseOrder());
				double Evector_Max = evects.get(0);
				double MA = Input_Parameters_VIZ.mode_amplitude;

				/* The mode amplitude should be adjusted to maximally stretch bonds without overly distorting them, about 1 Angstrom */
				if (number_of_atoms < 30) MA = ((MA / Evector_Max) + (number_of_atoms / 100));
				if (number_of_atoms >= 30) MA = ((MA / Evector_Max) + (number_of_atoms / 50));

				// System.out.println("\tMode " + (mode + 1) + " Mode Amplitude: " + df.format(mode_amplitude));

				String f_index = "";
				int frame_index = 0;

				for (int t = 0; t < Input_Parameters_VIZ.numberOfFrames; t++) // FRAME LOOP:
					{
						// System.out.println("Computing Frame: " + (frame_index + 1));

						double omega = ((2 * Input_Parameters_VIZ.numberOfCycles * pi / Input_Parameters_VIZ.numberOfFrames));
						/* A Sine function ensures the first structure is unperturbed; Preserves chain connectivity in PyMol movies... */
						double weight = Math.sin(omega * t);

						List<Atom> shifted_atoms = new ArrayList<Atom>(atoms.size()); // create a new list to hold the shifted atoms
						for (Atom a : atoms)
							{
								Atom sa = new Atom(a);
								shifted_atoms.add(sa);
							}

						int index = 0;
						for (Atom a : shifted_atoms) // ATOM LOOP: Iterates through the list of atoms shifting the atomic coordinates;
							{
								double x_coord = a.getX();
								double y_coord = a.getY();
								double z_coord = a.getZ();

								double v_x = evector.get(index, 0);
								double v_y = evector.get((index + number_of_atoms), 0);
								double v_z = evector.get((index + 2 * number_of_atoms), 0);

								double shift_x = (v_x * weight * MA);
								double shift_y = (v_y * weight * MA);
								double shift_z = (v_z * weight * MA);

								a.setX(x_coord + shift_x);
								a.setY(y_coord + shift_y);
								a.setZ(z_coord + shift_z);

								double sq_mode_val = (sq_mode.get(index, 0));
								if (sq_mode_val < MODE_MIN) sq_mode_val = MODE_MIN;
								if (sq_mode_val > MODE_MAX) sq_mode_val = MODE_MAX;
								double bf_ratio = (sq_mode_val / MODE_MIN);
								double log_bff = Math.log10(bf_ratio);
								double bf = (slope * log_bff);

								a.setB_factor(bf);

								index++;
							}

						f_index = String.format("%03d", (frame_index + 1));
						String output_file = OUT_INDIVIDUAL + "Mode_" + (mode + 1) + "_frame_" + f_index + ".pdb";
						PDB_IO.Write_PDB(output_file, shifted_atoms);

						frame_index++;
						shifted_atoms.clear();
						write_Pymol_Script_Individual_AA(mode);
					}
			}
	}

	private static void get_Essential_Visualization_AA()
	{
		int start = Input_Parameters_VIZ.MODE_START;
		int end = Input_Parameters_VIZ.MODE_END;

		Matrix sum = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, 0, 0);
		double eval_abs = Math.abs(eigenvalues.get(0));
		Matrix w_sum = sum.times(eval_abs);
		for (int i = 1; i < number_of_modes; i++)
			{
				Matrix plus = top_square_pca_modes.getMatrix(0, ROWS_Modes - 1, i, i);
				eval_abs = Math.abs(eigenvalues.get(i));
				Matrix w_plus = plus.times(eval_abs);
				w_sum = w_sum.plus(w_plus);
			}

		/* Norm and sort the vector */
		Matrix w_sum_normed = Projector.get_Normed_arrayF(w_sum);
		double[] sorted_w_sum_normed = w_sum_normed.getColumnPackedCopy();
		Arrays.sort(sorted_w_sum_normed);

		/* Establish the Log coloring scheme */
		double MODE_MIN = sorted_w_sum_normed[0];
		double MODE_MAX = sorted_w_sum_normed[ROWS_Modes - 1];
		if (MODE_MIN < Input_Parameters_VIZ.FLOOR) MODE_MIN = Input_Parameters_VIZ.FLOOR;
		double LOG_MODE_MIN = Math.log10(MODE_MIN);
		double LOG_MODE_MAX = Math.log10(MODE_MAX);
		double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
		double slope = (delta_y / delta_x);

		final double eigenvalue_max = eigenvalues.get(0); // The largest eigenvalue sets the wave number for the first mode

		int frame_index = 0;
		String f_index = "";
		for (int t = 0; t < Input_Parameters_VIZ.numberOfFramesEssential; t++) // FRAME LOOP:
			{
				List<Atom> shifted_atoms = new ArrayList<Atom>(atoms.size()); // create a new list to hold the shifted atoms
				for (Atom a : atoms)
					{
						Atom sa = new Atom(a);
						shifted_atoms.add(sa);
					}

				f_index = String.format("%03d", frame_index + 1);
				frame_index++;

				for (int mode = start - 1; mode < end - 1; mode++) // MODE LOOP: iterates over the modes
					{
						/* Get the MAX eigenvector component. This controls the maximum deflection... */
						Matrix evector = top_evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
						ArrayList<Double> evects = new ArrayList<Double>();
						for (int i = 0; i < ROWS_Evectors; i++)
							{
								double component = evector.get(i, 0);
								double compAbs = Math.abs(component);
								evects.add(compAbs);
							}
						Collections.sort(evects, Collections.reverseOrder());
						double Evector_Max = evects.get(0);

						/* Calculate the weight, this determines the amount of deflection of the atoms... */
						/* A Sine function ensures the first structure (t=0) is unperturbed; this preserves chain connectivity in PyMol movies... */
						double eval = eigenvalues.get(mode);
						double Amplitude = Math.sqrt(eval / eigenvalue_max);
						double omega = ((2 * Input_Parameters_VIZ.numberOfCyclesEssential * pi / Input_Parameters_VIZ.numberOfFramesEssential) * Math.sqrt(eigenvalue_max / eval));
						double weight = Amplitude * Math.sin(omega * t);

						/* The mode amplitude should be adjusted to maximally stretch the bonds without overly distorting them, about 1 Angstrom */
						double MA = Input_Parameters_VIZ.mode_amplitude;
						if (number_of_atoms < 30) MA = ((MA / Evector_Max) + (number_of_atoms / 100));
						if (number_of_atoms >= 30) MA = ((MA / Evector_Max) + (number_of_atoms / 50));

						int index = 0;
						for (Atom a : shifted_atoms) // ATOM LOOP: Iterates through the list of atoms, cumulatively shifting the atomic coordinates;
							{
								double x_coord = a.getX();
								double y_coord = a.getY();
								double z_coord = a.getZ();

								double v_x = evector.get(index, 0);
								double v_y = evector.get((index + number_of_atoms), 0);
								double v_z = evector.get((index + 2 * number_of_atoms), 0);

								double shift_x = (v_x * weight * MA);
								double shift_y = (v_y * weight * MA);
								double shift_z = (v_z * weight * MA);

								a.setX(x_coord + shift_x);
								a.setY(y_coord + shift_y);
								a.setZ(z_coord + shift_z);

								double bff = (w_sum_normed.get(index, 0));
								if (bff < MODE_MIN) bff = MODE_MIN;
								if (bff > MODE_MAX) bff = MODE_MAX;
								double BFF = (bff / MODE_MIN);
								double log_bff = Math.log10(BFF);
								double bf = (slope * log_bff);

								a.setB_factor(bf);

								index++;
							}
					}
				String output_file = OUT_ESSENTIAL + "Essential_Modes_" + number_of_modes + "_frame_" + f_index + ".pdb";
				PDB_IO.Write_PDB(output_file, shifted_atoms);
				write_Pymol_Script_Essential_AA(number_of_modes);
			}
	}

	private static void write_Pymol_Script_Individual_AA(int mode_number)
	{
		try
			{
				String out_dir_PS = "";
				if (OUT_INDIVIDUAL.endsWith(File.separator))
					{
						out_dir_PS = OUT_INDIVIDUAL.substring(0, OUT_INDIVIDUAL.length() - 1);
					}
				else
					out_dir_PS = OUT_INDIVIDUAL;

				File pymol_script_file = new File(OUT_INDIVIDUAL + "Mode_" + (mode_number + 1) + ".pml");
				BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));

				script_file_writer.write("from pymol import cmd" + "\n");
				script_file_writer.write("from pymol.cgo import *" + "\n");
				script_file_writer.write("bg_color white" + "\n");
				script_file_writer.write("from glob import glob" + "\n");
				script_file_writer.write("cd " + out_dir_PS + "\n");
				script_file_writer.write("filelist = glob (" + " \"" + "Mode_" + (mode_number + 1) + "_frame*.pdb\" )" + "\n");
				script_file_writer.write("for file in filelist: cmd.load( file, " + "\"" + "Mode_" + (mode_number + 1) + "\" )" + "\n");
				script_file_writer.write("show sticks, " + "Mode_" + (mode_number + 1) + "\n");
				script_file_writer.write("util.performance(0)" + "\n");
				script_file_writer.write("cmd.spectrum(\"b\",selection=\"(all)\",quiet=0)" + "\n");
				script_file_writer.write("set depth_cue,0" + "\n");
				script_file_writer.write("set two_sided_lighting,1" + "\n");
				script_file_writer.write("set stick_ball,1" + "\n");
				script_file_writer.write("set stick_ball_ratio,-1.000" + "\n");
				script_file_writer.write("set stick_h_scale,1.000" + "\n");
				script_file_writer.write("set movie_fps,10.000" + "\n");
				script_file_writer.write("rebuild" + "\n");
				script_file_writer.write("orient" + "\n");
				script_file_writer.close();

			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not write the PyMOL (TM) script file: " + OUT_INDIVIDUAL + "Mode_" + (mode_number + 1) + ".pml");
				io.printStackTrace();
			}
	}

	private static void write_Pymol_Script_Essential_AA(int modes)
	{
		try
			{
				String out_dir_PS = "";
				if (OUT_ESSENTIAL.endsWith(File.separator))
					{
						out_dir_PS = OUT_ESSENTIAL.substring(0, OUT_ESSENTIAL.length() - 1);
					}
				else
					out_dir_PS = OUT_ESSENTIAL;

				File pymol_script_file = new File(OUT_ESSENTIAL + "Essential_Modes_" + number_of_modes + ".pml");
				BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));

				script_file_writer.write("from pymol import cmd" + "\n");
				script_file_writer.write("from pymol.cgo import *" + "\n");
				script_file_writer.write("bg_color white" + "\n");
				script_file_writer.write("from glob import glob" + "\n");
				script_file_writer.write("cd " + out_dir_PS + "\n");
				script_file_writer.write("filelist = glob (" + " \"" + "Essential_Modes_" + number_of_modes + "_frame*.pdb\" )" + "\n");
				script_file_writer.write("for file in filelist: cmd.load( file, " + "\"" + "Essential_Modes_" + (modes) + "\" )" + "\n");
				script_file_writer.write("show sticks, " + "Essential_Modes_" + (modes) + "\n");
				script_file_writer.write("util.performance(0)" + "\n");
				script_file_writer.write("cmd.spectrum(\"b\",selection=\"(all)\",quiet=0)" + "\n");
				script_file_writer.write("set two_sided_lighting,1" + "\n");
				script_file_writer.write("set depth_cue,0" + "\n");
				script_file_writer.write("set stick_ball,1" + "\n");
				script_file_writer.write("set stick_ball_ratio,-1.000" + "\n");
				script_file_writer.write("set stick_h_scale,1.000" + "\n");
				script_file_writer.write("set movie_fps,10.000" + "\n");
				script_file_writer.write("rebuild" + "\n");
				script_file_writer.write("orient" + "\n");
				script_file_writer.close();

			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not write the PyMOL(TM) script file: " + OUT_ESSENTIAL + "Essential_Modes_" + number_of_modes + ".pml");
				io.printStackTrace();
			}
	}

	/* ---------------------------------------------------------------------------------------------------------------- */

	public static void main(String[] args) throws IOException
	{

		rm = RoundingMode.HALF_UP;
		nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setMinimumFractionDigits(3);
		nf.setRoundingMode(rm);

		System.out.println("Running JEDi VIZ Driver: ");
		System.out.println("Getting the input file: ");
		try
			{

				input_path = "VIZ_Parameters.txt";
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_path;

				if (args.length >= 1)
					{
						input_path = args[0];
						System.out.println("User specified input file (must be the first program argument)");
						System.out.println("These are the specified program args:");
						for (int i = 0; i < args.length; i++)
							{
								System.out.println("Arg " + (i + 1) + " Value = " + args[i]);
							}
						in_path = input_path;
					}
				System.out.println("Working Directory = " + WD);
				System.out.println("Input File Path = " + in_path);
				boolean check = new File(in_path).exists();
				if (!check)
					{
						System.err.println("The entered Input File does not exist: " + in_path);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				input_reader = new BufferedReader(new FileReader(in_path));

			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the input file: " + input_path);
				System.err.println("Program terminating.\n");
				e.printStackTrace();
				System.exit(0);
			}

		System.out.println("Reading Input File... ");
		read_input_file();
		System.out.println("Reading Data Files... ");
		read_data_files();
		Output_Control();
		System.out.println("Construdting PCA modes... ");
		construct_PCA_Modes();

		if (do_individual)
			{
				System.out.println("Performing Individual All-Atom Mode Visualizations... ");
				startTime = System.nanoTime();
				get_Individual_Mode_Visualizations_AA();
				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("\tDone. (" + nf.format(totalTime / 1000000000.0) + " seconds)");

			}

		if (do_essential)
			{
				System.out.println("Performing All-Atom Essential Mode Visualization... ");
				startTime = System.nanoTime();
				get_Essential_Visualization_AA();
				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("\tDone. (" + nf.format(totalTime / 1000000000.0) + " seconds)");

			}

		date = DateUtils.now();
		write_Job_Log();
		System.out.println("JEDi Mode Visualization successfully completed: " + date);
	}
}
