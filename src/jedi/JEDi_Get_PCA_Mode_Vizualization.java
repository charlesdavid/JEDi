package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import Jama.Matrix;
import support.Atom;
import support.Projector;
import supportIO.FortranFormat;
import supportIO.Input_Parameters;
import supportIO.PDB_File_Parser;
import supportIO.PDB_IO;

/**
 * JED class JED_PCA_Mode_Visualization: Constructs sets of PDB files and a PyMol(TM) Script to animate the PCA modes derived from the Cartesian subset. Copyright (C) 2012 Dr.
 * Charles David
 *
 * @author Dr.Charles David
 */

/**
 * @author Charles
 *
 */
public class JEDi_Get_PCA_Mode_Vizualization
{
	boolean exist, success;
	double mode_amplitude;

	String out_dir, file_name_head, output_file;
	BufferedWriter output_file_writer;

	final double LOG_FLOOR, delta_y = 99, pi = Math.PI, MA, modeNorm = 50.000, modeNormRes = 100.000;

	int modes_viz, numberEssentialModes = Input_Parameters.numberModeComponents;
	final int number_of_atoms, ROWS_Evectors, ROWS_Modes, COLS;

	final int numberEssentialCycles, numberEssentialFrames, numberModeFrames, numberModeCycles;
	final String directory, description, modelPCA, typePCA, CA = "CA", C = "C", N = "N", O = "O";
	final double[] pca_mode_maxs, pca_mode_mins;
	final List<Double> eigenvalues;
	final Matrix evectors, square_pca_modes;
	final List<Atom> reference_atoms;
	final FortranFormat formatter;
	final PDB_File_Parser parser;
	final DecimalFormat df;
	final RoundingMode rm;

	// ************************************* CONSTRUCTOR ******************************************************************************

	/**
	 *
	 * Constructor for generating the sets of PDB files for visualizing the Cartesian PCA modes.
	 *
	 * @param dir           The working directory
	 * @param des           The job description
	 * @param residues      The list of original residues from the PDB file
	 * @param residueAtoms  The atoms that comprise the subset
	 * @param evects        The eigenvectors for the top cPCA modes
	 * @param sq_modes      The top cPCA modes
	 * @param maxs          The array containing the maximum component of each cPCA mode
	 * @param mins          The array containing the minimum component of each cPCA mode
	 * @param mode_amp      The mode amplitude, which controls how far the CA, C, N, and O atoms are displaced
	 * @param num_modes_viz The number of cPCA modes for which to generate structures
	 * @param pca_type      The PCA type: displacement, Cartesian,...
	 * @param pca_model     The PCA model: COV, CORR, or PCORR
	 */
	public JEDi_Get_PCA_Mode_Vizualization(List<Atom> residueAtoms, List<Double> evals, Matrix evects, Matrix sq_modes, double[] maxs, double[] mins, int num_modes_viz,
			String pca_type, String pca_model)
	{
		this.directory = Input_Parameters.DIRECTORY;
		this.description = Input_Parameters.DESCRIPTION;
		this.MA = Input_Parameters.VIZ_MODE_SCALE_FACTOR;
		this.LOG_FLOOR = Input_Parameters.LOG_FLOOR;
		this.numberEssentialCycles = Input_Parameters.numberEssentialCycles;
		this.numberEssentialFrames = Input_Parameters.numberEssentialFrames;
		this.numberModeFrames = Input_Parameters.numberModeFrames;
		this.numberModeCycles = Input_Parameters.numberModeCycles;

		this.reference_atoms = residueAtoms;
		this.eigenvalues = evals;
		this.evectors = evects;
		this.square_pca_modes = sq_modes;
		this.pca_mode_maxs = maxs;
		this.pca_mode_mins = mins;
		this.modes_viz = num_modes_viz;
		this.typePCA = pca_type;
		this.modelPCA = pca_model;

		this.ROWS_Evectors = evects.getRowDimension();
		this.ROWS_Modes = square_pca_modes.getRowDimension();
		this.number_of_atoms = reference_atoms.size();
		this.COLS = evects.getColumnDimension();

		if (modes_viz > COLS) modes_viz = COLS;
		if (numberEssentialModes > COLS) numberEssentialModes = COLS;

		this.formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		this.formatter.setAddReturn(true);
		this.parser = new PDB_File_Parser();

		this.rm = RoundingMode.HALF_UP;
		this.df = new DecimalFormat("#0.000");
		this.df.setRoundingMode(rm);
	}

	// **************************************************** METHODS *****************************************************************************************************

	public void get_Mode_Visualizations_BackBone()
	{
		// System.out.println("\tDoing Backbone PCA Mode Visualization: \n");

		for (int mode = 0; mode < modes_viz; mode++) // MODE LOOP
			{
				double MODE_MIN = pca_mode_mins[mode];
				double MODE_MAX = pca_mode_maxs[mode];
				if (MODE_MIN < LOG_FLOOR) MODE_MIN = LOG_FLOOR;
				double LOG_MODE_MIN = Math.log10(MODE_MIN);
				double LOG_MODE_MAX = Math.log10(MODE_MAX);
				double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
				double slope = (delta_y / delta_x);

				Matrix evector = evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
				Matrix sq_mode = square_pca_modes.getMatrix(0, ROWS_Modes - 1, mode, mode);

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

				/* The mode amplitude should be set to stretch the bonds as much as possible without breaking them, about 1 - 2 Angstroms */
				if (number_of_atoms < 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNormRes));
				if (number_of_atoms >= 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNorm));

				// System.out.println("\tMode " + (mode + 1) + " Mode Amplitude: " + df.format(mode_amplitude));

				String f_index = "";
				int frame_index = 0;

				for (double mc = 0; mc < (2 * Math.PI); mc += (Math.PI / 10)) // FRAME LOOP: Generates 20 frames
					{
						List<Atom> shifted_atoms = new ArrayList<Atom>(reference_atoms.size()); // create a new list to hold the shifted atoms
						for (Atom a : reference_atoms)
							{
								Atom sa = new Atom(a);
								shifted_atoms.add(sa);
							}

						/* Set the mode weight: */
						/* A Sine function ensures the first structure (mc=0) is unperturbed; preserves chain connectivity in PyMol movies... */
						double w = Math.sin(mc);

						int atom_count = (ROWS_Evectors / 3);
						int index = 0;
						int count = 0;
						for (int i = 0; i < shifted_atoms.size(); i++) // ATOM LOOP: Shifts all backbone atoms along the *Alpha Carbon Eigenvector*.
							{
								index = (count / 4); // there are 4 atoms per residue in the list of backbone atoms (N-CA-C-O)
								Atom a = shifted_atoms.get(i);

								double x_coord = a.getX();
								double y_coord = a.getY();
								double z_coord = a.getZ();

								double v_x = evector.get(index, 0);
								double v_y = evector.get((index + atom_count), 0);
								double v_z = evector.get((index + 2 * atom_count), 0);

								double shift_x = (v_x * w * mode_amplitude);
								double shift_y = (v_y * w * mode_amplitude);
								double shift_z = (v_z * w * mode_amplitude);

								a.setX(x_coord + shift_x);
								a.setY(y_coord + shift_y);
								a.setZ(z_coord + shift_z);

								double bff = (sq_mode.get(index, 0));
								if (bff < MODE_MIN) bff = MODE_MIN;
								if (bff > MODE_MAX) bff = MODE_MAX;
								double BFF = bff / MODE_MIN;
								double log_bff = Math.log10(BFF);
								double bf = (slope * log_bff);

								a.setB_factor(bf);

								count++;
							}

						f_index = String.format("%03d", (frame_index + 1));
						output_file = file_name_head + "_Mode_" + (mode + 1) + "_" + modelPCA + "_frame_" + f_index + ".pdb.bz2";
						PDB_IO.Write_BZ2_PDB(output_file, shifted_atoms);
						frame_index++;
						shifted_atoms.clear();
					}
				write_Pymol_Script_AA(mode);
			}
	}

	public void get_Mode_Visualizations_All_Atom()
	{
		// System.out.println("\tDoing All Atom PCA Mode Visualization: \n");

		for (int mode = 0; mode < modes_viz; mode++) // MODE LOOP
			{
				Matrix evector = evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
				Matrix sq_mode = square_pca_modes.getMatrix(0, ROWS_Modes - 1, mode, mode);

				/* Establish the log coloring scheme using logs, Square Mode components as B-factors */
				double MODE_MIN = pca_mode_mins[mode];
				double MODE_MAX = pca_mode_maxs[mode];
				if (MODE_MIN < LOG_FLOOR) MODE_MIN = LOG_FLOOR;
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

				/* The mode amplitude should be adjusted to maximally stretch bonds without overly distorting them, about 1 Angstrom */
				if (number_of_atoms < 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNormRes));
				if (number_of_atoms >= 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNorm));

				// System.out.println("\tMode " + (mode + 1) + " Mode Amplitude: " + df.format(mode_amplitude));

				String f_index = "";
				int frame_index = 0;

				for (int t = 0; t < numberModeFrames; t++) // FRAME LOOP:
					{
						// System.out.println("Computing Frame: " + (frame_index + 1));

						double omega = ((2 * numberModeCycles * pi / numberModeFrames));
						/* A Sine function ensures the first structure is unperturbed; Preserves chain connectivity in PyMol movies... */
						double weight = Math.sin(omega * t);

						List<Atom> shifted_atoms = new ArrayList<Atom>(reference_atoms.size()); // create a new list to hold the shifted atoms
						for (Atom a : reference_atoms)
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

								double shift_x = (v_x * weight * mode_amplitude);
								double shift_y = (v_y * weight * mode_amplitude);
								double shift_z = (v_z * weight * mode_amplitude);

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
						output_file = file_name_head + "_Mode_" + (mode + 1) + "_" + modelPCA + "_frame_" + f_index + ".pdb.bz2";
						PDB_IO.Write_BZ2_PDB(output_file, shifted_atoms);

						// output_file_writer = new BufferedWriter(new FileWriter(output_file));
						// parser.write_PDB(output_file_writer, shifted_atoms, formatter);
						// output_file_writer.close();

						frame_index++;
						shifted_atoms.clear();
						write_Pymol_Script_AA(mode);
					}
			}
	}

	/* ---------------------------------------------------------------------------------------------------------------- */

	public void get_Essential_Visualization_BackBone()
	{
		// System.out.println("\tDoing Backbone Essential Mode Visualization: \n");

		/* Get the top combined square mode */
		Matrix sum = square_pca_modes.getMatrix(0, ROWS_Modes - 1, 0, 0);
		double eval_abs = Math.abs(eigenvalues.get(0));
		Matrix w_sum = sum.times(eval_abs);
		for (int i = 1; i < modes_viz; i++)
			{
				Matrix plus = square_pca_modes.getMatrix(0, ROWS_Modes - 1, i, i);
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
		if (MODE_MIN < LOG_FLOOR) MODE_MIN = LOG_FLOOR;
		double LOG_MODE_MIN = Math.log10(MODE_MIN);
		double LOG_MODE_MAX = Math.log10(MODE_MAX);
		double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
		double slope = (delta_y / delta_x);

		final double eigenvalue_max = eigenvalues.get(0); // The largest eigenvalue sets the wave number for the first mode
		if (modes_viz < numberEssentialModes) numberEssentialModes = modes_viz;
		int frame_index = 0;
		String f_index = "";
		for (int t = 0; t < numberEssentialFrames; t++) // FRAME LOOP
			{
				List<Atom> shifted_atoms = new ArrayList<Atom>(reference_atoms.size()); // create a new list to hold the shifted atoms
				for (Atom a : reference_atoms)
					{
						Atom sa = new Atom(a);
						shifted_atoms.add(sa);
					}

				f_index = String.format("%03d", frame_index + 1);
				frame_index++;
				double omega = 0, weight = 0;

				for (int mode = 0; mode < numberEssentialModes; mode++) // MODE LOOP
					{
						/* Get the MAX eigenvector component. This controls the maximum deflection... */
						Matrix evector = evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
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
						omega = ((2 * numberEssentialCycles * pi / numberEssentialFrames) * Math.sqrt(eigenvalue_max / eval));
						weight = Amplitude * Math.sin(omega * t);

						/* The mode amplitude should be set to stretch the bonds as much as possible without overly distorting them, about 1 Angstrom */
						if (number_of_atoms < 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNormRes));
						if (number_of_atoms >= 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNorm));

						// System.out.println("\n\tEssential Modes " + (mode + 1) + " Mode Amplitude: " + df.format(mode_amplitude));

						int atom_count = (ROWS_Evectors / 3);
						int index = 0, count = 0;
						for (int i = 0; i < shifted_atoms.size(); i++) // ATOM LOOP: Shifts all backbone atoms along the *Alpha Carbon Eigenvector*.
							{
								index = (count / 4); // there are 4 atom symbols in the vector of backbone atoms
								Atom a = shifted_atoms.get(i);

								double x_coord = a.getX();
								double y_coord = a.getY();
								double z_coord = a.getZ();

								double v_x = evector.get(index, 0);
								double v_y = evector.get((index + atom_count), 0);
								double v_z = evector.get((index + 2 * atom_count), 0);

								double shift_x = (v_x * weight * mode_amplitude);
								double shift_y = (v_y * weight * mode_amplitude);
								double shift_z = (v_z * weight * mode_amplitude);

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

								count++;
							}
					}

				output_file = file_name_head + "_Essential_Modes_" + modes_viz + "_" + modelPCA + "_frame_" + f_index + ".pdb.bz2";
				PDB_IO.Write_BZ2_PDB(output_file, shifted_atoms);
				write_Pymol_Script_Essential_AA(modes_viz);
			}
	}

	public void get_Essential_Visualization_All_Atom()
	{
		// System.out.println("\t Doing All Atom Essential Mode Visualization: \n");

		/* Get the top combined square mode */
		Matrix sum = square_pca_modes.getMatrix(0, ROWS_Modes - 1, 0, 0);
		double eval_abs = Math.abs(eigenvalues.get(0));
		Matrix w_sum = sum.times(eval_abs);
		for (int i = 1; i < modes_viz; i++)
			{
				Matrix plus = square_pca_modes.getMatrix(0, ROWS_Modes - 1, i, i);
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
		if (MODE_MIN < LOG_FLOOR) MODE_MIN = LOG_FLOOR;
		double LOG_MODE_MIN = Math.log10(MODE_MIN);
		double LOG_MODE_MAX = Math.log10(MODE_MAX);
		double delta_x = (LOG_MODE_MAX - LOG_MODE_MIN);
		double slope = (delta_y / delta_x);

		final double eigenvalue_max = eigenvalues.get(0); // The largest eigenvalue sets the wave number for the first mode
		if (modes_viz < numberEssentialModes) numberEssentialModes = modes_viz;

		int frame_index = 0;
		String f_index = "";
		for (int t = 0; t < numberEssentialFrames; t++) // FRAME LOOP:
			{
				List<Atom> shifted_atoms = new ArrayList<Atom>(reference_atoms.size()); // create a new list to hold the shifted atoms
				for (Atom a : reference_atoms)
					{
						Atom sa = new Atom(a);
						shifted_atoms.add(sa);
					}

				f_index = String.format("%03d", frame_index + 1);
				frame_index++;

				for (int mode = 0; mode < numberEssentialModes; mode++) // MODE LOOP: iterates over the modes
					{
						/* Get the MAX eigenvector component. This controls the maximum deflection... */
						Matrix evector = evectors.getMatrix(0, ROWS_Evectors - 1, mode, mode);
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
						double omega = ((2 * numberEssentialCycles * pi / numberEssentialFrames) * Math.sqrt(eigenvalue_max / eval));
						double weight = Amplitude * Math.sin(omega * t);

						/* The mode amplitude should be adjusted to maximally stretch the bonds without overly distorting them, about 1 Angstrom */
						if (number_of_atoms < 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNormRes));
						if (number_of_atoms >= 30) mode_amplitude = ((MA / Evector_Max) + (number_of_atoms / modeNorm));

						int index = 0;
						for (Atom a : shifted_atoms) // ATOM LOOP: Iterates through the list of atoms, cumulatively shifting the atomic coordinates;
							{
								double x_coord = a.getX();
								double y_coord = a.getY();
								double z_coord = a.getZ();

								double v_x = evector.get(index, 0);
								double v_y = evector.get((index + number_of_atoms), 0);
								double v_z = evector.get((index + 2 * number_of_atoms), 0);

								double shift_x = (v_x * weight * mode_amplitude);
								double shift_y = (v_y * weight * mode_amplitude);
								double shift_z = (v_z * weight * mode_amplitude);

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

				output_file = file_name_head + "_Essential_Modes_" + modes_viz + "_" + modelPCA + "_frame_" + f_index + ".pdb.bz2";
				PDB_IO.Write_BZ2_PDB(output_file, shifted_atoms);
				write_Pymol_Script_Essential_AA(modes_viz);
			}
	}

	/* ---------------------------------------------------------------------------------------------------------------- */

	private void write_Pymol_Script_AA(int mode_number)
	{
		try
			{
				String out_dir_PS = "";
				if (out_dir.endsWith(File.separator))
					{
						out_dir_PS = out_dir.substring(0, out_dir.length() - 1);
					}
				else
					out_dir_PS = out_dir;

				String name = "ss_" + number_of_atoms;
				File pymol_script_file = new File(file_name_head + "_Mode_" + (mode_number + 1) + "_" + modelPCA + ".pml");
				BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
				script_file_writer.write("from pymol import cmd" + "\n");
				script_file_writer.write("from pymol.cgo import *" + "\n");
				script_file_writer.write("bg_color white" + "\n");
				script_file_writer.write("from glob import glob" + "\n");
				script_file_writer.write("cd " + out_dir_PS + "\n");
				script_file_writer.write("filelist = glob (" + " \"" + name + "_Mode_" + (mode_number + 1) + "_" + modelPCA + "_frame*.pdb.bz2\" )" + "\n");
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
				System.err.println("IOException thrown. Could not write the PyMOL (TM) script file: " + file_name_head + "_Mode_" + (mode_number + 1) + "_" + modelPCA + ".pml");
				io.printStackTrace();
			}
	}

	private void write_Pymol_Script_Essential_AA(int modes)
	{
		String out_dir_PS = "";
		if (out_dir.endsWith(File.separator))
			{
				out_dir_PS = out_dir.substring(0, out_dir.length() - 1);
			}
		else
			out_dir_PS = out_dir;

		String name = "ss_" + number_of_atoms;
		String filename = file_name_head + "_Essential_Modes_" + modes_viz + "_" + modelPCA + ".pml";
		try
			{
				File pymol_script_file = new File(filename);
				BufferedWriter script_file_writer = new BufferedWriter(new FileWriter(pymol_script_file));
				script_file_writer.write("from pymol import cmd" + "\n");
				script_file_writer.write("from pymol.cgo import *" + "\n");
				script_file_writer.write("bg_color white" + "\n");
				script_file_writer.write("from glob import glob" + "\n");
				script_file_writer.write("cd " + out_dir_PS + "\n");
				script_file_writer.write("filelist = glob (" + " \"" + name + "_Essential_Modes_" + modes_viz + "_" + modelPCA + "_frame*.pdb.bz2\" )" + "\n");
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
				System.err.println("IOException thrown. Could not write the PyMOL(TM) script file: " + filename);
				io.printStackTrace();
			}
	}

	/* ****************************************************** SETTERS ************************************************ */

	public void set_Output_Directory(String dir)
	{
		this.out_dir = dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();

		this.file_name_head = out_dir + "ss_" + number_of_atoms;
	}

	/* ***************************************************** GETTERS ************************************************ */

	double get_Mode_Amplitude()
	{
		return mode_amplitude;
	}
}
