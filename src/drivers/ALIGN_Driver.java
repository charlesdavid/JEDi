package drivers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;

import Jama.Matrix;
import jedi.JEDi_Get_Transformed_Coordinates;
import supportIO.Input_Parameters_ALIGN;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportPlot.Plot_Line_Chart;
import supportPlot.RMSD_Plot;
import supportPlot.STATS_Plot;

public class ALIGN_Driver
{
	static boolean exist, success;
	static int reference_column, ROWS, COLS, number_Of_Input_Lines;
	static double RMSIP, MAD_Score, Z_Score;
	static String key, value, delim = "=", line, input_path, out_dir, description, date;
	static Matrix stats;
	static Matrix original_PDB_Coordinates, aligned_PDB_Coordinates, aligned_PDB_Coordinates_Outliers_REMOVED, aligned_PDB_Coordinates_Outliers_SELECTED, REF_COL,
			adjustments_per_variable_REMOVE_Outliers, adjustments_per_variable_SELECT_Outliers;
	static ArrayList<Double> atomic_RMSFs, conformation_RMSDs;
	static ArrayList<String> lines;
	static File log;
	static StringTokenizer sToken;
	static BufferedReader input_reader;
	static BufferedWriter log_writer;
	static Hashtable<String, String> parameters;
	static Input_Parameters_ALIGN ip;

	private static void read_input_file()
	{
		number_Of_Input_Lines = 0;
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
				ip = new Input_Parameters_ALIGN(parameters);
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
				log = new File(out_dir + "ALIGN_LOG.txt");
				log_writer = new BufferedWriter(new FileWriter(log));
				log_writer.write("JEDi Alignment Driver version 1.0" + "\n");
				log_writer.write("Description: " + Input_Parameters_ALIGN.DESCRIPTION + "\n");
				log_writer.write("Output directory: " + Input_Parameters_ALIGN.OUT_DIR + "\n");
				log_writer.write("PDB Coordinates file: " + Input_Parameters_ALIGN.PDB_COORDS + "\n");
				log_writer.write("Reference Column: " + Input_Parameters_ALIGN.REFERENCE_COL + "\n");
				log_writer.write("MAD_Score cutoff: " + Input_Parameters_ALIGN.MAD_Score + "\n");
				log_writer.write("Z_Score cutoff: " + Input_Parameters_ALIGN.Z_Score + "\n");
				log_writer.write("PROCESS_OUTLIERS = " + Input_Parameters_ALIGN.PROCESS_OUTLIERS + "\n");
				log_writer.write("COMPRESS Matrix files = " + Input_Parameters_ALIGN.COMPRESS + "\n");
				if (Input_Parameters_ALIGN.COMPRESS) log_writer.write("COMPRESS METHOD " + Input_Parameters_ALIGN.COMPRESS_METHOD + "\n");
				log_writer.write("Alignment to Reference Column Completed: " + date);
				log_writer.close();
			}

		catch (IOException e)
			{
				System.err.println("Could not write the LOG file: " + out_dir + log);
				System.err.println("Program terminating.\n");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void Output_Control()
	{
		description = Input_Parameters_ALIGN.DESCRIPTION;
		out_dir = Input_Parameters_ALIGN.OUT_DIR + "JEDi_ALIGN_Driver_" + description + File.separatorChar;
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
		System.out.println("Output Directory: " + out_dir);

		String path, path1, path2, path3;

		if (Input_Parameters_ALIGN.COMPRESS)
			{
				path1 = out_dir + "Aligned_PDB_Coordinates.txt." + Input_Parameters_ALIGN.COMPRESS_METHOD;
				path2 = out_dir + "Aligned_PDB_Coordinates_Outliers_REMOVED.txt." + Input_Parameters_ALIGN.COMPRESS_METHOD;
				path3 = out_dir + "Aligned_PDB_Coordinates_Outliers_SELECTED.txt." + Input_Parameters_ALIGN.COMPRESS_METHOD;
			}

		else
			{
				path1 = out_dir + "Aligned_PDB_Coordinates.txt";
				path2 = out_dir + "Aligned_PDB_Coordinates_Outliers_REMOVED.txt";
				path3 = out_dir + "Aligned_PDB_Coordinates_Outliers_SELECTED.txt";
			}

		if (Input_Parameters_ALIGN.PROCESS_OUTLIERS)
			{
				Matrix_IO.write_Matrix(aligned_PDB_Coordinates, path1, 12, 3);
				Matrix_IO.write_Matrix(aligned_PDB_Coordinates_Outliers_REMOVED, path2, 12, 3);
				Matrix_IO.write_Matrix(aligned_PDB_Coordinates_Outliers_SELECTED, path3, 12, 3);

				if (Input_Parameters_ALIGN.MAD_Score > 0)
					{
						path = out_dir + "Adjustments_per_Variable_REMOVE_Outliers_MAD.txt";
						Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_Outliers, path, 12, 0);
						Plot_Line_Chart.createChart_One_Series(out_dir, "Adjustments per Variable REMOVE Outliers MAD", "Variable Index", "Counts",
								adjustments_per_variable_REMOVE_Outliers);

						path = out_dir + "Adjustments_per_Variable_SELECT_Outliers_MAD.txt";
						Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_Outliers, path, 12, 0);
						Plot_Line_Chart.createChart_One_Series(out_dir, "Adjustments per Variable SELECT Outliers MAD", "Variable Index", "Counts",
								adjustments_per_variable_SELECT_Outliers);
					}
				if (Input_Parameters_ALIGN.Z_Score > 0)
					{
						path = out_dir + "Adjustments_per_Variable_REMOVE_Outliers_Z.txt";
						Matrix_IO.write_Matrix(adjustments_per_variable_REMOVE_Outliers, path, 12, 0);
						Plot_Line_Chart.createChart_One_Series(out_dir, "Adjustments per Variable REMOVE Outliers Z", "Variable Index", "Counts",
								adjustments_per_variable_REMOVE_Outliers);

						path = out_dir + "Adjustments_per_Variable_SELECT_Outliers_Z.txt";
						Matrix_IO.write_Matrix(adjustments_per_variable_SELECT_Outliers, path, 12, 0);
						Plot_Line_Chart.createChart_One_Series(out_dir, "Adjustments per Variable SELECT Outliers Z", "Variable Index", "Counts",
								adjustments_per_variable_SELECT_Outliers);
					}
			}
		else
			{
				Matrix_IO.write_Matrix(aligned_PDB_Coordinates, path1, 12, 3);
			}

		path = out_dir + "Atomic_RMSFs.txt";
		List_IO.write_Double_List(atomic_RMSFs, path, 12);

		path = out_dir + "Conformation_RMSDs.txt";
		List_IO.write_Double_List(conformation_RMSDs, path, 6);

		RMSD_Plot.create_RMSD_XY_Chart(out_dir, "Atomic_RMSFs", "Atom Number", atomic_RMSFs);
		RMSD_Plot.create_RMSD_XY_Chart(out_dir, "Conformation_RMSDs", "Conformation Number", conformation_RMSDs);

		STATS_Plot.create_Variables_Stat_XY_Chart(out_dir, "Variable_Means", "Atom Index", "Mean", stats, 1);
		STATS_Plot.create_Variables_Stat_XY_Chart(out_dir, "Variable_Variances", "Atom Index", "Variance", stats, 2);
		STATS_Plot.create_Variables_Stat_XY_Chart(out_dir, "Variable_Skews", "Atom Index", "Skew", stats, 3);
		STATS_Plot.create_Variables_Stat_XY_Chart(out_dir, "Variable_Kurtosis", "Atom Index", "Kurtosis", stats, 4);
	}

	public static void main(String[] args)
	{
		System.out.println("Running JEDi Alignment Driver:\n");
		System.out.println("Getting the input file: ");
		try
			{
				input_path = "ALIGN_Parameters.txt";
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_path;

				if (args.length >= 1)
					{
						input_path = args[0];
						System.out.println("\tUser specified input file (must be the first program argument)");
						System.out.println("\tThese are the specified program args:");
						for (int i = 0; i < args.length; i++)
							{
								System.out.println("\t\tArg " + (i + 1) + " Value = " + args[i]);
							}
						in_path = input_path;
					}
				// System.out.println("\tWorking Directory = " + WD);
				System.out.println("\tInput File Path = " + in_path);
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

		read_input_file();

		original_PDB_Coordinates = Matrix_IO.read_Matrix(Input_Parameters_ALIGN.PDB_COORDS);
		ROWS = original_PDB_Coordinates.getRowDimension();
		COLS = original_PDB_Coordinates.getColumnDimension();

		reference_column = Input_Parameters_ALIGN.REFERENCE_COL;
		if (reference_column > COLS)
			{
				System.err.println("Sorry, the reference column does not exist. Using the first column instead.");
				reference_column = 1;
			}
		REF_COL = original_PDB_Coordinates.getMatrix(0, ROWS - 1, reference_column - 1, reference_column - 1);

		System.out.println("Using Coordinate File: " + Input_Parameters_ALIGN.PDB_COORDS + "\n");
		System.out.println("\t" + ROWS + " Rows by " + COLS + " Columns\n");
		System.out.println("\tReference Column = " + reference_column + "\n");
		System.out.println("\tMAD_Score Cutoff = " + Input_Parameters_ALIGN.MAD_Score + "\n");
		System.out.println("\tZ_Score Cutoff = " + Input_Parameters_ALIGN.Z_Score + "\n");
		if (Input_Parameters_ALIGN.COMPRESS) System.out.println("\tOutput coordinate matrices will be compressed using " + Input_Parameters_ALIGN.COMPRESS_METHOD + "\n");
		System.out.println("\tAligning Coordinates...\n");

		JEDi_Get_Transformed_Coordinates tf_coords = new JEDi_Get_Transformed_Coordinates(original_PDB_Coordinates, REF_COL);
		tf_coords.set_MAD_Score_Cutoff(Input_Parameters_ALIGN.MAD_Score);
		tf_coords.set_Z_Score_Cutoff(Input_Parameters_ALIGN.Z_Score);

		aligned_PDB_Coordinates = tf_coords.get_SS_Transformed_coords();
		aligned_PDB_Coordinates_Outliers_REMOVED = tf_coords.get_SS_transformed_coordinates_OUTLIERS_REMOVED();
		aligned_PDB_Coordinates_Outliers_SELECTED = tf_coords.get_SS_transformed_coordinates_OUTLIERS_SELECTED();

		conformation_RMSDs = (ArrayList<Double>) tf_coords.get_SS_Conformation_RMSDs();
		atomic_RMSFs = (ArrayList<Double>) tf_coords.get_SS_RMSF();
		stats = tf_coords.get_SS_coordinate_STATS();

		if (Input_Parameters_ALIGN.MAD_Score > 0)
			{
				adjustments_per_variable_REMOVE_Outliers = tf_coords.getAdjustments_per_variable_REMOVE_MAD();
				adjustments_per_variable_SELECT_Outliers = tf_coords.getAdjustments_per_variable_SELECT_MAD();
			}
		if (Input_Parameters_ALIGN.Z_Score > 0)
			{
				adjustments_per_variable_REMOVE_Outliers = tf_coords.getAdjustments_per_variable_REMOVE_Z();
				adjustments_per_variable_SELECT_Outliers = tf_coords.getAdjustments_per_variable_SELECT_Z();
			}

		Output_Control();

		date = supportIO.DateUtils.now();
		write_Job_Log();

		System.out.println("Alignment Finished: " + date + "\n");
		System.out.println("------------------------------------------------------------------------------------------------------------\n");
	}
}
