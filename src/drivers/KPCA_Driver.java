package drivers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;

import Jama.Matrix;
import jedi.JEDi_Do_Kernel_PCA;
import supportIO.DateUtils;
import supportIO.Matrix_IO;
import supportPlot.Plot_XY_Scatter;

/**
 * JEDi class KPCA_Driver: Driver program for performing kernel PCA analysis.
 * 
 * Input file is "KPCA.txt" Copyright (C) 2018 Dr. Charles David
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>.
 *
 * @author Dr. Charles David
 */
public class KPCA_Driver
{
	static String line, key, value, PROJECTIONS, OUT_DIRECTORY, DESCRIPTION, file_name_head, path, date;
	static String input_file = "KPCA_Parameters.txt", delim = "=";
	static int number_Of_Input_Lines, line_count, ROWS, COLS, NUMBER_KPCs, MAX_KERNEL_FRAMES, KDE_CELL, KDE_RESOLUTION, MULTIPLIER;
	static double SLOPE, SIGMA, KDE_MARGIN, KERNEL_SHRINKAGE, FLOOR;
	static long startTime, endTime, totalTime;
	static Matrix delta_vector_projections;
	static ArrayList<String> lines;
	static File Job_log;
	static StringTokenizer sToken;
	static BufferedReader input_reader;
	static PrintWriter Job_log_writer;
	static NumberFormat nf;
	static DateUtils now;
	static Hashtable<String, String> parameters;
	static Plot_XY_Scatter plot;

	static boolean check, isDir, exist, success;

	static boolean do_linear_kernel;
	static boolean do_Degree_2_Poly;
	static boolean do_Degree_3_Poly;
	static boolean do_Degree_4_Poly;
	static boolean do_XY_Poly;
	static boolean do_Euclidean;
	static boolean do_Mahalanobis;
	static boolean do_Gaussian;
	static boolean do_Cauchy;
	static boolean do_Circular;
	static boolean do_Log;
	static boolean do_MI;
	static boolean do_MI_KDE;
	static boolean do_Neural_Net;


	private static void read_input_file()
	{
		number_Of_Input_Lines = 0;
		lines = new ArrayList<String>();
		parameters = new Hashtable<>();

		System.out.println("Reading the input file: " + input_file);
		System.out.println("Below are the parameters to be processed:");
		System.out.println("------------------------------------------------------------------------------------------------------------------------------------------");
		try
			{
				while ((line = input_reader.readLine()) != null && line.length() >= 1)
					{
						if (!(line.startsWith("#")))
							{
								lines.add(line);
								System.out.println(line);
								sToken = new StringTokenizer(line);
								key = sToken.nextToken(delim);
								value = sToken.nextToken(delim);
								parameters.put(key, value);
								number_Of_Input_Lines++;
							}
					}

				input_reader.close();
				System.out.println("-----------------------------------------------------------------------------------------------------------------------------------");
				System.out.println("The number of KEY=VALUE Pairs in the input file is: " + number_Of_Input_Lines + "\n");
			}
		catch (IOException e)
			{
				System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
				e.printStackTrace();
				System.exit(0);
			}
		System.gc();
	}

	private static void assign_parameters()
	{
		DESCRIPTION = parameters.get("DESCRIPTION");
		PROJECTIONS = parameters.get("PROJECTIONS");
		NUMBER_KPCs = Integer.valueOf(parameters.get("NUMBER_KPCs"));

		MAX_KERNEL_FRAMES = Integer.valueOf(parameters.get("MAX_KERNEL_FRAMES"));
		KDE_CELL = Integer.valueOf(parameters.get("KDE_CELL"));
		KDE_RESOLUTION = Integer.valueOf(parameters.get("KDE_RESOLUTION"));
		MULTIPLIER = Integer.valueOf(parameters.get("MULTIPLIER"));

		if (parameters.get("do_linear_kernel").equals("true")) do_linear_kernel = true;
		if (parameters.get("do_Degree_2_Poly").equals("true")) do_Degree_2_Poly = true;
		if (parameters.get("do_Degree_3_Poly").equals("true")) do_Degree_3_Poly = true;
		if (parameters.get("do_Degree_4_Poly").equals("true")) do_Degree_4_Poly = true;
		if (parameters.get("do_XY_Poly").equals("true")) do_XY_Poly = true;
		if (parameters.get("do_Cauchy").equals("true")) do_Cauchy = true;
		if (parameters.get("do_Circular").equals("true")) do_Circular = true;
		if (parameters.get("do_Euclidean").equals("true")) do_Euclidean = true;
		if (parameters.get("do_Gaussian").equals("true")) do_Gaussian = true;
		if (parameters.get("do_Log").equals("true")) do_Log = true;
		if (parameters.get("do_Mahalanobis").equals("true")) do_Mahalanobis = true;
		if (parameters.get("do_Neural_Net").equals("true")) do_Neural_Net = true;
		if (parameters.get("do_MI").equals("true")) do_MI = true;
		if (parameters.get("do_MI_KDE").equals("true")) do_MI_KDE = true;

		SIGMA = Double.valueOf(parameters.get("SIGMA"));
		SLOPE = Double.valueOf(parameters.get("SLOPE"));

		KDE_MARGIN = Double.valueOf(parameters.get("KDE_MARGIN"));
		KERNEL_SHRINKAGE = Double.valueOf(parameters.get("KERNEL_SHRINKAGE"));
		FLOOR = Double.valueOf(parameters.get("FLOOR"));

		OUT_DIRECTORY = parameters.get("OUT_DIRECTORY");

		exist = new File(OUT_DIRECTORY).exists();
		if (!exist)
			{
				success = (new File(OUT_DIRECTORY)).mkdirs();
				if (!success)
					{
						System.err.println("Sorry, unable to create the output directory. Please check permissions and format: " + OUT_DIRECTORY);
						System.exit(0);
					}
			}
		/* **************************************************************************************************************************************** */

		if (!(OUT_DIRECTORY.endsWith(File.separator)))
			{
				System.err.println("Expected the directory to end with " + File.separator + ", but got: " + OUT_DIRECTORY);
				System.err.println("Attempting to fix...");
				OUT_DIRECTORY = OUT_DIRECTORY + File.separator;
			}
		isDir = new File(OUT_DIRECTORY).isDirectory();
		if (!isDir)
			{
				System.err.println("Sorry, but the entered directory is not recognized as a proper directory: " + OUT_DIRECTORY);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		/* **************************************************************************************************************************************** */
		boolean check = new File(PROJECTIONS).exists();
		if (!check)
			{
				System.err.println("Sorry, the entered PDB Reference File can not be found: " + PROJECTIONS);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		/* **************************************************************************************************************************************** */
		delta_vector_projections = Matrix_IO.read_Matrix(PROJECTIONS);
	}

	private static void initialize_Job_Log()
	{
		Job_log = new File(OUT_DIRECTORY + "KPCA_Job_Log_" + DESCRIPTION + ".txt");
		try
			{
				Job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));

			}
		catch (IOException e)
			{
				System.err.println("Could not create the Job Log file: " + "KPCA_Job_Log_" + DESCRIPTION + ".txt");
				System.err.println("Program terminating");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void write_Log()
	{
		Job_log_writer.write("Principal components: " + PROJECTIONS + "\n");
		Job_log_writer.write("Number of Conformations in the DVPs: " + COLS + "\n");
		Job_log_writer.write("Number of Modes available in the DVPs: " + ROWS + "\n");
		Job_log_writer.write("The number of KPCA modes to process is: " + (NUMBER_KPCs) + "\n");
		Job_log_writer.flush();
		date = DateUtils.now();
		Job_log_writer.write("\nAnalysis completed: " + date);
		Job_log_writer.close();
	}

	public static void get_kPCA()
	{
		JEDi_Do_Kernel_PCA kpca = new JEDi_Do_Kernel_PCA(delta_vector_projections, OUT_DIRECTORY, NUMBER_KPCs, DESCRIPTION, MAX_KERNEL_FRAMES);

		kpca.set_Sigma(SIGMA);
		kpca.set_Slope(SLOPE);

		kpca.set_MI_kernelResolution(KDE_RESOLUTION);
		kpca.setMI_Kernel_Cell(KDE_CELL);
		kpca.setShrinkage(KERNEL_SHRINKAGE);
		kpca.setMultiplier(MULTIPLIER);
		kpca.setFLOOR(FLOOR);

		kpca.setDo_linear_kernel(do_linear_kernel);
		kpca.setDo_Degree_2_Poly(do_Degree_2_Poly);
		kpca.setDo_Degree_3_Poly(do_Degree_3_Poly);
		kpca.setDo_Degree_4_Poly(do_Degree_4_Poly);
		kpca.setDo_XY_Poly(do_XY_Poly);
		kpca.setDo_XY_Poly(do_Cauchy);
		kpca.setDo_Circular(do_Circular);
		kpca.setDo_Euclidean(do_Euclidean);
		kpca.setDo_Gaussian(do_Gaussian);
		kpca.setDo_Log(do_Log);
		kpca.setDo_Mahalanobis(do_Mahalanobis);
		kpca.setDo_Sigmoid(do_Neural_Net);
		kpca.setDo_MI(do_MI);
		kpca.setDo_MI_KDE(do_MI_KDE);
		kpca.setDo_XY_Poly(do_XY_Poly);

		kpca.kPCA_Driver();
	}

	/* ********************************************************************************************************************************** */
	public static void main(String[] args)
	{
		nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setMinimumFractionDigits(3);
		nf.setRoundingMode(RoundingMode.HALF_UP);

		System.out.println("Running the KPCA Driver: ");
		System.out.println("Getting the input file: ");
		try
			{
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_file;

				if (args.length >= 1)
					{
						input_file = args[0];
						System.out.println("User specified input file (must be the first program argument)");
						System.out.println("These are the specified program args:");
						for (int i = 0; i < args.length; i++)
							{
								System.out.println("Arg " + (i + 1) + " Value = " + args[i]);
							}
						in_path = input_file;
					}
				System.out.println("Working Directory = " + WD);
				System.out.println("Input File Path = " + in_path);
				check = new File(in_path).exists();
				if (!check)
					{
						System.err.println("Sorry, but the specified Input File could not be found: " + in_path);
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

		System.out.println("Reading Input File... ");
		read_input_file();

		System.out.println("Assigning parameters... ");
		assign_parameters();

		initialize_Job_Log();

		System.out.println("Doing KPCA... ");

		get_kPCA();

		write_Log();
	}
}