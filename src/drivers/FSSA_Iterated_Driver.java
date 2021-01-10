package drivers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;

import Jama.Matrix;
import support.Subspace_Analysis2;
import supportIO.DateUtils;
import supportIO.Matrix_IO;
import supportIO.Test_Numeric_Type;

/**
 * JED class FSSA_Iterated_Driver: Driver program for the Iterated Fast Subspace Analysis class. Reads input file: "SSA.txt" Copyright (C) 2012 Dr. Charles David
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
 * 
 */
public class FSSA_Iterated_Driver
{

	static int number_Of_Input_Lines, line_count, num_of_jobs, job_number, ROWS1, COLS1, ROWS2, COLS2, iterations;
	static String line, input_path, directory1, directory2, name1, name2, out_dir, description, batch_description, date = DateUtils.now();
	static double RMSIP;
	static Matrix matrix1, matrix2, PAs, Avg_RMSIP_Score_Matrix, RMSIP_Std_Dev_Matrix;
	static ArrayList<Double> principle_angles_svd, RMSIPs;
	static ArrayList<String> lines;
	static File file_1, file_2, Job_log, Batch_log;
	static StringTokenizer sToken;
	static BufferedReader input_reader;
	static PrintWriter Batch_log_Writer, Job_log_writer;
	static NumberFormat nf;
	static boolean check, isDir, exist;

	private static void read_input_file()
	{
		number_Of_Input_Lines = 0;
		lines = new ArrayList<>();
		System.out.println("Below is the input file that was read: " + input_path);
		System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
		try
			{
				while ((line = input_reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						System.out.println(line);
						number_Of_Input_Lines++;
					}

				input_reader.close();
				System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
				System.out.println("The number of lines of parameters in the input file is: " + number_Of_Input_Lines + "\n");
				if (number_Of_Input_Lines < 7)
					{
						System.err.println("INSUFFICIENT DATA IN THE INPUT FILE:");
						System.err.println("THERE MUST BE 7 LINES OF PARAMETERS FOR ONE COMPARISON.");
						System.err.println("Terminating program execution.");
						System.exit(0);
					}

			}
		catch (IOException e)
			{
				System.err.println("IOException thrown. Could not read the input file. Program will terminate.\n");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void read_batch_parameters()
	{
		line_count = 0;
		System.out.println("Reading line " + (line_count + 1)); // Reads line 1, the number of jobs
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		String test = sToken.nextToken();
		boolean OK = Test_Numeric_Type.test_Integer(test);
		if (OK) num_of_jobs = Integer.parseInt(test);
		System.out.println("\tThe number of jobs =  " + num_of_jobs);
		if (!OK || num_of_jobs < 1)
			{
				System.err.println("Expected Number of Jobs to be a positive integer, but got: " + test);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1)); // Reads the Output Directory for the batch
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		out_dir = sToken.nextToken();
		System.out.println("\tOutput Directory = " + out_dir);
		if (!(out_dir.endsWith("/") || out_dir.endsWith("\\")))
			{
				System.err.println("Expected the Output Directory to end with / or \\\\, but got: " + line);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		exist = new File(out_dir).exists();
		if (!exist)
			{
				System.out.println("\tThe output directory does not exist.");
				System.out.println("\tAttempting to create it:");
				boolean success = (new File(out_dir)).mkdirs();
				if (success) System.out.println("\t\tSuccess.");
				if (!success)
					{
						System.out.println("Failed to create the output directory:  " + out_dir);
						System.exit(0);
					}
			}
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1)); // Reads the Batch Description
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		batch_description = sToken.nextToken();
		System.out.println("\tBatch Description = " + batch_description);
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1));
		line = lines.get(line_count); // Reads the divider line between the batch parameters and the first job
		System.out.println(line + "\n");
		line_count++;
		/* *********************************************************************************************************************************** */
	}

	private static void read_job_parameters()
	{
		if (line_count >= lines.size())
			{
				System.err.println("No more data in input file for remaining jobs in batch.");
				System.err.println("User specified too many jobs.");
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1) + ", the FIRST LINE of job: " + (job_number + 1));
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		description = sToken.nextToken();
		System.out.println("\tJob Description = " + description);
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1) + ", the SECOND LINE of job: " + (job_number + 1));
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		directory1 = sToken.nextToken();
		System.out.println("\tDirectory 1 = " + directory1);
		if (!(directory1.endsWith("/") || directory1.endsWith("\\")))
			{
				System.err.println("Expected the working Directory to end with / or \\\\, but got: " + directory1);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		isDir = new File(directory1).isDirectory();
		if (!isDir)
			{
				System.err.println("The entered directory does not exist as a directory: " + directory1);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		name1 = sToken.nextToken();
		System.out.println("\tFirst Eigenvector File = " + name1);
		exist = new File(directory1 + name1).exists();
		if (!exist)
			{
				System.err.println("The entered Eigenvector File does not exist: " + directory1 + name1);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1) + ", the THIRD LINE of job: " + (job_number + 1));
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		directory2 = sToken.nextToken();
		System.out.println("\tDirectory 2 = " + directory2);
		if (!(directory2.endsWith("/") || directory2.endsWith("\\")))
			{
				System.err.println("Expected the working Directory to end with / or \\\\, but got: " + directory2);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		isDir = new File(directory2).isDirectory();
		if (!isDir)
			{
				System.err.println("The entered directory does not exist as a directory: " + directory2);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		name2 = sToken.nextToken();
		System.out.println("\tFirst Eigenvector File = " + name2);
		exist = new File(directory2 + name2).exists();
		if (!exist)
			{
				System.err.println("The entered Eigenvector File does not exist: " + directory2 + name2);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
		line_count++;
		/* ************************************************************************************************************************************************** */
		System.out.println("Reading line " + (line_count + 1));
		line = lines.get(line_count); // Reads the divider line between jobs
		System.out.println(line + "\n");
		line_count++;
		/* ************************************************************************************************************************************************** */
	}

	private static void initialize_Batch_Log()
	{
		Batch_log = new File(out_dir + "FSSA_Iterated_Batch_Log_" + batch_description + ".txt");
		try
			{
				Batch_log_Writer = new PrintWriter(new BufferedWriter(new FileWriter(Batch_log)));

			}
		catch (IOException e)
			{
				System.err.println("Could not create the Batch Log file. Program terminating");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void initialize_Job_Log()
	{
		Job_log = new File(out_dir + "FSSA_Iterated_Job_Log_" + description + "_dim_" + COLS1 + ".txt");
		try
			{
				Job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
			}
		catch (IOException e)
			{
				System.err.println("Could not create the Job Log file for Job Number " + job_number + 1 + ". Program terminating");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void write_Logs()
	{
		Batch_log_Writer.write(description + "\t" + "Max_SS_dim = " + COLS1 + ", VS_dim = " + ROWS1 + "\n");
		Batch_log_Writer.write("RMSIPs:\n");
		Job_log_writer.write("\nThe RMSIPs for the iterated subspace comparison are:\n");

		nf.setMaximumIntegerDigits(1);
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);

		int dim = 1;
		for (double p : RMSIPs)
			{
				Job_log_writer.write(String.format("%-6s%-6s%-6s%n", "Dim", dim, nf.format(p)));
				Batch_log_Writer.write(String.format("%-6s%-6s%-6s%n", "Dim", dim, nf.format(p)));
				dim++;
			}
		Job_log_writer.write("\nThe Principle Angle Spectra for the subspaces (in degrees) are:\n");
		PAs.print(Job_log_writer, 6, 0);

		Job_log_writer.write("\nAnalysis completed: " + date);
		Job_log_writer.close();

		Batch_log_Writer.write("\nThe Principle Angle Spectra for the subspaces (in degrees) are:\n");
		PAs.print(Batch_log_Writer, 6, 0);
		Batch_log_Writer.write("\n");
		Batch_log_Writer.flush();
	}

	public static void main(String[] args)
	{
		nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);

		System.out.println("Running FSSA Iterated Driver: ");
		System.out.println("Getting the input file: ");
		try
			{

				input_path = "SSA.txt";
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_path;

				if (args.length >= 1)
					{
						input_path = args[0];
						System.out.println("The path to the input file must be the first argument:");
						System.out.println("These are the command args:");
						for (int i = 0; i < args.length; i++)
							{
								System.out.println("Arg " + (i + 1) + " Value = " + args[i]);
							}
						in_path = input_path;
					}
				System.out.println("Working Directory = " + WD);
				System.out.println("Input File Path = " + in_path);
				check = new File(in_path).exists();
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
		read_batch_parameters();
		initialize_Batch_Log();

		for (job_number = 0; job_number < (num_of_jobs); job_number++)
			{
				read_job_parameters();

				file_1 = new File(directory1 + name1);
				file_2 = new File(directory2 + name2);
				matrix1 = Matrix_IO.read_Matrix(directory1, name1);
				ROWS1 = matrix1.getRowDimension();
				COLS1 = matrix1.getColumnDimension();
				matrix2 = Matrix_IO.read_Matrix(directory2, name2);
				ROWS2 = matrix2.getRowDimension();
				COLS2 = matrix2.getColumnDimension();

				if (ROWS1 != ROWS2)
					{
						System.err.println("FATAL ERROR: The subspaces do not come from the same vector space. Program will terminate.");
						System.exit(0);
					}
				if (COLS1 != COLS2)
					{
						System.err.println("FATAL ERROR: The subspaces do not have the same dimensions. Program will terminate");
						System.exit(0);
					}

				Subspace_Analysis2 ssa = new Subspace_Analysis2(matrix1, matrix2);
				ssa.get_fast_SSA_iterated();
				ssa.get_random_SSA();

				RMSIPs = ssa.getRMSIPs();
				PAs = ssa.getPAs();

				Avg_RMSIP_Score_Matrix = ssa.getAvg_RMSIP_Score_Matrix();
				RMSIP_Std_Dev_Matrix = ssa.getRMSIP_Std_Dev_Matrix();
				iterations = ssa.getIterations();

				initialize_Job_Log();
				write_Logs();
				RMSIPs.clear();
			}
		Batch_log_Writer.close();
	}
}
