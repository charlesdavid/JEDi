package jedi;

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

/**
 * JED class SSA Driver: Driver program for the Subspace Analysis class. Reads input file: "SSA.txt" Copyright (C) 2012 Dr. Charles David
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

public class SSA_Driver
	{

		static int number_Of_Input_Lines, line_count, num_of_jobs, job_number, ROWS1, COLS1, ROWS2, COLS2;
		static double RMSIP, max_angle;
		static String line, input_path, path1, path2, out_dir, description, batch_description, date = DateUtils.now();
		static Matrix matrix1, matrix2, projections, CO_matrix_1_2, CO_matrix_2_1;
		static ArrayList<Double> cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products, vectorial_sum_of_angles, RMSIPs;
		static ArrayList<String> lines;
		static File Job_log, Batch_log;
		static BufferedReader input_reader;
		static StringTokenizer sToken;
		static PrintWriter Batch_log_Writer, Job_log_writer;
		static NumberFormat nf, nff;
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

					} catch (IOException e)
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
				if (OK)
					num_of_jobs = Integer.parseInt(test);
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
				if (!(out_dir.endsWith(File.separator)))
					{
						System.err.println("Expected the Output Directory to end with: " + File.separatorChar + ", but got: " + out_dir);
						System.err.println("Attempting to fix...");
						out_dir = out_dir + File.separator;
					}
				exist = new File(out_dir).exists();
				if (!exist)
					{
						System.out.println("\tThe output directory does not exist.");
						System.out.println("\tAttempting to create it:");
						boolean success = (new File(out_dir)).mkdirs();
						if (success)
							System.out.println("\t\tSuccess.");
						if (!success)
							{
								System.err.println("Failed to create the output directory:  " + out_dir);
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
						System.err.println("User may have specified too many jobs.");
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
				path1 = sToken.nextToken();
				System.out.println("\tEigenvector File 1 = " + path1);
				exist = new File(path1).exists();
				if (!exist)
					{
						System.err.println("The entered Eigenvector File does not exist: " + path1);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				line_count++;
				/* ************************************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1) + ", the THIRD LINE of job: " + (job_number + 1));
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				path2 = sToken.nextToken();
				System.out.println("\tEigenvector File 2 = " + path2);
				exist = new File(path2).exists();
				if (!exist)
					{
						System.err.println("The entered Eigenvector File does not exist: " + path2);
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
				Batch_log = new File(out_dir + "SSA_Batch_Log_" + batch_description + ".txt");
				try
					{
						Batch_log_Writer = new PrintWriter(new BufferedWriter(new FileWriter(Batch_log)));

					} catch (IOException e)
					{
						System.err.println("Could not create the Batch Log file. Program terminating");
						e.printStackTrace();
						System.exit(0);
					}
			}

		private static void initialize_Job_Log()
			{
				Job_log = new File(out_dir + "SSA_Job_Log_" + description + "_dim_" + COLS1 + ".txt");
				try
					{
						Job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
					} catch (IOException e)
					{
						System.err.println("Could not create the Job Log file for Job Number " + job_number + 1 + ". Program terminating");
						e.printStackTrace();
						System.exit(0);
					}
			}

		private static void write_Logs()
			{
				Job_log_writer.write("Matrix 1: " + path1 + "\n");
				Job_log_writer.write("Rows: " + ROWS1 + "\n");
				Job_log_writer.write("Cols: " + COLS1 + "\n");
				Job_log_writer.write("Matrix 2: " + path2 + "\n");
				Job_log_writer.write("Rows: " + ROWS2 + "\n");
				Job_log_writer.write("Cols: " + COLS2 + "\n\n");
				Job_log_writer.write("The projections of each vector in subspace 1 with each vector in subspace 2 are:\n");
				projections.print(Job_log_writer, 6, 3);
				Job_log_writer.write("The cumulative overlaps CO_k for each vector in subspace 1 with all the vectors in subspace 2 are:\n");
				Batch_log_Writer.write(String.format("%-25s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\t", description, "SS_dim", COLS1, "VS_dim", ROWS1, "RMSIP", nf.format(RMSIP), "PAs"));
				int j = 1;
				for (double d : cumulative_overlaps_1_2)
					{
						Job_log_writer.write(String.format("%-8s%-4s%-8s%n", "Vector", j, nf.format(d)));
						j++;
					}
				Job_log_writer.write("\nThe cumulative overlaps CO_k for each vector in subspace 2 with all the vectors in subspace 1 are:\n");
				j = 1;
				for (double d : cumulative_overlaps_2_1)
					{
						Job_log_writer.write(String.format("%-8s%-4s%-8s%n", "Vector", j, nf.format(d)));
						j++;
					}
				Job_log_writer.write("\nThe RMSIP for the two subspaces is " + nf.format(RMSIP) + "\n");
				Job_log_writer.write("\nThe principle angles between these spaces (in degrees) are:\n");
				int i = 1;
				for (double p : principle_angles_svd)
					{
						Job_log_writer.write(String.format("%-8s%-4s%-8s%n", "Angle", i, nf.format(p)));
						Batch_log_Writer.write(String.format("%-4s\t", nff.format(p)));
						i++;
					}
				Job_log_writer.write("\nThe products of cosines (CPs) for these two subspaces (equivalents in degrees) are:\n");
				i = 1;
				for (double c : cosine_products)
					{
						Job_log_writer.write(String.format("%-4s%-4s%-8s%n", "CP", i, nf.format(Math.acos(c) * 180 / Math.PI)));
						i++;
					}

				Job_log_writer.write("\nThe vectorial sums of angles (VS, in degrees) are:\n");
				i = 1;
				for (double p : vectorial_sum_of_angles)
					{
						Job_log_writer.write(String.format("%-8s%-4s%-8s%n", "VS", i, nf.format(p)));
						i++;
					}
				Job_log_writer.write("\nMaximum possible angle between the two subspaces = " + nf.format(max_angle) + " degrees.\n\n");
				Job_log_writer.write("Analysis completed: " + date);
				Job_log_writer.close();
				Batch_log_Writer.write(String.format("%-12s%-6s%n", "Max_Angle ", nff.format(max_angle)));
				Batch_log_Writer.flush();

			}

		public static void main(String[] args)
			{

				nf = NumberFormat.getInstance();
				nf.setMaximumFractionDigits(2);
				nf.setMinimumFractionDigits(2);
				nff = NumberFormat.getInstance();
				nff.setMaximumIntegerDigits(3);
				nff.setMaximumFractionDigits(0);
				nff.setMinimumFractionDigits(0);

				System.out.println("Running SSA Driver: ");
				System.out.println("Getting the SSA input file: ");
				try
					{

						input_path = "SSA.txt";
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
						check = new File(in_path).exists();
						if (!check)
							{
								System.err.println("The entered Input File does not exist: " + in_path);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
						input_reader = new BufferedReader(new FileReader(in_path));

					} catch (FileNotFoundException e)
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

						matrix1 = Matrix_IO.read_Matrix(path1);
						ROWS1 = matrix1.getRowDimension();
						COLS1 = matrix1.getColumnDimension();
						matrix2 = Matrix_IO.read_Matrix(path2);
						ROWS2 = matrix2.getRowDimension();
						COLS2 = matrix2.getColumnDimension();

						if (ROWS1 != ROWS2)
							{
								System.err.println("ERROR: The subspaces do not come from the same vector space. Program will terminate.");
								System.exit(0);
							}
						if (COLS1 != COLS2)
							{
								System.err.println("ERROR: The subspaces do not have the same dimensions. Program will terminate");
								System.exit(0);
							}

						cumulative_overlaps_1_2 = new ArrayList<>();
						cumulative_overlaps_2_1 = new ArrayList<>();
						principle_angles_svd = new ArrayList<>();
						cosine_products = new ArrayList<>();
						vectorial_sum_of_angles = new ArrayList<>();

						Subspace_Analysis ssa = new Subspace_Analysis(matrix1, matrix2);
						ssa.get_SSA();

						RMSIP = ssa.getRMSIP();
						max_angle = ssa.getMax_angle();
						principle_angles_svd = ssa.getPrinciple_angles_svd();
						cumulative_overlaps_1_2 = ssa.getCumulative_overlaps_1_2();
						cumulative_overlaps_2_1 = ssa.getCumulative_overlaps_2_1();
						cosine_products = ssa.getCosine_products();
						CO_matrix_1_2 = ssa.getCO_matrix_1_2();
						CO_matrix_2_1 = ssa.getCO_matrix_2_1();
						projections = ssa.getProjections();
						vectorial_sum_of_angles = ssa.getVectorial_sum_of_angles();

						initialize_Job_Log();
						write_Logs();

						cumulative_overlaps_1_2.clear();
						cumulative_overlaps_2_1.clear();
						cosine_products.clear();
						principle_angles_svd.clear();
						vectorial_sum_of_angles.clear();

						ssa = null;
						System.gc();
					}
				Batch_log_Writer.close();
			}
	}
