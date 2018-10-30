package jedi;

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
import java.util.Arrays;
import java.util.StringTokenizer;

import Jama.Matrix;

/**
 * JED class FES_Driver: Driver program for getting the free energy of two order parameters using 2D KDE This class imports bits/kde developed by Philip DeCamp.
 * 
 * Input file is "FES.txt" Copyright (C) 2012 Dr. Charles David
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
public class FES_Driver
	{
		static String line, input_path, batch_description, directory, out_directory, description, principal_components, file_name_head, path, date;
		static int number_Of_Input_Lines, line_count, num_of_jobs, job_number, ROWS, COLS;
		static Integer OP1, OP2, number_of_points, OP_offset, KDE_offset;
		static double cellsize;
		static long startTime, endTime, totalTime;
		static double[] KDE_bounds, KDE_bandwidths;
		static Matrix FE, delta_vector_series, order_param1, order_param2;
		static ArrayList<String> lines;
		static File dvp_file, Job_log, Batch_log;
		static StringTokenizer sToken;
		static BufferedReader input_reader;
		static PrintWriter Batch_log_Writer, Job_log_writer;
		static NumberFormat nf;
		static boolean check, isDir, exist, success;
		private static KernelDensityEstimate2d KDE;
		static DateUtils now;

		/**
		 * Class to create outputs for generating Free Energy Surfaces based on two order parameters.
		 * 
		 * @param       DVPs: The DVPs (PCs) from a JED PCA run
		 * @param       op1: The first order parameter (DVP1) (INT: the first column number)
		 * @param       op2: The Second order parameter (DVP2) (INT: the second column number)
		 * @param       num_points: The number of points (The number of data points that are to be extracted from the DVP to use for the KDE)
		 * @param       offset: The offset for reading points from the DVP (important if there are multiple runs combined)
		 * @param       size: The cellsize for the 2D grid in the KDE. Set to zero to determine automatically.
		 * @param       dir: The job input directory.
		 * @param       out_dir: The output directory.
		 * @param descr
		 */
		public FES_Driver(String DVPs, int op1, int op2, int num_points, int offset, double size, String dir, String out_dir, String descr)
			{
				principal_components = DVPs;
				OP1 = op1;
				OP2 = op2;
				number_of_points = num_points;
				OP_offset = offset;
				cellsize = size;
				directory = dir;
				description = descr;

				delta_vector_series = Matrix_IO.read_Matrix(dir, principal_components);
				ROWS = delta_vector_series.getRowDimension();
				COLS = delta_vector_series.getColumnDimension();

				int row_index1 = OP_offset;
				int row_index2 = OP_offset + number_of_points - 1;

				if (row_index2 > ROWS - 1)
					row_index2 = ROWS - 1;

				order_param1 = delta_vector_series.getMatrix(row_index1, row_index2, OP1, OP1);
				order_param2 = delta_vector_series.getMatrix(row_index1, row_index2, OP2, OP2);

				nf = NumberFormat.getInstance();
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
				nf.setRoundingMode(RoundingMode.HALF_UP);
				out_directory = out_dir;
				file_name_head = out_dir + description + "_FES_" + (OP1 + 1) + "_" + (OP2 + 1) + ".txt";

			}

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
				/* ******************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1)); // Reads line 2, the Output Directory for the batch
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				out_directory = sToken.nextToken();
				System.out.println("\tOutput Directory = " + out_directory);
				if (!(out_directory.endsWith(File.separator)))
					{
						System.err.println("Expected the Output Directory to end with " + File.separator + ", but got: " + out_directory);
						System.err.println("Attempting to fix...");
						out_directory = out_directory + File.separator;
					}
				exist = new File(out_directory).exists();
				if (!exist)
					{
						System.out.println("\tThe output directory does not exist.");
						System.out.println("\tAttempting to create it:");
						boolean success = (new File(out_directory)).mkdirs();
						if (success)
							System.out.println("\t\tSuccess.");
						if (!success)
							{
								System.err.println("Failed to create the output directory:  " + out_directory);
								System.exit(0);
							}
					}
				line_count++;
				/* *********************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1)); // Reads line 3, the Batch Description
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				batch_description = sToken.nextToken();
				System.out.println("\tBatch Description = " + batch_description);
				line_count++;
				/* *********************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1)); // Reads the divider line between the batch parameters and the first job
				line = lines.get(line_count);
				System.out.println(line + "\n");
				line_count++;
				/* *********************************************************************************************************************************** */
			}

		private static void read_job_parameters()
			{
				int index = 1;
				if (line_count >= lines.size())
					{
						System.err.println("No more data in input file for remaining jobs in batch.");
						System.err.println("User specified too many jobs.");
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				/* *********************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1) + ",  LINE " + index + " of job: " + (job_number + 1));
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				description = sToken.nextToken();
				System.out.println("\tJob Description = " + description);
				line_count++;
				index++;
				/* ************************************************************************************************************************************ */
				System.out.println("Reading line " + (line_count + 1) + ",  LINE " + index + " of job: " + (job_number + 1));
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				directory = sToken.nextToken();
				System.out.println("\tJob Directory = " + directory);
				if (!(directory.endsWith(File.separator)))
					{
						System.err.println("Expected the Directory to end with " + File.separator + ", but got: " + directory);
						System.err.println("Attempting to fix...");
						directory = directory + File.separator;
					}
				isDir = new File(directory).isDirectory();
				if (!isDir)
					{
						System.err.println("The entered directory does not exist as a directory: " + directory);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				line_count++;
				index++;
				/* ************************************************************************************************************************************ */
				System.out.println("Reading line " + (line_count + 1) + ",  LINE " + index + " of job: " + (job_number + 1));
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				principal_components = sToken.nextToken();
				System.out.println("\tThe principal components = " + principal_components);
				exist = new File(directory + principal_components).exists();
				if (!exist)
					{
						System.err.println("The entered principal components file does not exist: " + directory + principal_components);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				line_count++;
				index++;
				/* *************************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1) + ",  LINE " + index + " of job: " + (job_number + 1));
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				OP1 = Integer.parseInt(sToken.nextToken()) - 1; // to adjust for Java numbering starting with 0
				System.out.println("\tOP1 = " + OP1);
				OP2 = Integer.parseInt(sToken.nextToken()) - 1; // to adjust for Java numbering starting with 0
				System.out.println("\tOP2 = " + OP2);
				number_of_points = Integer.parseInt(sToken.nextToken());
				System.out.println("\tNumber of points to use = " + number_of_points);
				OP_offset = Integer.parseInt(sToken.nextToken());
				System.out.println("\tOffset = " + OP_offset);
				cellsize = Double.parseDouble(sToken.nextToken());
				System.out.println("\tCellsize = " + cellsize);
				if (cellsize < 0)
					cellsize = 0;
				if (cellsize == 0)
					System.out.println("Program will automatically calculate a cellsize.");
				if (OP1 == null || OP1 < 0 || OP1 > COLS)
					OP1 = 0; // default value
				if (OP2 == null || OP2 < 0 || OP2 > COLS)
					OP2 = 1; // default value
				if (number_of_points == null || number_of_points <= 0)
					number_of_points = ROWS; // use all the points
				if (OP_offset == null || OP_offset < 0)
					OP_offset = 0;
				line_count++;
				/* *************************************************************************************************************************************** */
				System.out.println("Reading line " + (line_count + 1)); // Reads the divider line between jobs
				line = lines.get(line_count);
				System.out.println(line + "\n");
				line_count++;
				/* *************************************************************************************************************************************** */
			}

		private static void initialize_Batch_Log()
			{
				Batch_log = new File(out_directory + "FES_Batch_Log_" + batch_description + ".txt");
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
				Job_log = new File(out_directory + "FES_Job_Log_" + description + ".txt");
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
				Job_log_writer.write("Principal components: " + directory + principal_components + "\n");
				Job_log_writer.write("Total Number of Conformations (rows) used to create the DVPs: " + ROWS + "\n");
				Job_log_writer.write("Total Number of Modes (cols) used to create the DVPs: " + COLS + "\n");
				Job_log_writer.write("Order Parameter 1 (Col# 1): " + (OP1 + 1) + "\n");
				Job_log_writer.write("Order Parameter 2 (Col# 2): " + (OP2 + 1) + "\n");
				Job_log_writer.write("Number of points to extract from the OPs to use for KDE: " + number_of_points + "\n");
				Job_log_writer.write("Offset in selecting points from the order parameters: " + OP_offset + "\n");
				Job_log_writer.write("The cellsize used for the 2D grid for KDE: " + nf.format(cellsize) + "\n");
				Job_log_writer.write("The bounds for the KDE are: OP1(" + nf.format(KDE_bounds[0]) + "," + nf.format(KDE_bounds[2]) + "); OP2(" + nf.format(KDE_bounds[1]) + ","
						+ nf.format(KDE_bounds[3]) + ")" + "\n");
				Job_log_writer.write("The kernel bandwidths were: " + nf.format(KDE_bandwidths[0]) + "   " + nf.format(KDE_bandwidths[3]) + "\n\n");
				Job_log_writer.write(String.format("%-12s%-12s%-12s%-12s", " ", "OP1", "OP2", "FE"));
				Job_log_writer.flush();
				FE.print(Job_log_writer, 12, 3);
				Batch_log_Writer.write(String.format("%-1s%-16s%-8s%-8s%-8s%-8s%-8s%-12s%-12s%-12s%-8s%-4s", description, "Conformations =", order_param1.getRowDimension(),
						"OP1 =", OP1 + 1, "OP2 =", OP2 + 19, "Bounds are: OP1(",
						nf.format(KDE_bounds[0]) + "," + nf.format(KDE_bounds[2]) + "); OP2(" + nf.format(KDE_bounds[1]) + "," + nf.format(KDE_bounds[3]), ") ",
						"KDE Bandwidths are: ", nf.format(KDE_bandwidths[0])) + ", " + nf.format(KDE_bandwidths[3]) + "\n");
				date = DateUtils.now();
				Job_log_writer.write("\nAnalysis completed: " + date);
				Job_log_writer.close();
				Batch_log_Writer.flush();
			}

		/**
		 * Method to calculate the delta_G free energy from two order parameters (DVPx,DVPy)
		 */
		public static void get_FES()
			{
				/* Get the DVPs to use as order parameters in the deltaG free energy calculations */
				if (OP_offset + number_of_points > ROWS)
					System.err.println("ERROR! The OP offset plus the number of points exceeds the length of the PCs!");
				double[] order_parameter_1 = order_param1.getColumnPackedCopy();
				double[] order_parameter_2 = order_param2.getColumnPackedCopy();
				double[] order_parameter_1_sorted = order_param1.getColumnPackedCopy();
				double[] order_parameter_2_sorted = order_param2.getColumnPackedCopy();

				/* 2D KDE using Gaussian Functions */
				int length = order_parameter_1.length;
				double[] kde_array = new double[2 * length];
				for (int i = 0; i < length; i++)
					{
						kde_array[i + i] = order_parameter_1[i];
						kde_array[i + i + 1] = order_parameter_2[i];
					}
				Arrays.sort(order_parameter_1_sorted);
				Arrays.sort(order_parameter_2_sorted);
				double op1max = order_parameter_1_sorted[length - 1];
				double op2max = order_parameter_2_sorted[length - 1];
				double op1min = order_parameter_1_sorted[0];
				double op2min = order_parameter_2_sorted[0];
				KDE_bounds = new double[] { op1min, op2min, op1max, op2max };
				System.out.println("The bounds for the KDE are: OP1(" + KDE_bounds[0] + "," + KDE_bounds[2] + "); OP2(" + KDE_bounds[1] + "," + KDE_bounds[3] + ")");

				/* Define the KDE offset and number of points to use in calculating the KDE */
				KDE_offset = 0;
				number_of_points = length;
				if (cellsize == 0)
					{
						double maxDim = Math.max(KDE_bounds[2] - KDE_bounds[0], KDE_bounds[3] - KDE_bounds[1]);
						cellsize = (maxDim / 256);
						System.out.println("The effective cellsize is: " + cellsize);
					}
				try
					{
						KDE = KernelDensityEstimate2d.compute(kde_array, KDE_offset, number_of_points, KDE_bounds, cellsize, null);
					} catch (Exception e)
					{
						e.printStackTrace();
					}
				KDE_bandwidths = DiagonalBandwidthSelector2d.get_Bandwidths();
				System.out.println("Kernel Bandwidths are: " + nf.format(KDE_bandwidths[0]) + "\t" + nf.format(KDE_bandwidths[3]));
				double[] probabilities = new double[length];
				double[] probabilities_sorted = new double[length];
				for (int i = 0; i < length; i++)
					{
						double prob = KDE.apply(order_parameter_1[i], order_parameter_2[i]);
						probabilities[i] = prob;
						probabilities_sorted[i] = prob;
					}
				Arrays.sort(probabilities_sorted);
				final double prob_max = probabilities_sorted[length - 1];
				final double ln_prob_max = Math.log(prob_max);
				final double KBT = (-0.600); // Units are in kcal/mol, T = 300K (room temp)
				FE = new Matrix(length, 3);
				for (int i = 0; i < length; i++)
					{
						double prob = probabilities[i];
						double ln_prob = Math.log(prob);
						double delta_G = KBT * (ln_prob - ln_prob_max);
						if (delta_G <= 0)
							delta_G = 0.000; // solves the -0.0000000000000 problem...
						FE.set(i, 0, order_parameter_1[i]);
						FE.set(i, 1, order_parameter_2[i]);
						FE.set(i, 2, delta_G);
					}
				path = file_name_head;
				Matrix_IO.write_Matrix(FE, path, 12, 3);
			}

		public KernelDensityEstimate2d getKDE()
			{
				return KDE;
			}

		/* ********************************************************************************************************************************** */
		@SuppressWarnings("unused")
		public static void main(String[] args)
			{
				nf = NumberFormat.getInstance();
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
				nf.setRoundingMode(RoundingMode.HALF_UP);

				System.out.println("Running FES Driver: ");
				System.out.println("Getting the input file: ");
				try
					{

						input_path = "FES.txt";
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
								System.err.println("The entered Input File could not be found: " + in_path);
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
						initialize_Job_Log();
						FES_Driver fes = new FES_Driver(principal_components, OP1, OP2, number_of_points, OP_offset, cellsize, directory, out_directory, description);
						FES_Driver.get_FES();
						write_Logs();
						delta_vector_series = null;
						fes = null;
						System.gc();
					}
				Batch_log_Writer.close();
			}
	}