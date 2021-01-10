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
import java.util.List;
import java.util.StringTokenizer;

import Jama.Matrix;
import support.systemInfo;
import supportIO.List_IO;
import supportIO.Matrix_IO;

/**
 * JEDi class Pool_Driver: Constructs augmented matrices from input Coordinate Matrices:
 * 
 * Allows for multiple trajectories to be combined.
 * 
 * Options for which frames to use and to do down sampling are available.
 * 
 * The number of rows should be the same for all input matrices.
 * 
 * With no down sampling, the number of columns in the output matrix is equal to the sum of all the input matrix columns.
 * 
 * With down sampling, the number of columns in the output matrix is equal to the sum of all the input matrix columns divided by the Down Sample Factor.
 * 
 * Parameters are specified in an input file called "POOL.txt"
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */
public class POOL_Driver
{
	static String line, description, input_path, out_dir, path;
	static int line_count, number_of_input_matrices, ROWS, COLS, Col_sum, downsampleFactor;
	static double[][] first_array;
	static List<String> lines;
	static ArrayList<String> paths;
	static ArrayList<Integer> start_Indices;
	static ArrayList<Integer> end_Indices;
	static Matrix augmented_matrix;
	static StringTokenizer sToken;
	static File Job_Log;
	static BufferedReader input_reader;
	static PrintWriter Job_Log_Writer;
	static boolean doDownSample, OK, check, isDir, exist, compress;
	static NumberFormat nf3;
	static RoundingMode rm;

	private static void read_input_file()
	{
		lines = List_IO.read_Lines_From_File(input_path);

		System.out.println("Below is the input file that was read: " + input_path);
		System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
		for (String s : lines)
			{
				System.out.println(s);
			}
		System.out.println("--------------------------------------------------------------------------------------------------------------------------------");
	}

	private static void assign_parameters()
	{
		line_count = 0;
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		number_of_input_matrices = Integer.parseInt(sToken.nextToken());
		if (sToken.nextToken().equals("true")) compress = true;
		System.out.println("compress = " + compress);
		line_count++;
		/* *********************************************************************************************************************************** */

		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		description = sToken.nextToken();
		if (sToken.nextToken().equals("true")) doDownSample = true;
		System.out.println("doDownSample = " + doDownSample);
		if (doDownSample) downsampleFactor = Integer.parseInt(sToken.nextToken());
		System.out.println("downsampleFactor = " + downsampleFactor);
		line_count++;
		/* *********************************************************************************************************************************** */
		line = lines.get(line_count);
		sToken = new StringTokenizer(line);
		out_dir = sToken.nextToken();
		if (!(out_dir.endsWith("/") || out_dir.endsWith("\\")))
			{
				System.err.println("Expected the Output Directory to end with the File.Separator.Character");
				System.err.println("Attempting to fix...");
				out_dir = out_dir + File.separatorChar;
			}
		exist = new File(out_dir).exists();
		if (!exist)
			{
				System.err.println("\tThe output directory does not exist.");
				System.err.println("\tAttempting to create it:");
				boolean success = (new File(out_dir)).mkdirs();
				if (success) System.err.println("\t\tSuccess.");
				if (!success)
					{
						System.err.println("Failed to create the output directory:  " + out_dir);
						System.exit(0);
					}
			}
		line_count++;
		/* *********************************************************************************************************************************** */
		paths = new ArrayList<String>();
		start_Indices = new ArrayList<Integer>();
		end_Indices = new ArrayList<Integer>();

		for (int i = 0; i < number_of_input_matrices; i++)
			{
				line = lines.get(line_count);
				sToken = new StringTokenizer(line);
				path = sToken.nextToken();
				paths.add(path);
				int start = Integer.parseInt(sToken.nextToken());
				start_Indices.add(start);
				int end = Integer.parseInt(sToken.nextToken());
				end_Indices.add(end);

				System.out.println("\tMatrix " + (i + 1) + " = " + path);
				System.out.println("\tUsing colums from: " + start + " to " + end);
				exist = new File(path).exists();
				if (!exist)
					{
						System.err.println("The entered file can not be found: " + path);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				line_count++;
			}
	}

	private static void initialize_Log()
	{
		if (doDownSample) Job_Log = new File(out_dir + "POOL_Log_" + description + "_" + number_of_input_matrices + "_DSF_" + downsampleFactor + ".txt");
		else
			Job_Log = new File(out_dir + "POOL_Log_" + description + "_" + number_of_input_matrices + ".txt");
		try
			{
				Job_Log_Writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_Log)));
			}
		catch (IOException e)
			{
				System.err.println("Could not create the Job Log file. Check file/directory permissions. Program terminating");
				e.printStackTrace();
				System.exit(0);
			}
	}

	private static void write_Log()
	{
		Job_Log_Writer.write("The job description is: " + description + "\n");
		Job_Log_Writer.write("The number of matrices that were combined is: " + number_of_input_matrices + "\n");
		Job_Log_Writer.write("The output directory is: " + out_dir + "\n");
		Job_Log_Writer.write("The number of ROWS: " + ROWS + "\n");
		if (doDownSample) Job_Log_Writer.write("The input matrices were downsampled. DownSample Factor = " + downsampleFactor + "\n");
		Job_Log_Writer.write("The total number of columns (frames) in output matrix: " + Col_sum + "\n");
		if (compress) Job_Log_Writer.write("The BZIP2 Compressed Pooled Coordinate Matrix was written to: " + path);
		else
			Job_Log_Writer.write("The uncompressed Pooled Coordinate Matrix was written to: " + path);
		Job_Log_Writer.flush();
		Job_Log_Writer.close();
	}

	private static int[] do_DownSample(int start, int end, int dsf)
	{
		int cols = (end - start + 1);
		int size = (cols / downsampleFactor);
		System.out.println("Number of columns after down sampling = " + size);
		int[] column_indices = new int[size];
		for (int i = 0; i < size; i++)
			{
				column_indices[i] = (i * downsampleFactor + start);
			}
		return column_indices;
	}

	@SuppressWarnings("unused")
	private static Matrix do_Down_Sample(Matrix coords, int dsf)
	{
		int rows = coords.getRowDimension();
		int cols = coords.getColumnDimension();
		int size = (cols / downsampleFactor);
		int[] column_indices = new int[size];
		for (int i = 0; i < size; i++)
			{
				column_indices[i] = (i * downsampleFactor);
			}
		return coords.getMatrix(0, rows - 1, column_indices);
	}

	static double[][] appendArray2D(double[][] array1, double[][] array2)
	{
		int a = array1[0].length, b = array2[0].length;

		double[][] result = new double[Math.max(array1.length, array2.length)][a + b];

		int i;
		for (i = 0; i < array1.length && i < array2.length; i++)
			{
				if (array1[i].length != a || array2[i].length != b)
					{
						throw new IllegalArgumentException("The number of ROWS does not match at index: " + i);
					}
				System.arraycopy(array1[i], 0, result[i], 0, a);
				System.arraycopy(array2[i], 0, result[i], a, b);
			}
		return result;
	}

	public static void main(String[] args)
	{
		rm = RoundingMode.HALF_UP;
		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);

		long startTime, endTime, totalTime, poolStart, poolEnd;

		poolStart = System.nanoTime();

		System.out.println("Running the Pool Driver: ");
		System.out.println("Getting the input file: ");
		try
			{
				input_path = "POOL.txt";
				String WD = System.getProperty("user.dir");
				String in_path = WD + File.separator + input_path;

				if (args.length == 1)
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
		System.out.println("Assigning parameters... ");
		assign_parameters();

		Col_sum = 0;
		for (int i = 0; i < number_of_input_matrices; i++)
			{
				Col_sum += (end_Indices.get(i) - start_Indices.get(i) + 1);
			}

		if (doDownSample) Col_sum = (Col_sum / downsampleFactor);
		System.out.println("Col_sum = " + Col_sum);
		System.out.println("-------------------------------------------------------------------------------------");

		int size;
		double[][] appended = null;
		for (int i = 0; i < number_of_input_matrices; i++) // For every matrix to combine:
			{
				// Read the matrix from file:
				System.out.println("Reading Matrix number: " + (i + 1));

				startTime = System.nanoTime();

				Matrix input_matrix = Matrix_IO.read_Matrix(paths.get(i));
				ROWS = input_matrix.getRowDimension();
				COLS = input_matrix.getColumnDimension();

				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done: " + (totalTime / 1.000E9) + " seconds");

				// Get the specified frames:
				int colStart = (start_Indices.get(i) - 1);
				int colEnd = (end_Indices.get(i) - 1);
				// System.out.println("colStart = " + colStart);
				// System.out.println("colEnd = " + colEnd); (totalTime / 1.000E9) + " seconds)

				Matrix X = null;

				if (doDownSample)
					{
						System.out.println("Down-Sampling by factor : " + downsampleFactor);

						startTime = System.nanoTime();

						int[] columns = do_DownSample(colStart, colEnd, downsampleFactor);

						// for (int z : columns) System.out.println("column " + z);

						X = input_matrix.getMatrix(0, ROWS - 1, columns);
						size = X.getColumnDimension();

						endTime = System.nanoTime();
						totalTime = endTime - startTime;

						System.out.println("\tThe number of frames in the original matrix is: " + COLS);
						System.out.println("\tThe number of frames in the reduced matrix will be: " + size);
						System.out.println("Done: " + (totalTime / 1.000E9) + " seconds");
					}
				else
					{
						System.out.println("Getting matrix:");

						startTime = System.nanoTime();

						X = input_matrix.getMatrix(0, ROWS - 1, colStart, colEnd);
						size = X.getColumnDimension();

						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");
					}

				System.out.println("Getting matrix as 2D array : ");
				startTime = System.nanoTime();

				double[][] xArray = X.getArray();

				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");

				if (i == 0) first_array = xArray;

				if (i > 0)
					{
						System.out.println("Appending array : ");

						startTime = System.nanoTime();

						appended = appendArray2D(first_array, xArray);
						first_array = appended;

						endTime = System.nanoTime();
						totalTime = endTime - startTime;
						System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");
					}
				System.out.println("-------------------------------------------------------------------------------------");
				systemInfo sys = new systemInfo();
				System.out.println(sys.Info());
				System.out.println("-------------------------------------------------------------------------------------");
			}

		System.out.println("Constructing augmented matrix from the 2D array : ");

		startTime = System.nanoTime();

		augmented_matrix = new Matrix(appended, ROWS, Col_sum);

		endTime = System.nanoTime();
		totalTime = endTime - startTime;
		System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");

		if (doDownSample)
			{
				if (compress) path = out_dir + "Pooled_Coordinates_Matrix_" + description + "_" + number_of_input_matrices + "_DS_" + downsampleFactor + ".txt.bz2";
				else
					path = out_dir + "Pooled_Coordinates_Matrix_" + description + "_" + number_of_input_matrices + "_DS_" + downsampleFactor + ".txt";
			}
		else if (compress) path = out_dir + "Pooled_Coordinates_Matrix_" + description + "_" + number_of_input_matrices + ".txt.bz2";
		else
			path = out_dir + "Pooled_Coordinates_Matrix_" + description + "_" + number_of_input_matrices + ".txt";

		if (compress)
			{
				System.out.println("Writing Compressed BZ2 Output: ");
				startTime = System.nanoTime();

				Matrix_IO.write_BZ2_Matrix(augmented_matrix, path, 9, 3);

				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");
				systemInfo sys = new systemInfo();
				System.out.println(sys.Info());
			}
		else
			{
				System.out.println("Writing un-compressed output: ");
				startTime = System.nanoTime();

				Matrix_IO.write_Matrix(augmented_matrix, path, 9, 3);

				endTime = System.nanoTime();
				totalTime = endTime - startTime;
				System.out.println("Done: " + nf3.format(totalTime / (1.000E9)) + " seconds");
				systemInfo sys = new systemInfo();
				System.out.println(sys.Info());
			}

		poolEnd = System.nanoTime();
		totalTime = poolEnd - poolStart;

		// augmented_matrix = null;
		// System.gc();

		initialize_Log();
		write_Log();

		System.out.println("Done. ");
		System.out.println("Total Run Time = " + nf3.format(totalTime / (60.000E9)) + " minutes.");
		System.out.println("-------------------------------------------------------------------------------------------------------------------------------------");
	}
}
