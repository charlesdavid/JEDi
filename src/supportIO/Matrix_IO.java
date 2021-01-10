package supportIO;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.compress.compressors.CompressorException;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream;

import Jama.Matrix;

/**
 * JED class Matrix_IO: Handles reading and writing of matrices. Copyright (C) 2012 Dr. Charles David
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
 */

public class Matrix_IO
{
	/**
	 * Method to read in a matrix from a specified path:
	 * 
	 * Supports ZIP, BZIP2, and GZIP compression algorithms.
	 * 
	 * The file suffix must be ".txt"
	 * 
	 * @param path The path to the matrix
	 * @return The Jama matrix
	 */
	public static Matrix read_Matrix(String path)
	{
		Matrix data_in = null;
		try
			{
				if (path.endsWith("bz2") || path.endsWith(".gz"))
					{
						FileInputStream fin = new FileInputStream(path);
						BufferedInputStream bis = new BufferedInputStream(fin, 1024 * 1000000);
						CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
						BufferedReader reader = new BufferedReader(new InputStreamReader(input), 1024 * 1000000);
						data_in = Matrix.read(reader);
						reader.close();
						fin.close();
						bis.close();
						input.close();
					}
				if (path.endsWith(".zip"))
					{
						ZipFile zipFile = new ZipFile(path);
						Enumeration<? extends ZipEntry> entries = zipFile.entries();
						ZipEntry entry = entries.nextElement();
						InputStream zip = zipFile.getInputStream(entry);
						// System.out.println("\tProcessing file: " + entry.getName() + "\n");
						BufferedReader reader = new BufferedReader(new InputStreamReader(zip, "UTF-8"));
						reader.close();
						zip.close();
						zipFile.close();
					}
				if (path.endsWith(".txt"))
					{
						File matrix_in = new File(path);
						BufferedReader reader = new BufferedReader(new FileReader(matrix_in));
						data_in = Matrix.read(reader);
						reader.close();
					}
			}

		catch (IOException e)
			{
				System.err.println("The file: " + path + " could not be found.");
				e.printStackTrace();
			}
		catch (CompressorException e)
			{
				System.err.println("There was a compression error.");
				e.printStackTrace();
			}
		return data_in;
	}

	/**
	 * Method to read in a matrix from a specified directory and filename:
	 * 
	 * @param dir      The directory of the matrix
	 * 
	 * @param filename The name of the matrix to read
	 * @return The Jama matrix
	 */
	public static Matrix read_Matrix(String directory, String filename)
	{
		File in_file = new File(directory + filename);
		Matrix data_in = null;
		try
			{
				BufferedReader reader = new BufferedReader(new FileReader(in_file));
				data_in = Matrix.read(reader);
			}
		catch (IOException e)
			{
				System.err.println("Could not read the file: " + directory + filename);
				e.printStackTrace();
			}
		return data_in;
	}

	/**
	 * Method to write a matrix to a specified directory and filename:
	 * 
	 * @param data     The matrix to write to file
	 * 
	 * @param dir      The directory
	 * 
	 * @param filename The name of the matrix
	 */
	public static void write_Matrix_adaptive_spacing(Matrix data, String path)
	{
		int w, d;
		double test = data.get(0, 0);
		if (test <= 1)
			{
				w = 20;
				d = 12;
			}
		else if (test <= 10 && test > 1)
			{
				w = 20;
				d = 10;
			}
		else if (test <= 100 && test > 10)
			{
				w = 20;
				d = 8;
			}
		else if (test <= 1000 && test > 100)
			{
				w = 20;
				d = 6;
			}
		else
			{
				w = 20;
				d = 3;
			}

		try
			{
				if (path.endsWith(".txt"))
					{
						File matrix_out = new File(path);
						PrintWriter matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out), 1024 * 10000));
						data.print(matrix_writer, w, d);
						matrix_writer.flush();
						matrix_writer.close();
					}
				if (path.endsWith(".bz2"))
					{
						FileOutputStream output = new FileOutputStream(path);
						PrintWriter writer = new PrintWriter(
								new BufferedWriter(new OutputStreamWriter(new BZip2CompressorOutputStream(new BufferedOutputStream(output, 1024 * 900)), "UTF-8"), 1024 * 90000));
						data.print(writer, w, d);
						writer.flush();
						output.flush();
						writer.close();
						output.close();
					}
				if (path.endsWith(".gz"))
					{
						FileOutputStream output = new FileOutputStream(path);
						PrintWriter writer = new PrintWriter(
								new BufferedWriter(new OutputStreamWriter(new GzipCompressorOutputStream(new BufferedOutputStream(output, 1024 * 900)), "UTF-8"), 1024 * 90000));
						data.print(writer, w, d);
						writer.flush();
						output.flush();
						writer.close();
						output.close();
					}
			}
		catch (IOException e)
			{
				System.err.println("Could not write the file: " + path);
				e.printStackTrace();
			}
	}

	/**
	 * Method to write a matrix to a specified directory and filename using a specified field width and number of decimal places:
	 * 
	 * @param data       The matrix to write to file
	 * 
	 * @param directory  The directory
	 * 
	 * @param filename   The name of the matrix
	 * 
	 * @param width      The width of the matrix columns
	 * 
	 * @param dec_places The number of decimal places to use
	 */
	public static void write_Matrix(Matrix data, String directory, String filename, int width, int dec_places)
	{
		try
			{
				File matrix_out = new File(directory + filename);
				PrintWriter matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out), 1024 * 1000));
				data.print(matrix_writer, width, dec_places);
				matrix_writer.flush();
				matrix_writer.close();
			}
		catch (IOException e)
			{
				System.err.println("Could not write the file: " + directory + filename);
				e.printStackTrace();
			}
	}

	/**
	 * Method to write a matrix to a specified path using a specified field width and number of decimal places:
	 * 
	 * @param data       The matrix to write to file
	 * 
	 * @param path       The path to the matrix file
	 * 
	 * @param width      The width of the matrix columns
	 * 
	 * @param dec_places The number of decimal places to use
	 */
	public static void write_Matrix(Matrix data, String path, int width, int dec_places)
	{
		try
			{
				if (path.endsWith(".txt"))
					{
						File matrix_out = new File(path);
						PrintWriter matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out), 1024 * 10000));
						data.print(matrix_writer, width, dec_places);
						matrix_writer.flush();
						matrix_writer.close();
					}
				if (path.endsWith(".bz2"))
					{
						FileOutputStream output = new FileOutputStream(path);
						PrintWriter writer = new PrintWriter(
								new BufferedWriter(new OutputStreamWriter(new BZip2CompressorOutputStream(new BufferedOutputStream(output, 1024 * 900)), "UTF-8"), 1024 * 90000));
						data.print(writer, width, dec_places);
						writer.flush();
						output.flush();
						writer.close();
						output.close();
					}
				if (path.endsWith(".gz"))
					{
						FileOutputStream output = new FileOutputStream(path);
						PrintWriter writer = new PrintWriter(
								new BufferedWriter(new OutputStreamWriter(new GzipCompressorOutputStream(new BufferedOutputStream(output, 1024 * 900)), "UTF-8"), 1024 * 90000));
						data.print(writer, width, dec_places);
						writer.flush();
						output.flush();
						writer.close();
						output.close();
					}
			}
		catch (IOException e)
			{
				System.err.println("Could not write the file: " + path);
				e.printStackTrace();
			}
	}

	/**
	 * Method to write a BZIP2 Compressed Matrix to a specified path using a specified field width and number of decimal places:
	 * 
	 * @param data       The matrix to write to file
	 * 
	 * @param path       The path to the matrix file
	 * 
	 * @param width      The width of the matrix columns
	 * 
	 * @param dec_places The number of decimal places to use
	 */
	public static void write_BZ2_Matrix(Matrix data, String path, int width, int dec_places)
	{
		try
			{
				FileOutputStream output = new FileOutputStream(path);
				PrintWriter writer = new PrintWriter(
						new BufferedWriter(new OutputStreamWriter(new BZip2CompressorOutputStream(new BufferedOutputStream(output, 1024 * 900)), "UTF-8"), 1024 * 90000));
				data.print(writer, width, dec_places);
				writer.flush();
				output.flush();
				writer.close();
				output.close();
			}
		catch (IOException e)
			{
				System.err.println("Could not write the file: " + path);
				e.printStackTrace();
			}
	}
}
