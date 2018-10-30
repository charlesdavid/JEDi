package jedi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Jama.Matrix;

/**
 * JED class Matrix_IO: Handles reading and writing of matrices.
 * Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 */

public class Matrix_IO
{

	static String directory, name;
	static Matrix data_in, data_out;
	static File matrix_in, matrix_out;
	static BufferedReader reader;
	static PrintWriter matrix_writer;

	/**
	 * Method to read in a matrix from a specified path:
	 * 
	 * @param path
	 *            The path to the matrix
	 * @return The Jama matrix
	 */
	public static Matrix read_Matrix(String path)
		{

			matrix_in = new File(path);
			try
				{
					reader = new BufferedReader(new FileReader(matrix_in));
				}
			catch (FileNotFoundException e)
				{
					System.err.println("The file: " + directory + matrix_in + " could not be found.");
					e.printStackTrace();
				}
			try
				{
					data_in = Matrix.read(reader);
				}
			catch (IOException e)
				{
					System.err.println("Could not read the file: " + directory + matrix_in);
					e.printStackTrace();
				}
			return data_in;
		}

	/**
	 * Method to read in a matrix from a specified directory and filename:
	 * 
	 * @param dir
	 *            The directory of the matrix
	 * 
	 * @param filename
	 *            The name of the matrix to read
	 * @return The Jama matrix
	 */
	public static Matrix read_Matrix(String dir, String filename)
		{

			directory = dir;
			name = filename;
			matrix_in = new File(dir + name);
			try
				{
					reader = new BufferedReader(new FileReader(matrix_in));
				}
			catch (FileNotFoundException e)
				{
					System.err.println("The file: " + directory + matrix_in + " could not be found.");
					e.printStackTrace();
				}
			try
				{
					data_in = Matrix.read(reader);
				}
			catch (IOException e)
				{
					System.err.println("Could not read the file: " + directory + matrix_in);
					e.printStackTrace();
				}
			return data_in;
		}

	/**
	 * Method to write a matrix to a specified directory and filename:
	 * 
	 * @param data
	 *            The matrix to write to file
	 * 
	 * @param dir
	 *            The directory
	 * 
	 * @param filename
	 *            The name of the matrix
	 */
	public static void write_Matrix(Matrix data, String dir, String filename)
		{

			directory = dir;
			name = filename;
			data_out = data;
			matrix_out = new File(directory + name);
			try
				{
					matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out)));
				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + directory + matrix_out);
					e.printStackTrace();
				}
			data_out.print(matrix_writer, 12, 6);
			matrix_writer.flush();
			matrix_writer.close();
		}

	/**
	 * Method to write a matrix to a specified directory and filename using a specified field width and number of decimal places:
	 * 
	 * @param data
	 *            The matrix to write to file
	 * 
	 * @param dir
	 *            The directory
	 * 
	 * @param filename
	 *            The name of the matrix
	 * 
	 * @param width
	 *            The width of the matrix columns
	 * 
	 * @param dec_places
	 *            The number of decimal places to use
	 */
	public static void write_Matrix(Matrix data, String dir, String filename, int width, int dec_places)
		{

			directory = dir;
			name = filename;
			data_out = data;
			matrix_out = new File(directory + name);
			try
				{
					matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out)));
				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + directory + matrix_out);
					e.printStackTrace();
				}
			data_out.print(matrix_writer, width, dec_places);
			matrix_writer.flush();
			matrix_writer.close();
		}

	/**
	 * Method to write a matrix to a specified path using a specified field width and number of decimal places:
	 * 
	 * @param data
	 *            The matrix to write to file
	 * 
	 * @param path
	 *            The path to the matrix file
	 * 
	 * @param width
	 *            The width of the matrix columns
	 * 
	 * @param dec_places
	 *            The number of decimal places to use
	 */
	public static void write_Matrix(Matrix data, String path, int width, int dec_places)
		{

			data_out = data;
			matrix_out = new File(path);
			try
				{
					matrix_writer = new PrintWriter(new BufferedWriter(new FileWriter(matrix_out)));
				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + directory + matrix_out);
					e.printStackTrace();
				}
			data_out.print(matrix_writer, width, dec_places);
			matrix_writer.flush();
			matrix_writer.close();
		}
}
