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
import java.util.List;

import Jama.Matrix;

/**
 * JED class List_IO: Handles reading and writing of lists.
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

public class List_IO
{

	static String directory, name, full_path, type, line;
	@SuppressWarnings("rawtypes")
	static List data_in, data_out;
	static File list_in, list_out;
	static BufferedReader reader;
	static BufferedWriter list_writer;
	static NumberFormat nf;
	static RoundingMode rm;


	/**
	 * Method to read a list from a file.
	 * 
	 * @param path
	 *            The full path to the list file
	 * @param T
	 *            The data type (String, Integer, etc.)
	 * @return The list
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static List read_List(String path, String T)
		{

			full_path = path;
			type = T;
			list_in = new File(full_path);
			try
				{
					reader = new BufferedReader(new FileReader(list_in));
				}
			catch (FileNotFoundException e)
				{
					System.err.println("Could not find the file: " + full_path + list_in);
					e.printStackTrace();
				}
			String line;
			if (type.equals("String"))
				{
					try
						{
							data_in = new ArrayList<String>();
							while ((line = reader.readLine()) != null)
								{
									data_in.add(line);
								}
						}
					catch (IOException e)
						{
							System.err.println("Could not read the file: " + full_path + list_in);
							e.printStackTrace();
						}
				}
			else if (type.equals("Integer"))
				{
					try
						{
							data_in = new ArrayList<Integer>();
							while ((line = reader.readLine()) != null)
								{
									data_in.add(Integer.parseInt(line));
								}
						}
					catch (NumberFormatException e)
						{
							System.err.println("Expected list of integers: " + full_path + list_in);
							e.printStackTrace();
						}
					catch (IOException e)
						{
							System.err.println("Could not read the file: " + full_path + list_in);
							e.printStackTrace();
						}
				}
			else if (type.equals("Double"))
				{
					try
						{
							data_in = new ArrayList<Double>();
							while ((line = reader.readLine()) != null)
								{
									data_in.add(Double.parseDouble(line));
								}
						}
					catch (NumberFormatException e)
						{
							System.err.println("Expected list of decimals: " + full_path + list_in);
							e.printStackTrace();
						}
					catch (IOException e)
						{
							System.err.println("Could not read the file: " + full_path + list_in);
							e.printStackTrace();
						}
				}
			else try
				{
					data_in = new ArrayList<Object>();
					while ((line = reader.readLine()) != null)
						{
							data_in.add(line);
						}
				}
			catch (IOException e)
				{
					System.err.println("Could not read the list of objects: " + full_path + list_in);
					e.printStackTrace();
				}

			return data_in;

		}

	/**
	 * Method to write a list of String-Integer to file (For example, Residue ID pairs)
	 * 
	 * @param String_Data
	 * @param Integer_Data
	 * @param path
	 *            The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_String_Integer_List(List String_Data, List Integer_Data, String path)
		{
			try
				{
					list_out = new File(path);
					list_writer = new BufferedWriter(new FileWriter(list_out));
					int index = 0;
					for (Object o : String_Data)
						{
							String s = (String) o;
							Object oo = Integer_Data.get(index);
							int i = (Integer) oo;
							list_writer.write(s + "\t" + i + "\n");
							index++;
						}

					list_writer.flush();
					list_writer.close();

				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + path);
					e.printStackTrace();
				}
		}

	/**
	 * Method to write a list of Integers to file.
	 * 
	 * @param data
	 *            The list of type Integer
	 * @param path
	 *            The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_Integer_List(List data, String path)
		{

			data_out = data;
			list_out = new File(path);
			try
				{
					list_writer = new BufferedWriter(new FileWriter(list_out));
					for (Object o : data_out)
						{
							Integer i = (Integer) o;
							list_writer.write(i + "\n");
						}
					list_writer.flush();
					list_writer.close();
				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + path);
					e.printStackTrace();
				}
		}

	/**
	 * Method to write a list of decimal numbers to file
	 * 
	 * @param data
	 *            The list a numbers
	 * @param path
	 *            The full path to write the file
	 * @param dec_places
	 *            The number of decimal places to retain
	 */
	@SuppressWarnings("rawtypes")
	public static void write_Double_List(List data, String path, int dec_places)
		{

			data_out = data;
			list_out = new File(path);
			try
				{
					list_writer = new BufferedWriter(new FileWriter(list_out));
					nf = NumberFormat.getInstance();
					rm = RoundingMode.HALF_UP;
					nf.setRoundingMode(rm);
					nf.setMaximumFractionDigits(dec_places);
					nf.setMinimumFractionDigits(dec_places);
					for (Object o : data_out)
						{
							Double d = (Double) o;
							list_writer.write(nf.format(d) + "\n");
						}
					list_writer.flush();
					list_writer.close();

				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + path);
					e.printStackTrace();
				}
		}

	/**
	 * Method to write a list of strings to file.
	 * 
	 * @param data
	 *            The list of type string
	 * @param path
	 *            The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_String_List(List data, String path)
		{

			data_out = data;
			list_out = new File(path);
			try
				{
					list_writer = new BufferedWriter(new FileWriter(list_out));
					for (Object o : data_out)
						{
							String s = (String) o;
							list_writer.write(s + "\n");
						}
					list_writer.flush();
					list_writer.close();

				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + path);
					e.printStackTrace();
				}
		}

	/**
	 * Method to write a list of strings to file.
	 * 
	 * @param data
	 *            The list of type Matrix
	 * @param path
	 *            The full path to write the file
	 */
	public static void write_Matrix_List(List<Matrix> data, String path)
		{

			list_out = new File(path);
			try
				{
					PrintWriter matrix_writer = new PrintWriter(list_out);
					for (Matrix m : data)
						{
							m.print(matrix_writer, 9, 3);
						}
					matrix_writer.flush();
					matrix_writer.close();
				}
			catch (IOException e)
				{
					System.err.println("Could not write the file: " + path);
					e.printStackTrace();
				}
		}
}
