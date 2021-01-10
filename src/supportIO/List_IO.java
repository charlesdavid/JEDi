package supportIO;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.compress.compressors.CompressorException;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;

import Jama.Matrix;
import support.Residue_ID_Pair;

/**
 * JED class List_IO: Handles reading and writing of lists. Copyright (C) 2012 Dr. Charles David
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

public class List_IO
{
	/**
	 * Method to read a list from a file.
	 * 
	 * @param path The full path to the list file
	 * @param type The data type (String, Integer, etc.)
	 * @return The list of data Type 'type'
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static List read_List_From_File(String path, String type)
	{
		List data_in = null;
		File list_in = new File(path);
		BufferedReader reader = null;
		try
			{
				reader = new BufferedReader(new FileReader(list_in));
			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the file: " + path + list_in);
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
						System.err.println("Could not read the file: " + path + list_in);
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
						System.err.println("Expected list of integers: " + path + list_in);
						e.printStackTrace();
					}
				catch (IOException e)
					{
						System.err.println("Could not read the file: " + path + list_in);
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
						System.err.println("Expected list of decimals: " + path + list_in);
						e.printStackTrace();
					}
				catch (IOException e)
					{
						System.err.println("Could not read the file: " + path + list_in);
						e.printStackTrace();
					}
			}
		else
			try
				{
					data_in = new ArrayList<Object>();
					while ((line = reader.readLine()) != null)
						{
							data_in.add(line);
						}
				}
			catch (IOException e)
				{
					System.err.println("Could not read the list of objects: " + path + list_in);
					e.printStackTrace();
				}

		return data_in;
	}

	/**
	 * Method to read a file as a list of Strings.
	 * 
	 * @param path The full path to the list file
	 * @return The list of lines
	 */
	public static List<String> read_Lines_From_File(String path)
	{
		List<String> data_in = null;
		File list_in = new File(path);
		BufferedReader reader = null;
		try
			{
				reader = new BufferedReader(new FileReader(list_in));
				String line;
				data_in = new ArrayList<String>();
				while ((line = reader.readLine()) != null)
					{
						data_in.add(line);
					}
			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the file: " + path + list_in);
				e.printStackTrace();
			}
		catch (IOException e)
			{
				System.err.println("Could not read the file: " + path + list_in);
				e.printStackTrace();
			}
		return data_in;
	}

	/**
	 * Method to read a Compressed (GZ or BZ2) file as a list of Strings.
	 * 
	 * The files must be standard .bz2 or .gz files.
	 * 
	 * @param path The full path to the file
	 * @return The list of file lines
	 */
	public static List<String> read_Lines_From_BZ_GZ_File(String path)
	{
		List<String> data_in = new ArrayList<String>();
		String line;
		File list_in = new File(path);
		BufferedReader reader = null;
		try
			{
				FileInputStream fin = new FileInputStream(path);
				BufferedInputStream bis = new BufferedInputStream(fin);
				CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
				reader = new BufferedReader(new InputStreamReader(input));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						data_in.add(line);
						// System.out.println(line);
					}
				reader.close();
			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the file: " + path + list_in);
				e.printStackTrace();
			}
		catch (IOException e)
			{
				System.err.println("Could not read the file: " + path + list_in);
				e.printStackTrace();
			}
		catch (CompressorException e)
			{
				System.err.println("Could not decode the file... No valid decoder: " + path + list_in);
				e.printStackTrace();
			}
		return data_in;
	}

	/**
	 * Method to write a list of String-Integer to file (For example, Residue ID pairs)
	 * 
	 * @param String_Data
	 * @param Integer_Data
	 * @param path         The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_String_Integer_List(List String_Data, List Integer_Data, String path)
	{
		try
			{
				File list_out = new File(path);
				BufferedWriter writer = new BufferedWriter(new FileWriter(list_out));
				int index = 0;
				for (Object o : String_Data)
					{
						String s = (String) o;
						Object oo = Integer_Data.get(index);
						int i = (Integer) oo;
						writer.write(s + "\t" + i + "\n");
						index++;
					}
				writer.flush();
				writer.close();

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
	 * @param data The list of type Integer
	 * @param path The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_Integer_List(List data, String path)
	{
		try
			{
				File list_out = new File(path);
				BufferedWriter list_writer = new BufferedWriter(new FileWriter(list_out));
				for (Object o : data)
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
	 * @param data       The list a numbers
	 * @param path       The full path to write the file
	 * @param dec_places The number of decimal places to retain
	 */
	@SuppressWarnings("rawtypes")
	public static void write_Double_List(List data, String path, int dec_places)
	{
		try
			{
				File list_out = new File(path);
				BufferedWriter list_writer = new BufferedWriter(new FileWriter(list_out));
				NumberFormat nf = NumberFormat.getInstance();
				nf.setGroupingUsed(false); // This should eliminate the ',' thousands separators
				RoundingMode rm = RoundingMode.HALF_UP;
				nf.setRoundingMode(rm);
				nf.setMaximumFractionDigits(dec_places);
				nf.setMinimumFractionDigits(dec_places);
				for (Object o : data)
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
	 * @param data The list of type string
	 * @param path The full path to write the file
	 */
	@SuppressWarnings("rawtypes")
	public static void write_String_List(List data, String path)
	{
		try
			{
				File list_out = new File(path);
				BufferedWriter list_writer = new BufferedWriter(new FileWriter(list_out));
				for (Object o : data)
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
	 * Method to write a list of ResidueID Pairs to file.
	 * 
	 * @param data The list of type ResidueID_Pair
	 * @param path The full path to write the file
	 */
	public static void write_ResidueID_Pair_List(List<Residue_ID_Pair> data, String path)
	{
		try
			{
				File list_out = new File(path);
				BufferedWriter list_writer = new BufferedWriter(new FileWriter(list_out));
				for (Object o : data)
					{
						Residue_ID_Pair pair = (Residue_ID_Pair) o;
						list_writer.write(pair.toShortString() + "\n");
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
	 * @param data The list of type Matrix
	 * @param path The full path to write the file
	 */
	public static void write_Matrix_List(List<Matrix> data, String path, int width, int dec)
	{
		try
			{
				File list_out = new File(path);
				PrintWriter matrix_writer = new PrintWriter(list_out);
				for (Matrix m : data)
					{
						m.print(matrix_writer, width, dec);
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

	@SuppressWarnings("rawtypes")
	/**
	 * Method to print the elements of a list to the output stream.
	 * 
	 * @param list The list of elements
	 */
	public static void print_elements_in_List(List list)
	{
		for (Object o : list)
			{
				System.out.println(o);
			}
	}

	@SuppressWarnings("rawtypes")
	/**
	 * Method to print the elements of a list to the output stream. The elements are printed via their 'toString' method.
	 * 
	 * @param list The list of elements
	 */
	public static void print_elements_in_List_to_String(List list)
	{
		for (Object o : list)
			{
				System.out.println(o.toString());
			}
	}
}
