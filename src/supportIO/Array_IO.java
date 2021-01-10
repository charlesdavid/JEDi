package supportIO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;

/**
 * JED class Array_IO: Handles reading and writing of arrays.
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

public class Array_IO
{
	/**
	 * Writes an array of doubles
	 * 
	 * @param data
	 *            the double array
	 * @param path
	 *            the path to write the file to
	 * @param dec_places
	 *            the number of decimal places to use in the numbers
	 */
	public static void write_Double_Array(double[] data, String path, int dec_places)
	{
		try
		{
			File array_out = new File(path);
			BufferedWriter array_writer = new BufferedWriter(new FileWriter(array_out));

			NumberFormat nf = NumberFormat.getInstance();
			RoundingMode rm = RoundingMode.HALF_UP;
			nf.setRoundingMode(rm);
			nf.setMaximumFractionDigits(dec_places);
			nf.setMinimumFractionDigits(dec_places);

			for (double d : data)
			{
				array_writer.write(nf.format(d) + "\n");
			}
			array_writer.flush();
			array_writer.close();
		} catch (IOException io)
		{
			System.err.println("IOException thrown. Could not write to the file: " + path);
			io.printStackTrace();
		}
	}

	/**
	 * Writes an array of integers
	 * 
	 * @param data
	 *            the integer array
	 * @param path
	 *            the path to write the file to
	 */
	public static void write_Integer_Array(int[] data, String path)
	{
		try
		{
			File array_out = new File(path);
			BufferedWriter array_writer = new BufferedWriter(new FileWriter(array_out));
			for (int i : data)
			{
				array_writer.write(i + "\n");
			}
			array_writer.flush();
			array_writer.close();
		} catch (IOException io)
		{
			System.err.println("IOException thrown. Could not write to the file: " + path);
			io.getMessage();
			io.getStackTrace();
		}
	}

	/**
	 * Writes an array of strings
	 * 
	 * @param data
	 *            the string array
	 * @param path
	 *            the path to write the file to
	 */
	public static void write_String_Array(String[] data, String path)
	{
		try
		{
			File array_out = new File(path);
			BufferedWriter array_writer = new BufferedWriter(new FileWriter(array_out));
			for (String s : data)
			{
				array_writer.write(s + "\n");
			}
			array_writer.flush();
			array_writer.close();
		} catch (IOException io)
		{
			System.err.println("IOException thrown. Could not write to the file: " + path);
			io.getMessage();
			io.getStackTrace();
		}
	}
}
