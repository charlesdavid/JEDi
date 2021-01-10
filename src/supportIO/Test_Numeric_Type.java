package supportIO;

/**
 * JED class Test_numeric_Type: Support class for testing input data for proper format.
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
 * 
 */

public class Test_Numeric_Type
{
	/**
	 * Tests if the input is of type INTEGER.
	 * 
	 * @param input
	 *            The value to test
	 * @return TRUE if INTEGER, FALSE otherwise
	 */
	public static boolean test_Integer(String input)
		{
			try
				{
					Integer.parseInt(input);
					return true;
				} catch (Exception e)
				{
					return false;
				}
		}

	/**
	 * Tests if the input is of type DOUBLE (Java double precision decimal).
	 * 
	 * @param input
	 *            The value to test
	 * @return TRUE if DOUBLE, FALSE otherwise
	 */
	public static boolean test_Double(String input)
		{
			try
				{
					Double.parseDouble(input);
					return true;
				} catch (Exception e)
				{
					return false;
				}
		}
}
