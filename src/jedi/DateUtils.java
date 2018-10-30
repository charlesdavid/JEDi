package jedi;

import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * JED class DateUtils: Returns the current time as a date object. Copyright (C) 2012 Dr. Charles David
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
public class DateUtils
	{

		public static final String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";
		public static String current_date;

		/**
		 * Returns the current time as a date.
		 * 
		 * @return current_date
		 */
		public static String now()
			{
				Calendar cal = Calendar.getInstance();
				SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
				current_date = sdf.format(cal.getTime());
				return current_date;
			}

		public static void main(String arg[])
			{

				System.out.println("Now : " + DateUtils.now());
			}
	}
