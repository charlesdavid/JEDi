package jedi;

import Jama.Matrix;

/**
 * JED class JED_Get_Coordinates_from_Matrix: Gets the trajectory(s) coordinates from a matrix. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Coordinates_from_Matrix
	{

		String directory, input_file_name, description;
		int number_of_residues, number_of_conformations, ROWS, COLS;
		Matrix original_PDB_coordinates;

		/* ************************************** CONSTRUCTORS ******************************************************************************** */

		/**
		 * Constructor for reading in the coordinates matrix
		 * 
		 * @param dir  The working directory
		 * @param name The name of the Coordinates Matrix File
		 */
		public JEDi_Get_Coordinates_from_Matrix(String dir, String name)
			{
				directory = dir;
				input_file_name = name;
			}

		/* ************************************** METHODS ******************************************************************************** */

		/**
		 * Reads and returns the matrix of PDB coordinates
		 * 
		 * @return The matrix of PDB coordinates
		 */
		public Matrix get_Original_PDB_coordinates()
			{
				original_PDB_coordinates = Matrix_IO.read_Matrix(directory, input_file_name);
				ROWS = original_PDB_coordinates.getRowDimension();
				COLS = original_PDB_coordinates.getColumnDimension();
				number_of_conformations = COLS;
				number_of_residues = (ROWS / 3);
				return original_PDB_coordinates;
			}
	}
