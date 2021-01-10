package jedi;

import Jama.Matrix;

/**
 * JED class JED_Get_Subset: Constructs the subset of coordinates for a residue. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 * 
 */

public class JEDi_Get_Residue_Coords
{

	final int number_of_atoms, number_of_atoms_Residue, number_of_conformations, offset;
	final int ROWS, ROWS_Res, COLS, COLS_Res;
	Matrix Residue_coordinates;
	final Matrix subset_coordinates;

	// ***************************************** CONSTRUCTORS ***************************************************//

	JEDi_Get_Residue_Coords(Matrix coords, int num_atoms, int off_set)
	{
		this.subset_coordinates = coords;
		this.number_of_atoms_Residue = num_atoms;
		this.offset = off_set;

		this.ROWS = subset_coordinates.getRowDimension();
		this.COLS = subset_coordinates.getColumnDimension();
		this.number_of_atoms = (ROWS / 3);
		this.ROWS_Res = (number_of_atoms_Residue * 3);
		this.COLS_Res = COLS;
		this.number_of_conformations = COLS;
	}

	// ***************************************** METHODS ***************************************************//

	/**
	 * Calculates the matrix of PDB coordinates for a residue
	 * 
	 * @return Residue Coordinates Matrix
	 */
	Matrix get_residue_Coords()
	{
		Residue_coordinates = new Matrix(ROWS_Res, COLS_Res);

		for (int i = 0; i < number_of_atoms_Residue; i++)
		{
			for (int j = 0; j < number_of_conformations; j++)
			{
				double element_X = subset_coordinates.get((i + offset), j);
				double element_Y = subset_coordinates.get(((i + offset) + (number_of_atoms)), j);
				double element_Z = subset_coordinates.get(((i + offset) + (2 * number_of_atoms)), j);

				Residue_coordinates.set(i, j, element_X);
				Residue_coordinates.set((i + number_of_atoms_Residue), j, element_Y);
				Residue_coordinates.set((i + (2 * number_of_atoms_Residue)), j, element_Z);
			}
		}
		return Residue_coordinates;
	}
}
