package jedi;

import Jama.Matrix;

/**
 * JED class JED_Get_RMSD: Calculates the RMSD for 2 input conformations.
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
 * along with this program. If not, see <http://www.gnu.org/license>
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Conformation_RMSD
{

	int ROWS, COLS, number_of_atoms;
	double RMSD;
	Matrix ref_Structure, field_Structure;

	/**
	 * Constructor takes 2 conformations (frames) as single column matrices
	 * Note: The matrices must have the same ROW dimension or an exception will be thrown.
	 * 
	 * @param ref_conf
	 *            The reference conformation
	 * @param conf
	 *            The other conformation
	 */
	public JEDi_Get_Conformation_RMSD(Matrix ref_conf, Matrix conf)
		{
			ref_Structure = ref_conf;
			field_Structure = conf;
			ROWS = ref_Structure.getRowDimension();
			COLS = ref_Structure.getColumnDimension();
			number_of_atoms = (ROWS / 3);
		}

	/**
	 * This method calculates the RMSD between the two conformations
	 * 
	 * @return The RMSD
	 */
	public double get_RMSD()
		{
			RMSD = 0;
			double sum_of_squares = 0;
			for (int i = 0; i < ROWS; i++)
				{
					double val1 = ref_Structure.get(i, 0);
					double val2 = field_Structure.get(i, 0);
					double diff = (val1 - val2);
					double diff_sq = diff * diff;
					sum_of_squares += diff_sq;
				}
			RMSD = Math.sqrt((sum_of_squares) / number_of_atoms);
			return RMSD;
		}
}
