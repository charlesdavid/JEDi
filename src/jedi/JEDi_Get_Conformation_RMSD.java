package jedi;

import Jama.Matrix;

/**
 * JED class JED_Get_RMSD: Calculates the RMSD for 2 input conformations.
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Conformation_RMSD
{
	final int ROWS, COLS, number_of_atoms;
	double RMSD;
	final Matrix ref_Structure, field_Structure;

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
		this.ref_Structure = ref_conf;
		this.field_Structure = conf;
		this.ROWS = ref_Structure.getRowDimension();
		this.COLS = ref_Structure.getColumnDimension();
		this.number_of_atoms = (ROWS / 3);
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
