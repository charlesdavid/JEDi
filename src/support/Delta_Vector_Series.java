package support;

import Jama.Matrix;

/**
 * JEDi class Delta Vector Series computes the displacement vectors from a set of reference coordinates.
 * 
 * Copyright (C) 2012 Dr. Charles David
 */

public class Delta_Vector_Series
{
	int ROWS, COLS;
	Matrix delta_vector_series, coords, ref_coords;

	/* **************************************** CONSTRUCTOR ************************************************************************ */


	public Delta_Vector_Series(Matrix data, Matrix ref_data)
	{
		this.coords = data;
		this.ref_coords = ref_data;
		this.ROWS = coords.getRowDimension();
		this.COLS = coords.getColumnDimension();
		this.delta_vector_series = new Matrix(ROWS, COLS);
	}

	/* **************************************** METHODS ************************************************************************ */

	public Matrix get_Delta_Vector_Series()
	{
		Matrix ref_col = ref_coords;
		Matrix col = new Matrix(ROWS, 1);
		for (int i = 0; i < COLS; i++)
			{
				col = coords.getMatrix(0, ROWS - 1, i, i);
				Matrix delta = col.minus(ref_col);
				delta_vector_series.setMatrix(0, ROWS - 1, i, i, delta);
			}
		return delta_vector_series;
	}
}
