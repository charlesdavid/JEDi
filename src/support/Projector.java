package support;

import Jama.Matrix;

/**
 * Projector: Computes inner products, normed arrays, and correlations.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Projector
{
	/**
	 * Computes the NormF of the input vector and returns the normalized input vector
	 * 
	 * @param in_Vector:
	 *            The vector to be normalized (NormF)
	 * @return The NormF of the input vector
	 */
	public static Matrix get_Normed_arrayF(Matrix in_Vector)
	{
		double mag = in_Vector.normF();
		double mag_inv;
		if (mag < 1.0E-12)
		{
			mag_inv = 0;
		} else
		{
			mag_inv = Math.pow(mag, -1);
		}
		Matrix norm = new Matrix(in_Vector.getRowDimension(), in_Vector.getColumnDimension(), mag_inv);
		Matrix normed_array = in_Vector.arrayTimes(norm);
		return normed_array;
	}

	/**
	 * Computes the inner product of two equidimensional vectors.
	 * 
	 * @param vector1:
	 *            Input vector 1
	 * @param vector2:
	 *            Input vector 2
	 * @return The inner product of vector 1 and 2 (vector1 dot vector2)
	 */
	public static double get_InnerProduct(Matrix vector1, Matrix vector2)
	{
		int rows = vector1.getRowDimension();
		int rows2 = vector2.getRowDimension();
		if (!(rows == rows2)) System.err.println("Error: The row dimensions of the two vectors are not equal.");
		double inner_product = 0;
		for (int i = 0; i < rows; i++)
		{
			double value_A = vector1.get(i, 0);
			double value_B = vector2.get(i, 0);
			double product = (value_A * value_B);
			inner_product += product;
		}
		return inner_product;
	}

	/**
	 * Returns the correlation between the two vectors (lists of elements)
	 * 
	 * @param vector1:
	 *            Input vector 1
	 * @param vector2:
	 *            Input vector 2
	 * @return The Pearson correlation coefficient of vector 1 with vector 2
	 */
	public static double getCorrelation(Matrix vector1, Matrix vector2)
	{
		int index, n;
		double x_sum, y_sum, xx_sum, yy_sum, xy_sum, correlation;
		x_sum = y_sum = xx_sum = yy_sum = xy_sum = 0.0;
		index = n = vector1.getRowDimension();
		while (index-- != 0)
		{
			x_sum += vector1.get(index, 0);
			y_sum += vector2.get(index, 0);
			xx_sum += vector1.get(index, 0) * vector1.get(index, 0);
			yy_sum += vector2.get(index, 0) * vector2.get(index, 0);
			xy_sum += vector1.get(index, 0) * vector2.get(index, 0);
		}

		correlation = ((n * xy_sum - x_sum * y_sum) / (Math.sqrt((n * xx_sum - x_sum * x_sum) * (n * yy_sum - y_sum * y_sum))));

		return correlation;
	}
}
