package support;

import java.util.List;

import Jama.Matrix;

/**
 * JEDi class Delta_Vector_Projector: Constructs DVPs.
 * 
 * Note: DVPs are like Principle Components, except that the reference structure is an actual specified structure, not the mean structure.
 * 
 * Copyright (C) 2012 Dr. Charles David
 */

public class Delta_Vector_Projector
{
	final int number_of_modes, number_of_conformations, ROWS;
	Matrix projections, normed_projections, weighted_projections, weighted_normed_projections;
	final Matrix delta_vector_series, top_evectors;
	final List<Double> top_eigenvals;

	/* **************************************** CONSTRUCTOR ************************************************************************ */

	public Delta_Vector_Projector(Matrix dvs, Matrix evects, List<Double> top_eigenvalues)
	{
		this.delta_vector_series = dvs;
		this.top_evectors = evects;
		this.top_eigenvals = top_eigenvalues;
		this.ROWS = delta_vector_series.getRowDimension();
		this.number_of_conformations = delta_vector_series.getColumnDimension();
		this.number_of_modes = top_evectors.getColumnDimension();
	}

	/* ******************************************************************************************************************************* */

	public void get_DVPs()
	{
		projections = new Matrix(number_of_conformations, number_of_modes);
		normed_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_normed_projections = new Matrix(number_of_conformations, number_of_modes);

		for (int outer = 0; outer < number_of_modes; outer++)
			{
				for (int inner = 0; inner < number_of_conformations; inner++)
					{
						int row_index_1 = 0;
						int row_index_2 = ROWS - 1;
						Matrix data1 = top_evectors.getMatrix(row_index_1, row_index_2, outer, outer);
						Matrix data2 = delta_vector_series.getMatrix(row_index_1, row_index_2, inner, inner);
						double weight = Math.sqrt(top_eigenvals.get(outer)); // weight should have units of angstroms
						Matrix vector1 = Projector.get_Normed_arrayF(data1);
						Matrix vector2 = Projector.get_Normed_arrayF(data2);
						double dp = Projector.get_InnerProduct(data1, data2);
						double normed_dp = Projector.get_InnerProduct(vector1, vector2);
						if (dp == Double.NaN) dp = 0.000;
						if (normed_dp == Double.NaN) normed_dp = 0.000;
						double w_dp = weight * dp;
						double weighted_normed_dp = weight * normed_dp;
						projections.set(inner, outer, dp);
						normed_projections.set(inner, outer, normed_dp);
						weighted_projections.set(inner, outer, w_dp);
						weighted_normed_projections.set(inner, outer, weighted_normed_dp);
					}
			}
	}

	public Matrix get_Projections()
	{
		return projections;
	}

	public Matrix get_Normed_Projections()
	{
		return normed_projections;
	}

	public Matrix get_Weighted_Projections()
	{
		return weighted_projections;
	}

	public Matrix get_Weighted_Normed_Projections()
	{
		return weighted_normed_projections;
	}
}
