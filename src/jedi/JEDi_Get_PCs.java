package jedi;

import java.util.List;

import Jama.Matrix;
import support.Projector;

/**
 * JED class JED_Get_PCs: Constructs the principle components.
 * 
 * Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_PCs
{
	final double FLOOR = 1.000E-16;
	int number_of_atoms;
	final int number_of_modes, number_of_conformations, ROWS;
	final List<Double> eigenvalues;
	final Matrix centered_coords, top_evectors;
	Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;

	/**
	 * Constructor to create the PCs.
	 *
	 */
	JEDi_Get_PCs(Matrix centered_data, Matrix evects, List<Double> evals)
	{
		this.centered_coords = centered_data;
		this.top_evectors = evects;
		this.eigenvalues = evals;

		number_of_modes = top_evectors.getColumnDimension();
		number_of_conformations = centered_coords.getColumnDimension();
		ROWS = top_evectors.getRowDimension();
		number_of_atoms = (ROWS / number_of_modes);
	}

	/**
	 * Computes the Principle Components and writes them to files: Un-normed, normed, weighted, weighted-normed
	 */
	void get_PCs()
	{
		projections = new Matrix(number_of_conformations, number_of_modes);
		normed_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_normed_projections = new Matrix(number_of_conformations, number_of_modes);

		for (int mode = 0; mode < number_of_modes; mode++)
			{
				int row_index_1 = 0;
				int row_index_2 = ROWS - 1;
				Matrix data1 = top_evectors.getMatrix(row_index_1, row_index_2, mode, mode);
				Matrix vector1 = Projector.get_Normed_arrayF(data1);
				double weight = eigenvalues.get(mode);
				if (weight < FLOOR) weight = FLOOR;
				double sqrt_weight = Math.sqrt(weight);
				for (int conf = 0; conf < number_of_conformations; conf++)
					{
						Matrix data2 = centered_coords.getMatrix(row_index_1, row_index_2, conf, conf);
						Matrix vector2 = Projector.get_Normed_arrayF(data2);
						double dp = Projector.get_InnerProduct(data1, data2); // Units are Angstroms
						double normed_dp = Projector.get_InnerProduct(vector1, vector2); // Unitless
						double w_dp = (sqrt_weight * dp);  // Units are Angstroms Squared
						double weighted_normed_dp = (sqrt_weight * normed_dp);  // Units are Angstroms
						projections.set(conf, mode, dp);
						normed_projections.set(conf, mode, normed_dp);
						weighted_projections.set(conf, mode, w_dp);
						weighted_normed_projections.set(conf, mode, weighted_normed_dp);
					}
			}
	}

	/* ******************************************* GETTERS ***************************************************** */

	public Matrix getPCs()
	{
		return projections;
	}

	public Matrix getNormed_PCs()
	{
		return normed_projections;
	}

	public Matrix getWeighted_normed_PCs()
	{
		return weighted_normed_projections;
	}

	public Matrix getWeighted_PCs()
	{
		return weighted_projections;
	}
}
