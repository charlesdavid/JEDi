package jedi;

import java.util.List;

import Jama.Matrix;
import support.Projector;
import supportIO.Input_Parameters;

/**
 * JED class JED_Get_Cartesian_PCs: Constructs the Principle Components for the Cartesian PCA on the residues. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Residue_PCs
{
	boolean exist, success;

	final double FLOOR = Input_Parameters.FLOOR;
	final int number_of_modes, number_of_conformations, ROWS;
	final List<Double> eigenvalues;
	Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;
	final Matrix centered_coords, top_evectors;

	/* ************************************** CONSTRUCTOR ******************************************************************************** */

	/**
	 * Constructor to create Principle Components.
	 *
	 * @param centered_data The mean centered data
	 * @param evects        The top Cartesian eigenvectors
	 * @param evals         The top Cartesian eigenvalues
	 */
	JEDi_Get_Residue_PCs(Matrix centered_data, Matrix evects, List<Double> evals)
	{

		this.centered_coords = centered_data;
		this.top_evectors = evects;
		this.eigenvalues = evals;

		this.number_of_modes = top_evectors.getColumnDimension();
		this.number_of_conformations = centered_coords.getColumnDimension();
		this.ROWS = top_evectors.getRowDimension();
	}

	/* ************************************** Public Methods ******************************************************************************** */

	/**
	 * Computes the Principle Components: Un-normed, normed, weighted, weighted-normed
	 */
	public void get_PCs()
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
						double dp = Projector.get_InnerProduct(data1, data2);
						double normed_dp = Projector.get_InnerProduct(vector1, vector2);
						double w_dp = (sqrt_weight * dp);
						double weighted_normed_dp = (sqrt_weight * normed_dp);
						projections.set(conf, mode, dp);
						normed_projections.set(conf, mode, normed_dp);
						weighted_projections.set(conf, mode, w_dp);
						weighted_normed_projections.set(conf, mode, weighted_normed_dp);
					}
			}
	}

	/* ************************************** GETTERS ******************************************************************************** */

	public Matrix getProjections()
	{
		return projections;
	}

	public Matrix getNormed_projections()
	{
		return normed_projections;
	}

	public Matrix getWeighted_normed_projections()
	{
		return weighted_normed_projections;
	}

	public Matrix getWeighted_projections()
	{
		return weighted_projections;
	}
}
