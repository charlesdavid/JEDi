package jedi;

import java.util.List;

import Jama.Matrix;
import support.Projector;

/**
 * JED class JED_Get_Cartesian_PCs: Constructs the Principle Components for the Cartesian PCA on the residues. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Residue_Pair_PCs
{
	boolean exist, success;
	final double FLOOR = 1.000E-16;
	int Residue_Index;
	final int number_of_modes, number_of_conformations, ROWS, number_of_atoms;
	final List<Double> eigenvalues;
	Matrix projections;
	final Matrix centered_coords, top_evectors;

	/* ************************************** CONSTRUCTOR ******************************************************************************** */

	/**
	 * Constructor to create the delta vector series and the delta vector projections.
	 *
	 * @param centered_data The mean centered data
	 * @param evects        The top Cartesian eigenvectors
	 * @param evals         The top Cartesian eigenvalues
	 * @param dir           The working directory
	 * @param des           The job description
	 * @param pca_model     The type of PCA: COV, CORR, PCORR (Q, R, P)
	 */
	JEDi_Get_Residue_Pair_PCs(Matrix centered_data, Matrix evects, List<Double> evals, String res_pair_index)
	{

		this.centered_coords = centered_data;
		this.top_evectors = evects;
		this.eigenvalues = evals;

		this.number_of_modes = top_evectors.getColumnDimension();
		this.number_of_conformations = centered_coords.getColumnDimension();
		this.ROWS = top_evectors.getRowDimension();
		this.number_of_atoms = (ROWS / 3);
	}

	/* ************************************** Public Methods ******************************************************************************** */

	/**
	 * Computes the Principle Components and writes them to files: Un-normed, normed, weighted, weighted-normed
	 */
	public Matrix get_PCs()
	{
		projections = new Matrix(number_of_conformations, number_of_modes);

		for (int mode = 0; mode < number_of_modes; mode++)
			{
				int row_index_1 = 0;
				int row_index_2 = ROWS - 1;
				Matrix data1 = top_evectors.getMatrix(row_index_1, row_index_2, mode, mode);
				double weight = eigenvalues.get(mode);
				if (weight < FLOOR) weight = 0;
				for (int conf = 0; conf < number_of_conformations; conf++)
					{
						Matrix data2 = centered_coords.getMatrix(row_index_1, row_index_2, conf, conf);
						double dp = Projector.get_InnerProduct(data1, data2);
						projections.set(conf, mode, dp);
					}
			}
		return projections;
	}
}
