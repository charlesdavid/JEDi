package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;
import support.Projector;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;

/**
 * JED class JED_Get_Distance_Pair_DVPs: Constructs the DVPs for the Residue Distance Pairs subset.
 * 
 * Note: DVPs are like Principle Components, except that the reference structure is not the mean structure.
 * 
 * Copyright (C) 2012 Dr. Charles David
 */

public class JEDi_Get_Distance_Pair_DVPs
{
	boolean exist, success;
	String directory, description, model, Q = "COV", R = "CORR", P = "PCORR";
	String out_dir, file_name_head, path;
	int number_of_modes, number_of_conformations, number_of_pairs, ROWS, COLS;
	Matrix delta_vector_series, projections, normed_projections, weighted_projections, weighted_normed_projections;
	Matrix input_coords, reference_distances, top_evectors;
	Projector proj;
	List<Double> top_eigenvals;

	/* **************************************** CONSTRUCTOR ************************************************************************ */

	/**
	 * Constructor to calculate the delta vectors and the DVPs for the Distance Pair analysis.
	 *
	 * @param data            The distances for the residue pairs
	 * @param ref_distances   The reference distances to use to create the delta vectors
	 * @param evects          The eigenvectors from the Distance Pair PCA
	 * @param top_eigenvalues The top eigenvalues form the Distance Pair PCA
	 * @param dir             The working directory
	 * @param des             The job description
	 * @param model_pca       The model type for the PCA: COV or CORR (Q, R, P)
	 * @param num_modes       The number of PCA modes to process
	 */
	public JEDi_Get_Distance_Pair_DVPs(Matrix data, Matrix ref_dists, Matrix evects, List<Double> top_eigenvalues, String model_pca)
	{
		this.input_coords = data;
		this.reference_distances = ref_dists;
		this.top_evectors = evects;
		this.top_eigenvals = top_eigenvalues;
		this.model = model_pca;

		this.directory = Input_Parameters.DIRECTORY;
		this.description = Input_Parameters.DESCRIPTION;
		this.number_of_modes = Input_Parameters.MODES_DISTANCE_PAIRS;

		this.ROWS = input_coords.getRowDimension();
		this.COLS = input_coords.getColumnDimension();
		this.number_of_conformations = COLS;
		this.delta_vector_series = new Matrix(ROWS, COLS);
		this.number_of_pairs = ROWS;

		// this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Distance_Pair_PCA" + File.separatorChar + model + File.separatorChar;
		// exist = new File(out_dir).exists();
		// if (!exist) success = (new File(out_dir)).mkdirs();
		// this.file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";
	}

	public JEDi_Get_Distance_Pair_DVPs(Matrix evects, List<Double> top_eigenvalues, String model_pca)
	{
		this.top_evectors = evects;
		this.top_eigenvals = top_eigenvalues;
		this.model = model_pca;

		this.directory = Input_Parameters.DIRECTORY;
		this.description = Input_Parameters.DESCRIPTION;
		this.number_of_modes = Input_Parameters.MODES_DISTANCE_PAIRS;

		this.ROWS = top_evectors.getRowDimension();
		this.COLS = top_evectors.getColumnDimension();
		this.number_of_conformations = COLS;
		this.number_of_pairs = ROWS;

		// this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "Distance_Pair_PCA" + File.separatorChar + model + File.separatorChar;
		// exist = new File(out_dir).exists();
		// if (!exist) success = (new File(out_dir)).mkdirs();
		// this.file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";
	}

	/* ******************************************************************************************************************************* */

	/**
	 * Method to calculate the matrix of delta vectors using the reference frame
	 */
	public void get_Distance_Pair_DV_Series()
	{
		Matrix ref_col = reference_distances;
		Matrix col = new Matrix(ROWS, 1);
		for (int b = 0; b < COLS; b++)
			{
				col = input_coords.getMatrix(0, ROWS - 1, b, b);
				Matrix delta = col.minus(ref_col);
				delta_vector_series.setMatrix(0, ROWS - 1, b, b, delta);
			}
		path = file_name_head + "_Delta_Vectors.txt.bz2";
		Matrix_IO.write_BZ2_Matrix(delta_vector_series, path, 9, 3);
		get_DVPs();
	}

	/**
	 * Computes the Distance-Pair delta vector projections and writes them to files: Un-normed, normed, weighted, weighted-normed
	 */
	private void get_DVPs()
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
		delta_vector_series = null;
		path = file_name_head + "_top_" + number_of_modes + "_DVPs_" + model + ".txt.bz2";
		Matrix_IO.write_Matrix_adaptive_spacing(projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_normed_DVPs_" + model + ".txt.bz2";
		Matrix_IO.write_Matrix_adaptive_spacing(normed_projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_DVPs_" + model + ".txt.bz2";
		Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_DVPs_" + model + ".txt.bz2";
		Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections, path);
	}

	public void setOut_dir(String out)
	{
		this.out_dir = out;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
		this.file_name_head = out_dir + "ss_" + number_of_pairs + "_Atom_Pairs";
	}


	/**
	 * @return the projections
	 */
	public Matrix getProjections()
	{
		return projections;
	}

	/**
	 * @return the normed_projections
	 */
	public Matrix getNormed_projections()
	{
		return normed_projections;
	}

	/**
	 * @return the weighted_projections
	 */
	public Matrix getWeighted_projections()
	{
		return weighted_projections;
	}

	/**
	 * @return the weighted_normed_projections
	 */
	public Matrix getWeighted_normed_projections()
	{
		return weighted_normed_projections;
	}
}
