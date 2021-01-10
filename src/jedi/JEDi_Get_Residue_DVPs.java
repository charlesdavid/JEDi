package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;
import support.Projector;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;
import supportPlot.PC_Plot;

/**
 * JED class JED_Get_Cartesian_DVPs: Constructs the DVPs for the Cartesian subset. Note: DVPs are like Principle Components, except that the reference structure is not the mean
 * structure. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Residue_DVPs
{
	boolean exist, success, verbose;
	String out_dir, file_name_head, path;
	final String model;
	int Residue_Index;
	final int number_of_modes, number_of_conformations, ROWS, number_of_atoms;
	final double FLOOR = Input_Parameters.FLOOR;
	final List<Double> eigenvalues;
	final Matrix reference_PDB_coordinates, input_coords, delta_vector_series, top_evectors;
	Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;

	/* ************************************** CONSTRUCTOR ******************************************************************************** */

	/**
	 * Constructor to create the delta vector series and the delta vector projections.
	 *
	 * @param data       The transformed Coordinates Matrix
	 * @param ref_coords The Reference coordinates from the reference PDB file
	 * @param evects     The top Cartesian eigenvectors
	 * @param evals      The top Cartesian eigenvalues
	 * @param dir        The working directory
	 * @param des        The job description
	 * @param pca_model  The type of PCA: COV, CORR, PCORR (Q, R, P)
	 */
	JEDi_Get_Residue_DVPs(Matrix data, Matrix ref_coords, Matrix evects, List<Double> evals, String pca_model)
	{
		this.input_coords = data;
		this.reference_PDB_coordinates = ref_coords;
		this.top_evectors = evects;
		this.eigenvalues = evals;
		this.model = pca_model;

		this.number_of_modes = top_evectors.getColumnDimension();
		this.number_of_conformations = input_coords.getColumnDimension();
		this.ROWS = top_evectors.getRowDimension();
		this.number_of_atoms = (ROWS / 3);

		this.verbose = Input_Parameters.verbose;

		delta_vector_series = new Matrix(ROWS, number_of_conformations);
	}
	/* ************************************** PUBLIC METHODS ******************************************************************************** */

	/**
	 * Computes the Cartesian delta vector series and writes it to file.
	 */
	public Matrix get_Cartesian_DV_Series()
	{
		Matrix ref_col = reference_PDB_coordinates;
		Matrix col = new Matrix(ROWS, 1);
		for (int b = 0; b < number_of_conformations; b++)
			{
				col = input_coords.getMatrix(0, ROWS - 1, b, b);
				Matrix delta = col.minus(ref_col);
				delta_vector_series.setMatrix(0, ROWS - 1, b, b, delta);
			}
		get_DVPs();
		return delta_vector_series;
	}

	/* ************************************** PRIVATE METHODS ******************************************************************************** */

	/**
	 * Computes the Cartesian delta vector projections and writes them to files: Un-normed, normed, weighted, weighted-normed
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
						double weight = eigenvalues.get(outer);
						if (weight < FLOOR) weight = FLOOR;
						double sqrt_weight = Math.sqrt(weight);
						Matrix vector1 = Projector.get_Normed_arrayF(data1);
						Matrix vector2 = Projector.get_Normed_arrayF(data2);
						double dp = Projector.get_InnerProduct(data1, data2);
						double normed_dp = Projector.get_InnerProduct(vector1, vector2);
						double w_dp = sqrt_weight * dp;
						double weighted_normed_dp = sqrt_weight * normed_dp;
						projections.set(inner, outer, dp);
						normed_projections.set(inner, outer, normed_dp);
						weighted_projections.set(inner, outer, w_dp);
						weighted_normed_projections.set(inner, outer, weighted_normed_dp);
					}
			}
		String res_index = String.format("%03d", (Residue_Index));
		file_name_head = out_dir + "Residue_" + res_index;

		path = file_name_head + "_top_" + number_of_modes + "_DVPs_" + model + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(projections, path);

		path = file_name_head + "_top_" + number_of_modes + "_normed_DVPs_" + model + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(normed_projections, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_DVPs_" + model + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections, path);

		path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_DVPs_" + model + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections, path);

		PC_Plot.createChart4Series(file_name_head + "_", "Weighted_Projections_" + model, weighted_projections);
	}

	private void create_Directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();

	}

	/* ************************************** SETTERS ******************************************************************************** */

	public void set_Residue_Index(int residue_Index)
	{
		Residue_Index = residue_Index;
	}

	public void set_Output_Directory(String dir)
	{
		out_dir = dir;
		create_Directory(out_dir);
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

	public int get_Residue_Index()
	{
		return Residue_Index;
	}
}
