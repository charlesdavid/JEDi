package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;
import support.Projector;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;

/**
 * JED class JED_Get_Cartesian_DVPs: Constructs the DVPs for the Cartesian subset. Note: DVPs are like Principle Components, except that the reference structure is not the mean
 * structure. Copyright (C) 2012 Dr. Charles David
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/license>
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Cartesian_DVPs
{
	boolean exist, success, verbose;
	String out_dir, out, type, model, file_name_head, path;
	final String AA = "All_Atom_cPCA", HA = "Heavy_Atom_cPCA", AC = "Alpha_Carbon_cPCA", BB = "Backbone_cPCA", Res = "Residue_cPCA", Q = "COV", R = "CORR", P = "PCORR";
	final int number_of_modes, number_of_conformations, ROWS, number_of_atoms;
	int number_of_residues;
	List<Double> eigenvalues;
	final Matrix reference_PDB_coordinates, input_coords, delta_vector_series, top_evectors;
	Matrix projections, normed_projections, weighted_normed_projections, weighted_projections;
	long startTime, endTime, totalTime;

	/**
	 * Constructor to create the delta vector series and the delta vector projections.
	 *
	 * @param data       The transformed Coordinates Matrix
	 * @param ref_coords The Reference coordinates from the reference PDB file
	 * @param evects     The top Cartesian eigenvectors
	 * @param evals      The top Cartesian eigenvalues
	 * @param dir        The working directory
	 * @param des        The job description
	 * @param type_pca   The type of PCA: All-Atom, Alpha-Carbon, Residue (AA, AC, Res)
	 * @param model_pca  The PCA model: COV, CORR, PCORR (Q, R, P)
	 */
	JEDi_Get_Cartesian_DVPs(Matrix data, Matrix ref_coords, Matrix evects, List<Double> evals, String type_pca, String model_pca)
	{

		this.input_coords = data;
		this.reference_PDB_coordinates = ref_coords;
		this.top_evectors = evects;
		this.eigenvalues = evals;
		this.type = type_pca;
		this.model = model_pca;

		this.verbose = Input_Parameters.verbose;

		this.number_of_modes = top_evectors.getColumnDimension();
		this.number_of_conformations = input_coords.getColumnDimension();
		this.ROWS = top_evectors.getRowDimension();
		this.number_of_atoms = (ROWS / 3);
		this.delta_vector_series = new Matrix(ROWS, number_of_conformations);
	}

	/**
	 * Computes the Cartesian delta vector series and writes it to file.
	 */
	public void get_Cartesian_DV_Series()
	{
		Matrix col = new Matrix(ROWS, 1);

		for (int i = 0; i < number_of_conformations; i++)
			{
				col = input_coords.getMatrix(0, ROWS - 1, i, i);
				Matrix delta = col.minus(reference_PDB_coordinates);
				delta_vector_series.setMatrix(0, ROWS - 1, i, i, delta);
			}
		get_DVPs();
	}

	/**
	 * Computes the Cartesian delta vector projections and writes them to files: Un-normed, normed, weighted, weighted-normed
	 */
	private void get_DVPs()
	{
		out = out_dir + model + File.separatorChar;
		create_Directory(out);

		projections = new Matrix(number_of_conformations, number_of_modes);
		normed_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_projections = new Matrix(number_of_conformations, number_of_modes);
		weighted_normed_projections = new Matrix(number_of_conformations, number_of_modes);

		for (int mode = 0; mode < number_of_modes; mode++)
			{
				double val = eigenvalues.get(mode);
				if (val < 0) val = 0;
				double weight = Math.sqrt(val); // Weight has units of Angstroms

				Matrix data1 = top_evectors.getMatrix(0, ROWS - 1, mode, mode);
				Matrix vector1 = Projector.get_Normed_arrayF(data1);

				for (int conf = 0; conf < number_of_conformations; conf++)
					{
						Matrix data2 = delta_vector_series.getMatrix(0, ROWS - 1, conf, conf);
						Matrix vector2 = Projector.get_Normed_arrayF(data2);

						double dp = Projector.get_InnerProduct(data1, data2); // dp has units of Angstroms
						double normed_dp = Projector.get_InnerProduct(vector1, vector2);  // normed_dp is unitless
						double w_dp = weight * dp;  // w_dp has units of Angstroms Squared
						double weighted_normed_dp = weight * normed_dp;  // Weight has units of Angstroms

						projections.set(conf, mode, dp);
						normed_projections.set(conf, mode, normed_dp);
						weighted_projections.set(conf, mode, w_dp);
						weighted_normed_projections.set(conf, mode, weighted_normed_dp);
					}
			}
		path = file_name_head + "_top_" + number_of_modes + "_DVPs_" + type + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_normed_DVPs_" + type + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(normed_projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_DVPs_" + type + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections, path);
		path = file_name_head + "_top_" + number_of_modes + "_weighted_normed_DVPs_" + type + ".txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections, path);
	}

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

	private void create_Directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
		file_name_head = out + "ss_" + number_of_residues + "_" + number_of_atoms;
	}

	// ********************************** SETTERS ***********************************************************

	public void set_Out_dir(String out)
	{
		this.out_dir = out;
		create_Directory(out_dir);
		file_name_head = out_dir + "ss_" + number_of_residues + "_" + number_of_atoms;
	}

	public void setNumber_of_residues(int numberOfResidues)
	{
		this.number_of_residues = numberOfResidues;
	}
}
