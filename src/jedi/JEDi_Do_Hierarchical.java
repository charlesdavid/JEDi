package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;
import support.Collectivity;
import support.Projector;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;

/**
 * JED class JED_Do_Hierarchical: Top class for implementing the Generalized Cartesian analysis. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Do_Hierarchical
{
	boolean exist, success, verbose;
	final int number_of_residues, number_of_eigenesidues, number_of_modes_Hierarchical, number_of_Conformations;
	int number_of_Atoms, rank;
	final double FLOOR = Input_Parameters.FLOOR;
	double trace, cond, det, shrinkage;
	String out_dir, file_name_head;
	double[] pca_mode_maxes, pca_mode_mins, convoluted_pca_mode_maxes, convoluted_pca_mode_mins;
	Matrix convoluted_Eigenvectors, convoluted_DVs, centered_data;
	Matrix top_evectors, square_pca_modes, weighted_square_pca_modes, weighted_pca_modes, pca_modes, convoluted_square_pca_modes, convoluted_weighted_square_pca_modes,
			convoluted_weighted_pca_modes, convoluted_pca_modes;
	Matrix normed_projections_rgc, projections_rgc, weighted_projections_rgc, weighted_normed_projections_rgc, normed_projections_pc, projections_pc, weighted_projections_pc,
			weighted_normed_projections_pc;
	final Matrix residuePCs;
	List<Double> top_eigenvalues, eigenvalues;
	final List<Matrix> residue_Eigenvectors, residue_DVs;
	JEDi_Get_FES fes;

	/* ****************************************************** CONSTRUCTOR ************************************************************************************** */

	public JEDi_Do_Hierarchical(int eigenresidues, int modes_hierarchical, Matrix resPCs, List<Matrix> Residue_Evects, List<Matrix> Residue_DVs)
	{
		this.number_of_eigenesidues = eigenresidues;
		this.residuePCs = resPCs;
		this.residue_Eigenvectors = Residue_Evects;
		this.residue_DVs = Residue_DVs;
		this.number_of_Conformations = residue_DVs.get(0).getColumnDimension();
		this.number_of_residues = (residuePCs.getRowDimension() / number_of_eigenesidues);
		this.number_of_modes_Hierarchical = modes_hierarchical;
		this.verbose = Input_Parameters.verbose;
	}

	/* *********************************************************** PUBLIC DRIVER METHOD *************************************************** */

	public void do_Hierarchical_PCA()
	{
		get_Hierarchical_PCA();
		get_Hierarchical_PCs();
		get_Convoluted_Eigenvectors();
		get_Convoluted_PCA_Modes();
		get_Convoluted_DVs();
		get_Convoluted_DVPs();
	}

	/* ************************************************************** PRIVATE METHODS *************************************************** */

	private void get_Hierarchical_PCA()
	{
		JEDi_Get_Hierarchical_PCA gcov_pca = new JEDi_Get_Hierarchical_PCA(residuePCs, number_of_eigenesidues, number_of_modes_Hierarchical);
		gcov_pca.setNumber_of_atoms(number_of_Atoms);
		gcov_pca.set_Out_Dir(out_dir);
		gcov_pca.get_Hierarchical_PCA();

		cond = gcov_pca.get_cond();
		trace = gcov_pca.get_trace();
		det = gcov_pca.get_det();
		rank = gcov_pca.get_rank();
		shrinkage = gcov_pca.get_shrinkage();

		eigenvalues = gcov_pca.getEigenvalues();
		top_eigenvalues = gcov_pca.getTop_eigenvalues();
		top_evectors = gcov_pca.getTop_evectors(); // This is the U matrix... will be convoluted with the V matrices (eigenresidues)

		pca_modes = gcov_pca.getPca_modes();
		square_pca_modes = gcov_pca.getSquare_pca_modes();
		weighted_square_pca_modes = gcov_pca.getWeighted_square_pca_modes();
		weighted_pca_modes = gcov_pca.getWeighted_pca_modes();

		pca_mode_mins = gcov_pca.get_pca_mode_mins();
		pca_mode_maxes = gcov_pca.get_pca_mode_maxes();

		centered_data = gcov_pca.getCentered_input_data();
	}

	private void get_Hierarchical_PCs()
	{
		JEDi_Get_PCs pcs_cov = new JEDi_Get_PCs(centered_data, top_evectors, eigenvalues);
		pcs_cov.get_PCs();

		projections_pc = pcs_cov.getPCs();
		normed_projections_pc = pcs_cov.getNormed_PCs();
		weighted_projections_pc = pcs_cov.getWeighted_PCs();
		weighted_normed_projections_pc = pcs_cov.getWeighted_normed_PCs();

		String path = file_name_head + "_Convoluted_PCs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(projections_pc, path);

		path = file_name_head + "_normed_Convoluted_PCs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(normed_projections_pc, path);

		path = file_name_head + "_weighted_Convoluted_PCs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections_pc, path);

		path = file_name_head + "_weighted_normed_Convoluted_PCs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections_pc, path);
	}

	private void get_Convoluted_Eigenvectors()
	{
		JEDi_Get_Convoluted_Eigenvectors getUV = new JEDi_Get_Convoluted_Eigenvectors(top_evectors, residue_Eigenvectors, number_of_eigenesidues);

		convoluted_Eigenvectors = getUV.get_Convoluted_Eigenvectors();

		String path = file_name_head + "_Convoluted_Eigenvectors.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix(convoluted_Eigenvectors, path, 15, 12);
	}

	private void get_Convoluted_PCA_Modes()
	{
		JEDi_Get_PCA_Modes gPCA = new JEDi_Get_PCA_Modes(convoluted_Eigenvectors, eigenvalues);

		convoluted_pca_modes = gPCA.get_PCA_modes();
		convoluted_weighted_pca_modes = gPCA.get_Weighted_PCA_modes();
		convoluted_square_pca_modes = gPCA.get_Square_PCA_modes();
		convoluted_weighted_square_pca_modes = gPCA.get_Weighted_Square_PCA_modes();
		convoluted_pca_mode_maxes = gPCA.get_PCA_mode_maxs();
		convoluted_pca_mode_mins = gPCA.get_PCA_mode_mins();

		Matrix collectivity = Collectivity.get_Collectivity(square_pca_modes);
		String path = file_name_head + "_Eigenvector_Collectivity.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(collectivity, path);

		path = file_name_head + "_Convoluted_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(convoluted_pca_modes, path);
		path = file_name_head + "_Convoluted_weighted_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(convoluted_weighted_pca_modes, path);
		path = file_name_head + "_Convoluted_square_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(convoluted_square_pca_modes, path);
		path = file_name_head + "_Convoluted_weighted_square_pca_modes.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(convoluted_weighted_square_pca_modes, path);
	}

	private void get_Convoluted_DVs()
	{
		convoluted_DVs = new Matrix(number_of_Atoms * 3, number_of_Conformations);

		int big_DV_offset = 0;
		for (int i = 0; i < number_of_residues; i++)
			{
				Matrix dv = residue_DVs.get(i);
				int num_of_atoms_Res = (dv.getRowDimension() / 3);
				for (int j = 0; j < number_of_Conformations; j++)
					{
						for (int k = 0; k < num_of_atoms_Res; k++)
							{
								double X = dv.get(k, j);
								double Y = dv.get(k + num_of_atoms_Res, j);
								double Z = dv.get(k + 2 * num_of_atoms_Res, j);

								convoluted_DVs.set(k + big_DV_offset, j, X);
								convoluted_DVs.set(k + big_DV_offset + number_of_Atoms, j, Y);
								convoluted_DVs.set(k + big_DV_offset + 2 * number_of_Atoms, j, Z);
							}
					}
				big_DV_offset += num_of_atoms_Res;
			}
	}

	private void get_Convoluted_DVPs()
	{
		projections_rgc = new Matrix(number_of_Conformations, number_of_modes_Hierarchical);
		weighted_projections_rgc = new Matrix(number_of_Conformations, number_of_modes_Hierarchical);
		normed_projections_rgc = new Matrix(number_of_Conformations, number_of_modes_Hierarchical);
		weighted_normed_projections_rgc = new Matrix(number_of_Conformations, number_of_modes_Hierarchical);

		for (int mode = 0; mode < number_of_modes_Hierarchical; mode++)
			{
				double val = eigenvalues.get(mode);
				if (val < FLOOR) val = FLOOR;
				double weight = Math.sqrt(val); // Weight has units of Angstroms

				int row_index_1 = 0;
				int row_index_2 = (3 * (number_of_Atoms - 1));

				Matrix data1 = convoluted_Eigenvectors.getMatrix(row_index_1, row_index_2, mode, mode);
				Matrix vector1 = Projector.get_Normed_arrayF(data1);

				for (int conf = 0; conf < number_of_Conformations; conf++)
					{
						Matrix data2 = convoluted_DVs.getMatrix(row_index_1, row_index_2, conf, conf);
						Matrix vector2 = Projector.get_Normed_arrayF(data2);

						double dp = Projector.get_InnerProduct(data1, data2); // dp has units of Angstroms
						double normed_dp = Projector.get_InnerProduct(vector1, vector2); // normed_dp is unitless
						double w_dp = weight * dp; // w_dp has units of Angstroms Squared
						double weighted_normed_dp = weight * normed_dp; // w_dp has units of Angstroms

						projections_rgc.set(conf, mode, dp);
						normed_projections_rgc.set(conf, mode, normed_dp);
						weighted_projections_rgc.set(conf, mode, w_dp);
						weighted_normed_projections_rgc.set(conf, mode, weighted_normed_dp);
					}
			}
		String path = file_name_head + "_Convoluted_DVPs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(projections_rgc, path);

		path = file_name_head + "_normed_Convoluted_DVPs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(normed_projections_rgc, path);

		path = file_name_head + "_weighted_Convoluted_DVPs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections_rgc, path);

		path = file_name_head + "_weighted_normed_Convoluted_DVPs.txt.bz2";
		if (verbose) Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections_rgc, path);
	}

	public void do_FES(String out)
	{
		fes = new JEDi_Get_FES(normed_projections_rgc, 0, 1, number_of_Conformations, 0, 0);
		fes.set_Out_dir(out);
		fes.get_FES();
		fes.write_FES_Log();
	}

	public void do_KPCA(String out)
	{
		JEDi_Do_Kernel_PCA kpca = new JEDi_Do_Kernel_PCA(normed_projections_rgc, out);
		kpca.kPCA_Driver();
	}

	/* ************************************************************** SETTERS ******************************************************************* */

	public void set_Out_Dir(String dir) // be sure that the number of atoms is set first!
	{
		this.out_dir = dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();

		file_name_head = out_dir + "ss_" + number_of_residues + "_" + number_of_Atoms + "_top_" + number_of_modes_Hierarchical;
	}

	public void setNumber_of_Atoms(int numberOfAtoms)
	{
		this.number_of_Atoms = numberOfAtoms;
	}

	/* ************************************************************** GETTERS ******************************************************************* */

	public double get_Shrinkage()
	{
		return shrinkage;
	}

	public double get_Trace()
	{
		return trace;
	}

	public double get_cond()
	{
		return cond;
	}

	public double get_Cond()
	{
		return cond;
	}

	public double get_Det()
	{
		return det;
	}

	public int get_Rank()
	{
		return rank;
	}

	public List<Double> getTop_eigenvalues()
	{
		return top_eigenvalues;
	}

	public Matrix getTop_evectors()
	{
		return top_evectors;
	}

	public Matrix getProjections()
	{
		return projections_rgc;
	}

	public Matrix getWeighted_normed_projections()
	{
		return weighted_normed_projections_rgc;
	}

	public int getNumber_of_residues()
	{
		return number_of_residues;
	}

	public double getTrace()
	{
		return trace;
	}

	public double getDet()
	{
		return det;
	}

	public double getRank()
	{
		return rank;
	}

	public List<Double> getEigenvalues()
	{
		return eigenvalues;
	}

	public Matrix getNormed_projections_rgc()
	{
		return normed_projections_rgc;
	}

	public Matrix getProjections_rgc()
	{
		return projections_rgc;
	}

	public Matrix getWeighted_projections_rgc()
	{
		return weighted_projections_rgc;
	}

	public Matrix getWeighted_normed_projections_rgc()
	{
		return weighted_normed_projections_rgc;
	}

	public Matrix getNormed_projections_pc()
	{
		return normed_projections_pc;
	}

	public Matrix getProjections_pc()
	{
		return projections_pc;
	}

	public Matrix getWeighted_projections_pc()
	{
		return weighted_projections_pc;
	}

	public Matrix getWeighted_normed_projections_pc()
	{
		return weighted_normed_projections_pc;
	}

	public double[] getConvoluted_pca_mode_max()
	{
		return convoluted_pca_mode_maxes;
	}

	public double[] getConvoluted_pca_mode_min()
	{
		return convoluted_pca_mode_mins;
	}

	public Matrix getConvoluted_Eigenvectors()
	{
		return convoluted_Eigenvectors;
	}

	public Matrix getConvoluted_DVs()
	{
		return convoluted_DVs;
	}

	public Matrix getConvoluted_square_pca_modes()
	{
		return convoluted_square_pca_modes;
	}

	public Matrix getConvoluted_weighted_square_pca_modes()
	{
		return convoluted_weighted_square_pca_modes;
	}

	public Matrix getConvoluted_weighted_pca_modes()
	{
		return convoluted_weighted_pca_modes;
	}

	public Matrix getConvoluted_pca_modes()
	{
		return convoluted_pca_modes;
	}
}
