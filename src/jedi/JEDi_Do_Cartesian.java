package jedi;

import java.io.File;
import java.util.List;

import Jama.Matrix;
import supportIO.Input_Parameters;

/**
 * JED class JED_Do_Cartesian: Top class for implementing the Cartesian analysis. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Do_Cartesian
{
	boolean exist, success;
	final boolean doCORR, doPCORR, doSPARSIFY;
	final String DIRECTORY, DESCRIPTION, typePCA, Q = "COV", R = "CORR", P = "PCORR", RS = "CORR_SPARSE", PS = "PCORR_SPARSE";
	int number_of_residues, rank_COV;
	String out_dir;
	final int number_of_modes_SS, number_of_conformations;
	double shrinkage, trace_COV, cond_COV, det_COV;
	double[] pca_mode_max_COV, pca_mode_min_COV, pca_mode_max_CORR, pca_mode_min_CORR, pca_mode_max_PCORR, pca_mode_min_PCORR, pca_mode_max_CORR_SPARSE, pca_mode_min_CORR_SPARSE,
			pca_mode_max_PCORR_SPARSE, pca_mode_min_PCORR_SPARSE;
	List<Double> top_cartesian_eigenvalues_COV, top_cartesian_eigenvalues_CORR, top_cartesian_eigenvalues_PCORR, eigenvalues_COV, eigenvalues_CORR, eigenvalues_PCORR;
	List<Double> top_cartesian_eigenvalues_CORR_SPARSE, top_cartesian_eigenvalues_PCORR_SPARSE;
	final Matrix aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates, adjusted_subset_PDB_coordinates;
	Matrix cov, corr, pcorr;
	Matrix top_cartesian_evectors_COV, top_cartesian_evectors_CORR, top_cartesian_evectors_PCORR, square_pca_modes_COV, weighted_square_pca_modes_COV, weighted_pca_modes_COV,
			square_pca_modes_CORR, weighted_square_pca_modes_CORR, weighted_pca_modes_CORR, pca_modes_COV, pca_modes_CORR, square_pca_modes_PCORR, weighted_square_pca_modes_PCORR,
			weighted_pca_modes_PCORR, pca_modes_PCORR;
	Matrix normed_projections_COV, normed_projections_CORR, projections_COV, projections_CORR, normed_projections_PCORR, projections_PCORR, weighted_projections_COV,
			weighted_projections_CORR, weighted_projections_PCORR, weighted_normed_projections_COV, weighted_normed_projections_CORR, weighted_normed_projections_PCORR;
	Matrix top_cartesian_evectors_CORR_SPARSE, top_cartesian_evectors_PCORR_SPARSE, square_pca_modes_CORR_SPARSE, square_pca_modes_PCORR_SPARSE, normed_projections_CORR_SPARSE,
			weighted_projections_CORR_SPARSE, normed_projections_PCORR_SPARSE, weighted_projections_PCORR_SPARSE, weighted_square_pca_modes_CORR_SPARSE,
			weighted_square_pca_modes_PCORR_SPARSE;

	/* ******************************************************************* CONSTRUCTOR **************************************************************** */

	/**
	 * Driver for the Direct Cartesian Analyses:
	 *
	 * @param number_of_modes     The number of Cartesian PCA modes to compute.
	 * @param PDB_coordinates     The Coordinates Matrix
	 * @param Ref_PDB_coordinates The Reference Coordinates Matrix
	 * @param type_pca            The type of PCA, e.g., all atom, heavy atom, etc.
	 */
	public JEDi_Do_Cartesian(int number_of_modes, Matrix Aligned_PDB_coordinates, Matrix Aligned_REF_PDB_coords, Matrix Adjusted_PDB_Coords, String type_pca)
	{
		this.number_of_modes_SS = number_of_modes;
		this.aligned_subset_PDB_coordinates = Aligned_PDB_coordinates;
		this.aligned_subset_REF_PDB_coordinates = Aligned_REF_PDB_coords;
		this.adjusted_subset_PDB_coordinates = Adjusted_PDB_Coords;
		this.typePCA = type_pca;
		this.number_of_conformations = aligned_subset_PDB_coordinates.getColumnDimension();

		this.DIRECTORY = Input_Parameters.DIRECTORY;
		this.DESCRIPTION = Input_Parameters.DESCRIPTION;
		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doSPARSIFY = Input_Parameters.doSPARSIFY;
	}
	/* ************************************************** DRIVER METHODS ******************************************************************************** */

	public void do_Cartesian()
	{
		get_Cartesian_PCA();
		get_Cartesian_DVPs();
	}

	/* ****************************************************** METHODS *********************************************************************************** */

	private void get_Cartesian_PCA()
	{
		JEDi_Get_Cartesian_PCA c_pca = new JEDi_Get_Cartesian_PCA(adjusted_subset_PDB_coordinates, number_of_modes_SS, typePCA);

		c_pca.setNumber_of_residues(number_of_residues);
		c_pca.set_Out_dir(out_dir);
		c_pca.Do_Cov_PCA();

		cond_COV = c_pca.get_cond_COV();
		trace_COV = c_pca.get_trace_COV();
		det_COV = c_pca.get_det_COV();
		rank_COV = c_pca.get_rank_COV();
		shrinkage = c_pca.getShrinkage();

		cov = c_pca.cov;
		corr = c_pca.corr;
		pcorr = c_pca.pcorr;

		eigenvalues_COV = c_pca.getEigenvalues_COV();
		top_cartesian_evectors_COV = c_pca.getTop_evectors_COV();
		top_cartesian_eigenvalues_COV = c_pca.getTop_eigenvalues_COV();
		pca_modes_COV = c_pca.getPca_modes_COV();
		square_pca_modes_COV = c_pca.getSquare_pca_modes_COV();
		weighted_square_pca_modes_COV = c_pca.getWeighted_square_pca_modes_COV();
		weighted_pca_modes_COV = c_pca.getWeighted_pca_modes_COV();
		pca_mode_min_COV = c_pca.get_pca_mode_mins_COV();
		pca_mode_max_COV = c_pca.get_pca_mode_maxes_COV();

		if (doCORR)
			{
				c_pca.Do_Corr_PCA();

				eigenvalues_CORR = c_pca.getEigenvalues_CORR();
				top_cartesian_evectors_CORR = c_pca.getTop_evectors_CORR();
				top_cartesian_eigenvalues_CORR = c_pca.getTop_eigenvalues_CORR();
				pca_modes_CORR = c_pca.getPca_modes_CORR();
				square_pca_modes_CORR = c_pca.getSquare_pca_modes_CORR();
				weighted_square_pca_modes_CORR = c_pca.getWeighted_square_pca_modes_CORR();
				weighted_pca_modes_CORR = c_pca.getWeighted_pca_modes_CORR();
				pca_mode_min_CORR = c_pca.get_pca_mode_mins_CORR();
				pca_mode_max_CORR = c_pca.get_pca_mode_maxes_CORR();

				if (doSPARSIFY)
					{
						c_pca.Do_CORR_SPARSE_PCA();

						top_cartesian_eigenvalues_CORR_SPARSE = c_pca.getTop_eigenvalues_CORR_SPARSE();
						top_cartesian_evectors_CORR_SPARSE = c_pca.getTop_evectors_CORR_SPARSE();
						square_pca_modes_CORR_SPARSE = c_pca.getSquare_pca_modes_CORR_SPARSE();
						weighted_square_pca_modes_CORR_SPARSE = c_pca.getWeighted_square_pca_modes_CORR_SPARSE();
						pca_mode_min_CORR_SPARSE = c_pca.get_pca_mode_mins_CORR_SPARSE();
						pca_mode_max_CORR_SPARSE = c_pca.get_pca_mode_maxes_CORR_SPARSE();
					}
			}

		if (doPCORR)
			{
				c_pca.Do_PCorr_PCA();

				eigenvalues_PCORR = c_pca.getEigenvalues_PCORR();
				top_cartesian_evectors_PCORR = c_pca.getTop_evectors_PCORR();
				top_cartesian_eigenvalues_PCORR = c_pca.getEigenvalues_PCORR();
				pca_modes_PCORR = c_pca.getPca_modes_PCORR();
				square_pca_modes_PCORR = c_pca.getSquare_pca_modes_PCORR();
				weighted_square_pca_modes_PCORR = c_pca.getWeighted_square_pca_modes_PCORR();
				weighted_pca_modes_PCORR = c_pca.getWeighted_pca_modes_PCORR();
				pca_mode_min_PCORR = c_pca.get_pca_mode_mins_PCORR();
				pca_mode_max_PCORR = c_pca.get_pca_mode_maxes_PCORR();

				if (doSPARSIFY)
					{
						c_pca.Do_PCORR_SPARSE_PCA();

						top_cartesian_eigenvalues_PCORR_SPARSE = c_pca.getTop_eigenvalues_PCORR_SPARSE();
						top_cartesian_evectors_PCORR_SPARSE = c_pca.getTop_evectors_PCORR_SPARSE();
						square_pca_modes_PCORR_SPARSE = c_pca.getSquare_pca_modes_PCORR_SPARSE();
						weighted_square_pca_modes_PCORR_SPARSE = c_pca.getWeighted_square_pca_modes_PCORR_SPARSE();
						pca_mode_min_PCORR_SPARSE = c_pca.get_pca_mode_mins_PCORR_SPARSE();
						pca_mode_max_PCORR_SPARSE = c_pca.get_pca_mode_maxes_PCORR_SPARSE();
					}
			}
	}

	private void get_Cartesian_DVPs()
	{
		JEDi_Get_Cartesian_DVPs pcs_cov = new JEDi_Get_Cartesian_DVPs(aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates, top_cartesian_evectors_COV,
				eigenvalues_COV, typePCA, Q);

		pcs_cov.set_Out_dir(out_dir);
		pcs_cov.setNumber_of_residues(number_of_residues);
		pcs_cov.get_Cartesian_DV_Series();
		projections_COV = pcs_cov.getProjections();
		normed_projections_COV = pcs_cov.getNormed_projections();
		weighted_projections_COV = pcs_cov.getWeighted_projections();
		weighted_normed_projections_COV = pcs_cov.getWeighted_normed_projections();

		if (doCORR)
			{
				JEDi_Get_Cartesian_DVPs pcs_corr = new JEDi_Get_Cartesian_DVPs(aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates, top_cartesian_evectors_CORR,
						top_cartesian_eigenvalues_CORR, typePCA, R);

				pcs_corr.set_Out_dir(out_dir);
				pcs_corr.setNumber_of_residues(number_of_residues);
				pcs_corr.get_Cartesian_DV_Series();
				projections_CORR = pcs_corr.getProjections();
				normed_projections_CORR = pcs_corr.getNormed_projections();
				weighted_projections_CORR = pcs_corr.getWeighted_projections();
				weighted_normed_projections_CORR = pcs_corr.getWeighted_normed_projections();

				if (doSPARSIFY)
					{
						JEDi_Get_Cartesian_DVPs pcs_corrS = new JEDi_Get_Cartesian_DVPs(aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates,
								top_cartesian_evectors_CORR_SPARSE, top_cartesian_eigenvalues_CORR_SPARSE, typePCA, R + File.separatorChar + "sparse");

						pcs_corrS.set_Out_dir(out_dir);
						pcs_corrS.setNumber_of_residues(number_of_residues);
						pcs_corrS.get_Cartesian_DV_Series();
						normed_projections_CORR_SPARSE = pcs_corr.getNormed_projections();
						weighted_projections_CORR_SPARSE = pcs_corrS.getWeighted_projections();
					}
			}

		if (doPCORR)
			{
				JEDi_Get_Cartesian_DVPs pcs_pcorr = new JEDi_Get_Cartesian_DVPs(aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates, top_cartesian_evectors_PCORR,
						top_cartesian_eigenvalues_PCORR, typePCA, P);

				pcs_pcorr.set_Out_dir(out_dir);
				pcs_pcorr.setNumber_of_residues(number_of_residues);
				pcs_pcorr.get_Cartesian_DV_Series();
				projections_PCORR = pcs_pcorr.getProjections();
				normed_projections_PCORR = pcs_pcorr.getNormed_projections();
				weighted_projections_PCORR = pcs_pcorr.getWeighted_projections();
				weighted_normed_projections_PCORR = pcs_pcorr.getWeighted_normed_projections();

				if (doSPARSIFY)
					{
						JEDi_Get_Cartesian_DVPs pcs_pcorrS = new JEDi_Get_Cartesian_DVPs(aligned_subset_PDB_coordinates, aligned_subset_REF_PDB_coordinates,
								top_cartesian_evectors_PCORR_SPARSE, top_cartesian_eigenvalues_PCORR_SPARSE, typePCA, P + File.separatorChar + "sparse");

						pcs_pcorrS.set_Out_dir(out_dir);
						pcs_pcorrS.setNumber_of_residues(number_of_residues);
						pcs_pcorrS.get_Cartesian_DV_Series();
						normed_projections_PCORR_SPARSE = pcs_pcorr.getNormed_projections();
						weighted_projections_PCORR_SPARSE = pcs_pcorr.getWeighted_projections();
					}
			}
	}

	public void do_FES(String out)
	{
		JEDi_Get_FES fes = new JEDi_Get_FES(normed_projections_COV, 0, 1, number_of_conformations, 0, 0);
		fes.set_Out_dir(out + Q + File.separatorChar);
		fes.get_FES();
		fes.write_FES_Log();

		if (doCORR)
			{
				fes = new JEDi_Get_FES(normed_projections_CORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(out + R + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();
				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_CORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(out + R + File.separatorChar + "sparse" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}

		if (doPCORR)
			{
				fes = new JEDi_Get_FES(normed_projections_PCORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(out + P + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();
				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_PCORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(out + P + File.separatorChar + "sparse" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}
	}

	public void do_KPCA(String outKPCA)
	{
		JEDi_Do_Kernel_PCA kpca_COV = new JEDi_Do_Kernel_PCA(normed_projections_COV, outKPCA + Q + File.separatorChar);
		kpca_COV.kPCA_Driver();

		if (doCORR)
			{
				JEDi_Do_Kernel_PCA kpca_CORR = new JEDi_Do_Kernel_PCA(normed_projections_CORR, outKPCA + R + File.separatorChar);
				kpca_CORR.kPCA_Driver();
				if (doSPARSIFY)
					{
						kpca_CORR = new JEDi_Do_Kernel_PCA(normed_projections_CORR_SPARSE, outKPCA + R + File.separatorChar + "sparse" + File.separatorChar);
						kpca_CORR.kPCA_Driver();
					}
			}

		if (doPCORR)
			{
				JEDi_Do_Kernel_PCA kpca_PCORR = new JEDi_Do_Kernel_PCA(normed_projections_PCORR, outKPCA + P + File.separatorChar);
				kpca_PCORR.kPCA_Driver();
				if (doSPARSIFY)
					{
						kpca_PCORR = new JEDi_Do_Kernel_PCA(normed_projections_PCORR_SPARSE, outKPCA + P + File.separatorChar + "sparse" + File.separatorChar);
						kpca_PCORR.kPCA_Driver();
					}
			}
	}

	/* ************************************************************ SETTERS ********************************************************************* */

	public void set_Out_dir(String out)
	{
		this.out_dir = out;
		exist = new File(out).exists();
		if (!exist) success = (new File(out)).mkdirs();
	}

	public void setNumber_of_residues(int numberOfResidues)
	{
		this.number_of_residues = numberOfResidues;
	}

	/* ************************************************************** GETTERS ******************************************************************* */

	public double get_Shrinkage()
	{
		return shrinkage;
	}

	public double get_Cond_COV()
	{
		return cond_COV;
	}

	public double get_Det_COV()
	{
		return det_COV;
	}

	public int get_Rank_COV()
	{
		return rank_COV;
	}

	public double get_Trace_COV()
	{
		return trace_COV;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public List<Double> getTop_cartesian_eigenvalues_COV()
	{

		return top_cartesian_eigenvalues_COV;
	}

	public List<Double> getTop_cartesian_eigenvalues_CORR()
	{

		return top_cartesian_eigenvalues_CORR;
	}

	public List<Double> getTop_cartesian_eigenvalues_PCORR()
	{
		return top_cartesian_eigenvalues_PCORR;
	}

	public List<Double> getEigenvalues_COV()
	{
		return eigenvalues_COV;
	}

	public List<Double> getEigenvalues_CORR()
	{
		return eigenvalues_CORR;
	}

	public List<Double> getEigenvalues_PCORR()
	{
		return eigenvalues_PCORR;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public double[] getPca_mode_max_COV()
	{
		return pca_mode_max_COV;
	}

	public double[] getPca_mode_min_COV()
	{
		return pca_mode_min_COV;
	}

	public double[] getPca_mode_max_CORR()
	{
		return pca_mode_max_CORR;
	}

	public double[] getPca_mode_min_CORR()
	{
		return pca_mode_min_CORR;
	}

	public double[] getPca_mode_max_PCORR()
	{
		return pca_mode_max_PCORR;
	}

	public double[] getPca_mode_min_PCORR()
	{
		return pca_mode_min_PCORR;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public Matrix getCov()
	{
		return cov;
	}

	public Matrix getCorr()
	{

		return corr;
	}

	public Matrix getPcorr()
	{

		return pcorr;
	}

	public Matrix getTop_cartesian_evectors_COV()
	{

		return top_cartesian_evectors_COV;
	}

	public Matrix getTop_cartesian_evectors_CORR()
	{
		return top_cartesian_evectors_CORR;
	}

	public Matrix getTop_cartesian_evectors_PCORR()
	{
		return top_cartesian_evectors_PCORR;
	}

	public Matrix getWeighted_pca_modes_COV()
	{

		return weighted_pca_modes_COV;
	}

	public Matrix getSquare_pca_modes_COV()
	{

		return square_pca_modes_COV;
	}

	public Matrix getWeighted_square_pca_modes_COV()
	{

		return weighted_square_pca_modes_COV;
	}

	public Matrix getNormed_projections_COV()
	{

		return normed_projections_COV;
	}

	public Matrix getProjections_COV()
	{

		return projections_COV;
	}

	public Matrix getWeighted_pca_modes_CORR()
	{

		return weighted_pca_modes_CORR;
	}

	public Matrix getSquare_pca_modes_CORR()
	{
		return square_pca_modes_CORR;
	}

	public Matrix getWeighted_square_pca_modes_CORR()
	{
		return weighted_square_pca_modes_CORR;
	}

	public Matrix getWeighted_square_pca_modes_PCORR()
	{
		return weighted_square_pca_modes_PCORR;
	}

	public Matrix getNormed_projections_CORR()
	{
		return normed_projections_CORR;
	}

	public Matrix getProjections_CORR()
	{
		return projections_CORR;
	}

	public Matrix getNormed_projections_PCORR()
	{
		return normed_projections_PCORR;
	}

	public Matrix getProjections_PCORR()
	{
		return projections_PCORR;
	}

	public Matrix getWeighted_projections_COV()
	{
		return weighted_projections_COV;
	}

	public Matrix getWeighted_projections_CORR()
	{
		return weighted_projections_CORR;
	}

	public Matrix getWeighted_projections_PCORR()
	{
		return weighted_projections_PCORR;
	}

	public Matrix getWeighted_normed_projections_COV()
	{
		return weighted_normed_projections_COV;
	}

	public Matrix getWeighted_normed_projections_CORR()
	{
		return weighted_normed_projections_CORR;
	}

	public Matrix getWeighted_normed_projections_PCORR()
	{
		return weighted_normed_projections_PCORR;
	}

	public Matrix getPca_modes_COV()
	{
		return pca_modes_COV;
	}

	public Matrix getPca_modes_CORR()
	{
		return pca_modes_CORR;
	}

	public Matrix getSquare_pca_modes_PCORR()
	{
		return square_pca_modes_PCORR;
	}

	public double[] getPca_mode_max_CORR_SPARSE()
	{
		return pca_mode_max_CORR_SPARSE;
	}

	public double[] getPca_mode_min_CORR_SPARSE()
	{
		return pca_mode_min_CORR_SPARSE;
	}

	public double[] getPca_mode_max_PCORR_SPARSE()
	{
		return pca_mode_max_PCORR_SPARSE;
	}

	public double[] getPca_mode_min_PCORR_SPARSE()
	{
		return pca_mode_min_PCORR_SPARSE;
	}

	public List<Double> getTop_cartesian_eigenvalues_CORR_SPARSE()
	{
		return top_cartesian_eigenvalues_CORR_SPARSE;
	}

	public List<Double> getTop_cartesian_eigenvalues_PCORR_SPARSE()
	{
		return top_cartesian_eigenvalues_PCORR_SPARSE;
	}

	public Matrix getTop_cartesian_evectors_CORR_SPARSE()
	{
		return top_cartesian_evectors_CORR_SPARSE;
	}

	public Matrix getTop_cartesian_evectors_PCORR_SPARSE()
	{
		return top_cartesian_evectors_PCORR_SPARSE;
	}

	public Matrix getSquare_pca_modes_CORR_SPARSE()
	{
		return square_pca_modes_CORR_SPARSE;
	}

	public Matrix getSquare_pca_modes_PCORR_SPARSE()
	{
		return square_pca_modes_PCORR_SPARSE;
	}

	public Matrix getWeighted_projections_CORR_SPARSE()
	{
		return weighted_projections_CORR_SPARSE;
	}

	public Matrix getWeighted_projections_PCORR_SPARSE()
	{
		return weighted_projections_PCORR_SPARSE;
	}

	public Matrix getWeighted_square_pca_modes_CORR_SPARSE()
	{
		return weighted_square_pca_modes_CORR_SPARSE;
	}

	public Matrix getWeighted_square_pca_modes_PCORR_SPARSE()
	{
		return weighted_square_pca_modes_PCORR_SPARSE;
	}
}
