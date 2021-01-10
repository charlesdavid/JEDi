package jedi;

import java.io.File;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;

import Jama.Matrix;
import support.Delta_Vector_Projector;
import support.Delta_Vector_Series;
import support.KMO_MSA;
import supportIO.Input_Parameters;
import supportPlot.PC_Plot;

public class JEDi_Do_Dist_Pairs
{
	final boolean doCORR, doPCORR, doSPARSIFY;
	boolean exist, success;
	final int number_of_pairs, number_of_modes_dist_pairs, number_of_conformations;
	int rank_cov;
	double trace_dist_COV, cond, det, shrinkage, KMO;
	final String OUT_DIR, type;
	String out, out_dir_DP, out_dir_SSA, out_dir_FES, out_dir_KPCA, path, line;
	List<Double> top_distance_eigenvalues_COV, top_distance_eigenvalues_CORR, top_distance_eigenvalues_PCORR, top_distance_eigenvalues_CORR_SPARSE,
			top_distance_eigenvalues_PCORR_SPARSE;
	Matrix msa, cov_dist, corr_dist, pcorr_dist;
	Matrix reference_distances, distance_matrix, delta_vector_series;;
	Matrix top_distance_evectors_COV, top_distance_evectors_CORR, top_distance_evectors_PCORR, top_distance_evectors_CORR_SPARSE, top_distance_evectors_PCORR_SPARSE;
	Matrix normed_projections_dist_COV, projections_dist_COV, weighted_normed_projections_dist_COV, weighted_projections_dist_COV;
	Matrix projections_dist_CORR, normed_projections_dist_CORR, weighted_projections_dist_CORR, weighted_normed_projections_dist_CORR;
	Matrix projections_dist_CORR_SPARSE, normed_projections_dist_CORR_SPARSE, weighted_projections_dist_CORR_SPARSE, weighted_normed_projections_dist_CORR_SPARSE;
	Matrix normed_projections_dist_PCORR, projections_dist_PCORR, weighted_projections_dist_PCORR, weighted_normed_projections_dist_PCORR;
	Matrix normed_projections_dist_PCORR_SPARSE, projections_dist_PCORR_SPARSE, weighted_projections_dist_PCORR_SPARSE, weighted_normed_projections_dist_PCORR_SPARSE;
	Matrix square_PCA_modes_COV, square_PCA_modes_CORR, square_PCA_modes_PCORR, weighted_Square_PCA_modes_COV, weighted_Square_PCA_modes_CORR, weighted_Square_PCA_modes_PCORR,
			weighted_Square_PCA_modes_CORR_SPARSE, weighted_Square_PCA_modes_PCORR_SPARSE;

	final NumberFormat nf3, nf6;
	final RoundingMode rm;

	/* ***************************************************** CONSTRUCTOR ********************************************************************** */

	public JEDi_Do_Dist_Pairs(Matrix distances, Matrix ref_distances, String outPCA, String outSSA, String outFES, String outKPCA)
	{
		super();

		this.OUT_DIR = Input_Parameters.OUT_DIR;
		this.doCORR = Input_Parameters.doCORR;
		this.doPCORR = Input_Parameters.doPCORR;
		this.doSPARSIFY = Input_Parameters.doSPARSIFY;
		this.number_of_modes_dist_pairs = Input_Parameters.MODES_DISTANCE_PAIRS;

		this.number_of_pairs = distances.getRowDimension();
		this.number_of_conformations = distances.getColumnDimension();
		this.type = "Distance_Pair_PCA";

		this.out_dir_DP = outPCA;
		this.out_dir_SSA = outSSA;
		this.out_dir_FES = outFES;
		this.out_dir_KPCA = outKPCA;

		this.distance_matrix = distances;
		this.reference_distances = ref_distances;


		this.rm = RoundingMode.HALF_UP;

		this.nf3 = NumberFormat.getInstance();
		this.nf3.setRoundingMode(rm);
		this.nf3.setMaximumFractionDigits(3);
		this.nf3.setMinimumFractionDigits(3);

		this.nf6 = NumberFormat.getInstance();
		this.nf6.setRoundingMode(rm);
		this.nf6.setMaximumFractionDigits(6);
		this.nf6.setMinimumFractionDigits(6);
	}

	/* ************************************************** DRIVER METHOD ************************************************************************* ******* */

	/**
	 * Method for performing dpPCA.
	 */
	public void do_Dist()
	{
		get_Distance_Pair_PCA();
		get_Distance_Pair_DVPs();
		do_Distance_Pair_SSA();
		if (Input_Parameters.doFES)
			{
				do_FES();
			}
		if (Input_Parameters.doKPCA)
			{
				do_KPCA();
			}
	}

	/* ******************************************************** METHODS ************************************************************************* ******* */

	private void get_Distance_Pair_PCA()
	{
		JEDi_Get_Distance_Pair_PCA dp_pca = new JEDi_Get_Distance_Pair_PCA(distance_matrix, number_of_pairs);

		dp_pca.set_Out_dir(out_dir_DP);
		dp_pca.get_Distance_Pair_PCA();

		shrinkage = dp_pca.getShrinkage();
		trace_dist_COV = dp_pca.getTrace_COV();
		cond = dp_pca.getCond_COV();
		det = dp_pca.getDet_COV();
		rank_cov = dp_pca.getRank_COV();
		top_distance_evectors_COV = dp_pca.getTop_evectors_dist_COV();
		top_distance_eigenvalues_COV = dp_pca.getTop_eigenvalues_COV();
		weighted_Square_PCA_modes_COV = dp_pca.getWeightedSquarePCAmodesCOV();

		corr_dist = dp_pca.getCORR_dist();
		pcorr_dist = dp_pca.getPcorr_dist();

		msa = KMO_MSA.get_MSA(corr_dist, pcorr_dist);
		KMO = KMO_MSA.getKMO();

		if (doCORR)
			{
				top_distance_evectors_CORR = dp_pca.getTop_evectors_dist_CORR();
				top_distance_eigenvalues_CORR = dp_pca.getTop_eigenvalues_CORR();
				weighted_Square_PCA_modes_CORR = dp_pca.getWeightedSquarePCAmodesCORR();
				if (doSPARSIFY)
					{
						top_distance_eigenvalues_CORR_SPARSE = dp_pca.getTop_eigenvalues_CORR_SPARSE();
						top_distance_evectors_CORR_SPARSE = dp_pca.getTop_evectors_dist_CORR_SPARSE();
						weighted_Square_PCA_modes_CORR_SPARSE = dp_pca.getWeighted_Square_PCA_modes_CORR_SPARSE();
					}
			}
		if (doPCORR)
			{
				top_distance_evectors_PCORR = dp_pca.getTop_evectors_dist_PCORR();
				top_distance_eigenvalues_PCORR = dp_pca.getTop_eigenvalues_PCORR();
				weighted_Square_PCA_modes_PCORR = dp_pca.getWeightedSquarePCAmodesPCORR();
				if (doSPARSIFY)
					{
						top_distance_eigenvalues_PCORR_SPARSE = dp_pca.getTop_eigenvalues_PCORR_SPARSE();
						top_distance_evectors_PCORR_SPARSE = dp_pca.getTop_evectors_dist_PCORR_SPARSE();
						weighted_Square_PCA_modes_PCORR_SPARSE = dp_pca.getWeighted_Square_PCA_modes_PCORR_SPARSE();
					}
			}
	}

	private void get_Distance_Pair_DVPs()
	{
		Delta_Vector_Series dvs = new Delta_Vector_Series(distance_matrix, reference_distances);
		delta_vector_series = dvs.get_Delta_Vector_Series();

		Delta_Vector_Projector dvp = new Delta_Vector_Projector(delta_vector_series, top_distance_evectors_COV, top_distance_eigenvalues_COV);
		dvp.get_DVPs();
		projections_dist_COV = dvp.get_Projections();
		normed_projections_dist_COV = dvp.get_Normed_Projections();
		weighted_projections_dist_COV = dvp.get_Weighted_Projections();
		weighted_normed_projections_dist_COV = dvp.get_Weighted_Normed_Projections();

		String out = out_dir_DP + "COV" + File.separatorChar;
		PC_Plot.createChart4Series(out, "Weighted_Projections_COV", weighted_projections_dist_COV);

		if (doCORR)
			{
				dvp = new Delta_Vector_Projector(delta_vector_series, top_distance_evectors_CORR, top_distance_eigenvalues_CORR);
				dvp.get_DVPs();
				projections_dist_CORR = dvp.get_Projections();
				normed_projections_dist_CORR = dvp.get_Normed_Projections();
				weighted_projections_dist_CORR = dvp.get_Weighted_Projections();
				weighted_normed_projections_dist_CORR = dvp.get_Weighted_Normed_Projections();

				out = out_dir_DP + "CORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_CORR", weighted_projections_dist_CORR);

				if (doSPARSIFY)
					{
						dvp = new Delta_Vector_Projector(delta_vector_series, top_distance_evectors_CORR_SPARSE, top_distance_eigenvalues_CORR_SPARSE);
						dvp.get_DVPs();
						projections_dist_CORR_SPARSE = dvp.get_Projections();
						normed_projections_dist_CORR_SPARSE = dvp.get_Normed_Projections();
						weighted_projections_dist_CORR_SPARSE = dvp.get_Weighted_Projections();
						weighted_normed_projections_dist_CORR_SPARSE = dvp.get_Weighted_Normed_Projections();

						out = out_dir_DP + "CORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_CORR_SPARSE", weighted_projections_dist_CORR_SPARSE);
					}

			}
		if (doPCORR)

			{
				dvp = new Delta_Vector_Projector(delta_vector_series, top_distance_evectors_PCORR, top_distance_eigenvalues_PCORR);
				dvp.get_DVPs();
				projections_dist_PCORR = dvp.get_Projections();
				normed_projections_dist_PCORR = dvp.get_Normed_Projections();
				weighted_projections_dist_PCORR = dvp.get_Weighted_Projections();
				weighted_normed_projections_dist_PCORR = dvp.get_Weighted_Normed_Projections();

				out = out_dir_DP + "PCORR" + File.separatorChar;
				PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR", weighted_projections_dist_PCORR);

				if (doSPARSIFY)
					{
						dvp = new Delta_Vector_Projector(delta_vector_series, top_distance_evectors_PCORR_SPARSE, top_distance_eigenvalues_PCORR_SPARSE);
						dvp.get_DVPs();
						projections_dist_PCORR_SPARSE = dvp.get_Projections();
						normed_projections_dist_PCORR_SPARSE = dvp.get_Normed_Projections();
						weighted_projections_dist_PCORR_SPARSE = dvp.get_Weighted_Projections();
						weighted_normed_projections_dist_PCORR_SPARSE = dvp.get_Weighted_Normed_Projections();

						out = out_dir_DP + "PCORR" + File.separatorChar + "sparse" + File.separatorChar;
						PC_Plot.createChart4Series(out, "Weighted_Projections_PCORR_SPARSE", weighted_projections_dist_PCORR_SPARSE);
					}
			}
	}

	public void do_FES()
	{
		JEDi_Get_FES fes = new JEDi_Get_FES(normed_projections_dist_COV, 0, 1, number_of_conformations, 0, 0);
		fes.set_Out_dir(out_dir_FES + "COV" + File.separatorChar);
		fes.get_FES();
		fes.write_FES_Log();

		if (doCORR)
			{
				fes = new JEDi_Get_FES(normed_projections_dist_CORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(out_dir_FES + "CORR" + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();

				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_dist_CORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(out_dir_FES + "CORR" + File.separatorChar + "sparse" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}

		if (doPCORR)
			{
				fes = new JEDi_Get_FES(normed_projections_dist_PCORR, 0, 1, number_of_conformations, 0, 0);
				fes.set_Out_dir(out_dir_FES + "PCORR" + File.separatorChar);
				fes.get_FES();
				fes.write_FES_Log();

				if (doSPARSIFY)
					{
						fes = new JEDi_Get_FES(normed_projections_dist_PCORR_SPARSE, 0, 1, number_of_conformations, 0, 0);
						fes.set_Out_dir(out_dir_FES + "PCORR" + File.separatorChar + "sparse" + File.separatorChar);
						fes.get_FES();
						fes.write_FES_Log();
					}
			}
	}

	public void do_KPCA()
	{
		JEDi_Do_Kernel_PCA kpca = new JEDi_Do_Kernel_PCA(normed_projections_dist_COV, out_dir_KPCA + "COV" + File.separatorChar);
		kpca.kPCA_Driver();

		if (doCORR)
			{
				kpca = new JEDi_Do_Kernel_PCA(normed_projections_dist_CORR, out_dir_KPCA + "CORR" + File.separatorChar);
				kpca.kPCA_Driver();

				if (doSPARSIFY)
					{
						kpca = new JEDi_Do_Kernel_PCA(normed_projections_dist_CORR_SPARSE, out_dir_KPCA + "CORR" + File.separatorChar + "sparse" + File.separatorChar);
						kpca.kPCA_Driver();
					}
			}

		if (doPCORR)
			{
				kpca = new JEDi_Do_Kernel_PCA(normed_projections_dist_PCORR, out_dir_KPCA + "PCORR" + File.separatorChar);
				kpca.kPCA_Driver();

				if (doSPARSIFY)
					{
						kpca = new JEDi_Do_Kernel_PCA(normed_projections_dist_PCORR_SPARSE, out_dir_KPCA + "PCORR" + File.separatorChar + "sparse" + File.separatorChar);
						kpca.kPCA_Driver();
					}
			}
	}

	private void do_Distance_Pair_SSA()
	{
		if (Input_Parameters.doCORR)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis("Distance_Pair_PCA", "COV_vs_CORR", top_distance_evectors_CORR, top_distance_evectors_COV);

				ssa.setOut_dir(out_dir_SSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis("Distance_Pair_PCA", "COV_vs_PCORR", top_distance_evectors_PCORR, top_distance_evectors_COV);

				ssa.setOut_dir(out_dir_SSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doPCORR)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis("Distance_Pair_PCA", "CORR_vs_PCORR", top_distance_evectors_PCORR, top_distance_evectors_CORR);

				ssa.setOut_dir(out_dir_SSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doCORR && Input_Parameters.doSPARSIFY)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis("Distance_Pair_PCA", "CORR_vs_CORR_SPARSE", top_distance_evectors_CORR_SPARSE,
						top_distance_evectors_CORR);

				ssa.setOut_dir(out_dir_SSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
		if (Input_Parameters.doPCORR && Input_Parameters.doSPARSIFY)
			{
				JEDi_Get_Subspace_Analysis ssa = new JEDi_Get_Subspace_Analysis("Distance_Pair_PCA", "PCORR_vs_PCORR_SPARSE", top_distance_evectors_PCORR_SPARSE,
						top_distance_evectors_PCORR);

				ssa.setOut_dir(out_dir_SSA);
				ssa.get_SSA();
				ssa.get_FSSA_Iterated();
			}
	}

	/* ************************************* SETTERS ************************************************************************* ** */

	public void setOut_dir_DP(String out_dir)
	{
		this.out_dir_DP = out_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}

	public void setOut_dir_SSA(String out_dir)
	{
		this.out_dir_SSA = out_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}

	public void setOut_dir_FES(String out_dir)
	{
		this.out_dir_FES = out_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}

	public void setOut_dir_KPCA(String out_dir)
	{
		this.out_dir_KPCA = out_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}

	/* ************************************* GETTERS ************************************************************************* ** */

	public int getRank_cov()
	{
		return rank_cov;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public double getShrinkage()
	{
		return shrinkage;
	}

	public double getCond_cov()
	{
		return cond;
	}

	public double getDet_cov()
	{
		return det;
	}

	public double getTrace_dist_COV()
	{
		return trace_dist_COV;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public List<Double> getTop_distance_eigenvalues_COV()
	{
		return top_distance_eigenvalues_COV;
	}

	public List<Double> getTop_distance_eigenvalues_CORR()
	{
		return top_distance_eigenvalues_CORR;
	}

	public List<Double> getTop_distance_eigenvalues_PCORR()
	{
		return top_distance_eigenvalues_PCORR;
	}

	// ----------------------------------------------------------------------------------------------------- //

	public Matrix getDistance_matrix()
	{
		return distance_matrix;
	}

	public Matrix getCov_dist()
	{
		return cov_dist;
	}

	public Matrix getCorr_dist()
	{
		return corr_dist;
	}

	public Matrix getPcorr_dist()
	{
		return pcorr_dist;
	}

	public Matrix getTop_distance_evectors_COV()
	{
		return top_distance_evectors_COV;
	}

	public Matrix getTop_distance_evectors_CORR()
	{
		return top_distance_evectors_CORR;
	}

	public Matrix getTop_distance_evectors_PCORR()
	{
		return top_distance_evectors_PCORR;
	}

	public Matrix getNormed_projections_dist_COV()
	{
		return normed_projections_dist_COV;
	}

	public Matrix getProjections_dist_COV()
	{
		return projections_dist_COV;
	}

	public Matrix getProjections_dist_CORR()
	{
		return projections_dist_CORR;
	}

	public Matrix getNormed_projections_dist_CORR()
	{
		return normed_projections_dist_CORR;
	}

	public Matrix getWeighted_normed_projections_dist_COV()
	{
		return weighted_normed_projections_dist_COV;
	}

	public Matrix getWeighted_normed_projections_dist_CORR()
	{
		return weighted_normed_projections_dist_CORR;
	}

	public Matrix getWeighted_projections_dist_COV()
	{
		return weighted_projections_dist_COV;
	}

	public Matrix getWeighted_projections_dist_CORR()
	{
		return weighted_projections_dist_CORR;
	}

	public Matrix getNormed_projections_dist_PCORR()
	{
		return normed_projections_dist_PCORR;
	}

	public Matrix getProjections_dist_PCORR()
	{
		return projections_dist_PCORR;
	}

	public Matrix getWeighted_projections_dist_PCORR()
	{
		return weighted_projections_dist_PCORR;
	}

	public Matrix getWeighted_normed_projections_dist_PCORR()
	{
		return weighted_normed_projections_dist_PCORR;
	}

	public Matrix getWeightedSquarePCAmodesCOV()
	{
		return weighted_Square_PCA_modes_COV;
	}

	public Matrix getWeightedSquarePCAmodesCORR()
	{
		return weighted_Square_PCA_modes_CORR;
	}

	public Matrix getWeightedSquarePCAmodesPCORR()
	{
		return weighted_Square_PCA_modes_PCORR;
	}

	public Matrix getMSA()
	{
		return msa;
	}

	public double getKMO()
	{
		return KMO;
	}

	public Matrix getNormed_projections_dist_CORR_SPARSE()
	{
		return normed_projections_dist_CORR_SPARSE;
	}

	public Matrix getNormed_projections_dist_PCORR_SPARSE()
	{
		return normed_projections_dist_PCORR_SPARSE;
	}

	public Matrix getWeighted_projections_dist_CORR_SPARSE()
	{
		return weighted_projections_dist_CORR_SPARSE;
	}

	public Matrix getWeighted_projections_dist_PCORR_SPARSE()
	{
		return weighted_projections_dist_PCORR_SPARSE;
	}

	public Matrix getWeighted_Square_PCA_modes_CORR()
	{
		return weighted_Square_PCA_modes_CORR;
	}

	public Matrix getWeighted_Square_PCA_modes_PCORR()
	{
		return weighted_Square_PCA_modes_PCORR;
	}

	public Matrix getWeighted_Square_PCA_modes_CORR_SPARSE()
	{
		return weighted_Square_PCA_modes_CORR_SPARSE;
	}

	public Matrix getWeighted_Square_PCA_modes_PCORR_SPARSE()
	{
		return weighted_Square_PCA_modes_PCORR_SPARSE;
	}
}
