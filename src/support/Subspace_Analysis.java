package support;

import java.util.ArrayList;

import Jama.Matrix;

/**
 * Subspace_Analysis: Class for performing comparative subspace analysis.
 * 
 * Note: The expected input is 2 matrices of eigenvectors representing 2 equidimensional subspaces from a vector space.
 * 
 * Equidimensional Orthonormal bases are expected:
 * 
 * This means that the number of rows must match, and number of columns must match.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Subspace_Analysis
{
	final int ROWS1, COLS1, ROWS2, COLS2, row_index_1, row_index_2;
	int iterations;
	double RMSIP, max_angle;
	Matrix Avg_RMSIP_Score_Matrix, RMSIP_Std_Dev_Matrix, Avg_CO_Score_Matrix, Avg_PA_Matrix, projections, projections_abs, projections_squared, CO_matrix_1_2, CO_matrix_2_1,
			input_matrix_1, input_matrix_2, PAs;
	ArrayList<Double> projections_squared_CO_k, projections_squared_RMSIP, cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products,
			vectorial_sum_of_angles, RMSIPs;
	final Descriptive_Stats ds;

	/* ***************************** CONSTRUCTOR ************************************************************* */

	/**
	 * Constructor to initiate the subspace analysis of two sets of eigenvectors. The sets of eigenvectors are equidimensional bases for subspaces of the eigenspaces of the
	 * trajectories.
	 * 
	 * @param basis1 The first set of eigenvectors.
	 * @param basis2 The second set of eigenvectors.
	 */
	public Subspace_Analysis(Matrix basis1, Matrix basis2)
	{
		input_matrix_1 = basis1;
		input_matrix_2 = basis2;

		ROWS1 = basis1.getRowDimension();
		COLS1 = basis1.getColumnDimension();
		ROWS2 = basis2.getRowDimension();
		COLS2 = basis2.getColumnDimension();

		this.ds = new Descriptive_Stats();

		if (ROWS1 != ROWS2)
			{
				System.err.println("ERROR: The subspaces do not come from the same vector space. Program will terminate.");
				System.err.println("ROWS1 = " + ROWS1);
				System.err.println("ROWS2 = " + ROWS2);
				// System.exit(0); program will crash and give stack trace.
			}
		if (COLS1 != COLS2)
			{
				System.err.println("ERROR: The subspaces do not have the same dimensions. Program will terminate");
				System.err.println("COLS1 = " + COLS1);
				System.err.println("COLS2 = " + COLS2);
				// System.exit(0); program will crash and give stack trace.
			}
		row_index_1 = 0;
		row_index_2 = ROWS1 - 1;
	}

	/* ******************************** METHODS **************************************************************** */

	/**
	 * Method to compute the subspace analysis: No iteration, all metrics.
	 */
	public synchronized void get_SSA()
	{
		CO_matrix_1_2 = new Matrix(COLS1, 1);
		CO_matrix_2_1 = new Matrix(COLS1, 1);

		projections_squared_CO_k = new ArrayList<>();
		projections_squared_RMSIP = new ArrayList<>();
		cumulative_overlaps_1_2 = new ArrayList<>();
		cumulative_overlaps_2_1 = new ArrayList<>();
		principle_angles_svd = new ArrayList<>();
		cosine_products = new ArrayList<>();
		vectorial_sum_of_angles = new ArrayList<>();

		projections = input_matrix_1.transpose().times(input_matrix_2);
		projections_squared = projections.arrayTimes(projections);

		get_CO_1_2();
		get_CO_2_1();
		get_RMSIP();

		Principal_Angles pa = new Principal_Angles(projections);
		principle_angles_svd = pa.get_Principle_Angles_Degrees();
		cosine_products = pa.get_Cosine_Products();
		max_angle = pa.get_Max_Angle();
		double sum_PA_squared = 0;
		double vector_sum_of_angles = 0;
		for (double d : principle_angles_svd)
			{
				sum_PA_squared += (d * d);
				vector_sum_of_angles = Math.sqrt(sum_PA_squared);
				vectorial_sum_of_angles.add(vector_sum_of_angles);
			}
	}

	/**
	 * Calculates the RMSIP.
	 */
	private synchronized void get_RMSIP()
	{
		int cols = projections_squared.getColumnDimension();

		double cumulative_overlap_RMSIP = 0;
		for (double d : projections_squared_RMSIP)
			{
				cumulative_overlap_RMSIP += d;
			}
		cumulative_overlap_RMSIP = (cumulative_overlap_RMSIP / cols);
		RMSIP = Math.sqrt(cumulative_overlap_RMSIP);
	}

	/**
	 * Calculates the CO_K for every vector in subspace 1 with subspace 2.
	 */
	private synchronized void get_CO_1_2()
	{
		int cols = projections_squared.getColumnDimension();

		for (int k = 0; k < cols; k++)
			{
				Matrix row = projections_squared.getMatrix(k, k, 0, cols - 1);
				double cumulative_overlap = 0;

				for (int j = 0; j < cols; j++)
					{
						double element = row.get(0, j);
						projections_squared_CO_k.add(element);
						projections_squared_RMSIP.add(element);
						cumulative_overlap += element;
					}
				double CO_k = Math.sqrt(cumulative_overlap);
				cumulative_overlaps_1_2.add(CO_k);
				CO_matrix_1_2.set(k, 0, CO_k);
				projections_squared_CO_k.clear();
			}
	}

	/**
	 * Calculates the CO_K for every vector in subspace 2 with subspace 1.
	 */
	private synchronized void get_CO_2_1()
	{
		int cols = projections_squared.getColumnDimension();

		for (int k = 0; k < COLS2; k++)
			{
				Matrix col = projections_squared.getMatrix(0, cols - 1, k, k);
				double cumulative_overlap = 0;

				for (int j = 0; j < cols; j++)
					{
						double element = col.get(j, 0);
						projections_squared_CO_k.add(element);
						cumulative_overlap += element;
					}
				double CO_k = Math.sqrt(cumulative_overlap);
				cumulative_overlaps_2_1.add(CO_k);
				CO_matrix_2_1.set(k, 0, CO_k);
				projections_squared_CO_k.clear();
			}
	}

	/**
	 * Method to compute the fast subspaces analysis: No iteration, only RMSIP and PA
	 */
	public synchronized void get_fast_SSA()
	{
		CO_matrix_1_2 = new Matrix(COLS1, 1);

		projections_squared_CO_k = new ArrayList<>();
		projections_squared_RMSIP = new ArrayList<>();
		cumulative_overlaps_1_2 = new ArrayList<>();
		principle_angles_svd = new ArrayList<>();

		projections = input_matrix_1.transpose().times(input_matrix_2);
		projections_squared = projections.arrayTimes(projections);

		get_CO_1_2();
		get_RMSIP();

		Principal_Angles pa = new Principal_Angles(projections);
		principle_angles_svd = pa.get_Principle_Angles_Degrees();
	}

	/**
	 * Method to compute the iterated fast subspace analysis. Iterates from dim 1 to the dim of the entered data. RMSIPs and PAs are calculated.
	 */
	public synchronized void get_fast_SSA_iterated()
	{
		RMSIPs = new ArrayList<>();
		PAs = new Matrix(COLS1, COLS1);
		CO_matrix_1_2 = new Matrix(COLS1, 1);

		Matrix orig_projections = input_matrix_1.transpose().times(input_matrix_2);
		Matrix orig_projections_squared = orig_projections.arrayTimes(orig_projections);

		for (int index = 0; index < COLS1; index++)
			{
				projections_squared_CO_k = new ArrayList<>();
				projections_squared_RMSIP = new ArrayList<>();
				cumulative_overlaps_1_2 = new ArrayList<>();
				principle_angles_svd = new ArrayList<>();

				projections = orig_projections.getMatrix(0, index, 0, index);
				projections_squared = orig_projections_squared.getMatrix(0, index, 0, index);

				get_CO_1_2();
				get_RMSIP();
				RMSIPs.add(RMSIP);

				Principal_Angles pa = new Principal_Angles(projections);
				principle_angles_svd = pa.get_Principle_Angles_Degrees();

				for (int i = 0; i <= index; i++)
					{
						double element = principle_angles_svd.get(i);
						PAs.set(index, i, element);
					}
			}
	}

	/**
	 * Method to compute an iterated subspace analysis of random bases. This can be a used to form a baseline for statistical significance of RMSIP and PA scores.
	 */
	public synchronized void get_random_SSA()
	{
		if (ROWS1 <= 1000) iterations = ROWS1;
		if (ROWS1 > 1000) iterations = 1000;

		final double scale = Math.pow(iterations, -1);

		Avg_RMSIP_Score_Matrix = new Matrix(1, COLS1);
		RMSIP_Std_Dev_Matrix = new Matrix(1, COLS1);
		Avg_CO_Score_Matrix = new Matrix(COLS1, COLS1);
		Avg_PA_Matrix = new Matrix(COLS1, COLS1);

		Matrix RMSIP_Score_Matrix = new Matrix(iterations, COLS1);

		for (int i = 0; i < iterations; i++)
			{
				Matrix CO_Score_Matrix = new Matrix(COLS1, COLS1);
				Matrix PA_Matrix = new Matrix(COLS1, COLS1);

				for (int SS_DIM = 1; SS_DIM < COLS1 + 1; SS_DIM++)
					{
						Matrix top_evects1_ss = input_matrix_1.getMatrix(row_index_1, row_index_2, 0, SS_DIM - 1);
						Matrix top_evects2_ss = input_matrix_2.getMatrix(row_index_1, row_index_2, 0, SS_DIM - 1);

						Matrix_Column_Permutation mcp = new Matrix_Column_Permutation(top_evects1_ss);
						Matrix rand1 = mcp.Get_Random_Permuted_Matrix();
						mcp = new Matrix_Column_Permutation(top_evects2_ss);
						Matrix rand2 = mcp.Get_Random_Permuted_Matrix();

						Subspace_Analysis rand_ssa = new Subspace_Analysis(rand1, rand2);
						rand_ssa.get_SSA();

						double RMSIP = rand_ssa.getRMSIP();
						RMSIP_Score_Matrix.set(i, (SS_DIM - 1), RMSIP);
						Matrix co = rand_ssa.getCO_matrix_1_2();
						principle_angles_svd = rand_ssa.getPrinciple_angles_svd();

						Matrix pa_col = new Matrix(SS_DIM, 1);
						for (int pa = 0; pa < SS_DIM; pa++)
							{
								pa_col.set(pa, 0, principle_angles_svd.get(pa));
							}
						CO_Score_Matrix.setMatrix(0, (SS_DIM - 1), (SS_DIM - 1), (SS_DIM - 1), co);
						PA_Matrix.setMatrix(0, (SS_DIM - 1), (SS_DIM - 1), (SS_DIM - 1), pa_col);
					}
				Avg_CO_Score_Matrix.plusEquals(CO_Score_Matrix);
				Avg_PA_Matrix.plusEquals(PA_Matrix);

			}

		Avg_CO_Score_Matrix = Avg_CO_Score_Matrix.times(scale);
		Avg_PA_Matrix = Avg_PA_Matrix.times(scale);

		for (int i = 0; i < COLS1; i++)
			{

				Matrix SS_Stats = RMSIP_Score_Matrix.getMatrix(0, iterations - 1, i, i);
				double[] stat_array = SS_Stats.getColumnPackedCopy();
				double avg = ds.get_mean(stat_array);
				double std_dev = ds.get_standard_deviation(stat_array, avg);
				Avg_RMSIP_Score_Matrix.set(0, i, avg);
				RMSIP_Std_Dev_Matrix.set(0, i, std_dev);
			}
	}

	/* ******************************** GETTERS **************************************************************** */

	/**
	 * @return the ROWS1
	 */
	public int getROWS1()
	{
		return ROWS1;
	}

	/**
	 * @return the COLS1
	 */
	public int getCOLS1()
	{
		return COLS1;
	}

	/**
	 * @return the ROWS2
	 */
	public int getROWS2()
	{
		return ROWS2;
	}

	/**
	 * @return the COLS2
	 */
	public int getCOLS2()
	{
		return COLS2;
	}

	/**
	 * @return the row_index_1
	 */
	public int getRow_index_1()
	{
		return row_index_1;
	}

	/**
	 * @return the row_index_2
	 */
	public int getRow_index_2()
	{
		return row_index_2;
	}

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}

	/**
	 * @return the RMSIP
	 */
	public double getRMSIP()
	{
		return RMSIP;
	}

	/**
	 * @return the max_angle
	 */
	public double getMax_angle()
	{
		return max_angle;
	}


	/**
	 * @return the avg_RMSIP_Score_Matrix
	 */
	public Matrix getAvg_RMSIP_Score_Matrix()
	{
		return Avg_RMSIP_Score_Matrix;
	}

	/**
	 * @return the RMSIP_Std_Dev_Matrix
	 */
	public Matrix getRMSIP_Std_Dev_Matrix()
	{
		return RMSIP_Std_Dev_Matrix;
	}

	/**
	 * @return the avg_CO_Score_Matrix
	 */
	public Matrix getAvg_CO_Score_Matrix()
	{
		return Avg_CO_Score_Matrix;
	}

	/**
	 * @return the avg_PA_Matrix
	 */
	public Matrix getAvg_PA_Matrix()
	{
		return Avg_PA_Matrix;
	}

	/**
	 * @return the projections
	 */
	public Matrix getProjections()
	{
		return projections;
	}

	/**
	 * @return the projections_abs
	 */
	public Matrix getProjections_abs()
	{
		return projections_abs;
	}

	/**
	 * @return the projections_squared
	 */
	public Matrix getProjections_squared()
	{
		return projections_squared;
	}

	/**
	 * @return the CO_matrix_1_2
	 */
	public Matrix getCO_matrix_1_2()
	{
		return CO_matrix_1_2;
	}

	/**
	 * @return the CO_matrix_2_1
	 */
	public Matrix getCO_matrix_2_1()
	{
		return CO_matrix_2_1;
	}

	/**
	 * @return the input_matrix_1
	 */
	public Matrix getInput_matrix_1()
	{
		return input_matrix_1;
	}

	/**
	 * @return the input_matrix_2
	 */
	public Matrix getInput_matrix_2()
	{
		return input_matrix_2;
	}

	/**
	 * @return the PAs
	 */
	public Matrix getPAs()
	{
		return PAs;
	}

	/**
	 * @return the projections_squared_CO_k
	 */
	public ArrayList<Double> getProjections_squared_CO_k()
	{
		return projections_squared_CO_k;
	}

	/**
	 * @return the projections_squared_RMSIP
	 */
	public ArrayList<Double> getProjections_squared_RMSIP()
	{
		return projections_squared_RMSIP;
	}

	/**
	 * @return the cumulative_overlaps_1_2
	 */
	public ArrayList<Double> getCumulative_overlaps_1_2()
	{
		return cumulative_overlaps_1_2;
	}

	/**
	 * @return the cumulative_overlaps_2_1
	 */
	public ArrayList<Double> getCumulative_overlaps_2_1()
	{
		return cumulative_overlaps_2_1;
	}

	/**
	 * @return the principle_angles_svd
	 */
	public ArrayList<Double> getPrinciple_angles_svd()
	{
		return principle_angles_svd;
	}

	/**
	 * @return the cosine_products
	 */
	public ArrayList<Double> getCosine_products()
	{
		return cosine_products;
	}

	/**
	 * @return the vectorial_sum_of_angles
	 */
	public ArrayList<Double> getVectorial_sum_of_angles()
	{
		return vectorial_sum_of_angles;
	}

	/**
	 * @return the rMSIPs
	 */
	public ArrayList<Double> getRMSIPs()
	{
		return RMSIPs;
	}
}
