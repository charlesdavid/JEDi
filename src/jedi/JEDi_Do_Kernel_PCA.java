package jedi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import Jama.Matrix;
import support.Descriptive_Stats;
import support.Mutual_Information;
import support.PCA;
import support.Row_Center_Data;
import supportIO.Input_Parameters;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportPlot.PC_Plot;

public class JEDi_Do_Kernel_PCA
{
	boolean exist, success, verbose, noSolution;

	boolean do_linear_kernel = false;
	boolean do_Degree_2_Poly = true;
	boolean do_Degree_3_Poly = false;
	boolean do_Degree_4_Poly = false;
	boolean do_XY_Poly = false;
	boolean do_Euclidean = false;
	boolean do_Mahalanobis = false;
	boolean do_Gaussian = true;
	boolean do_Sigmoid = true;
	boolean do_Log = false;
	boolean do_Circular = false;
	boolean do_Cauchy = false;
	boolean do_MI = false;
	boolean do_MI_KDE = false;

	double sigma;
	double slope;
	double MI_Kernel_Resolution;
	double MI_Kernel_Cell;
	double MI_Kernel_Margin;
	double FLOOR;
	double shrinkage;

	int NUMBER_PCs;
	int MAX_FRAMES;
	int multiplier;

	int ROWS, COLS, KERNEL_FRAMES, DOWNSAMPLE, num_top_evals, num_non_zero_evals;
	double[] evals, means;

	String path, directory, out_dir, description, name;
	Matrix data, data_row_centered, data_reduced, data_reduced_transposed_centered, data_reduced_transposed, K, S, S_inv;
	ArrayList<Double> non_zero_evals, top_evals;

	Jama.EigenvalueDecomposition evd;
	Mutual_Information MI;

	final PCA pca;
	final Row_Center_Data rcd;
	final Descriptive_Stats ds;

	final double SQ_RT_2 = Math.sqrt(2);
	final double SQ_RT_3 = Math.sqrt(3);
	final double SQ_RT_6 = Math.sqrt(6);

	PrintWriter writer;
	File input;
	BufferedReader reader;

	NumberFormat nf3;
	RoundingMode rm;

	/* ************************************** CONSTRUCTORS ***************************************************************************** */

	public JEDi_Do_Kernel_PCA(String output_dir, Matrix projections, int num_kpca_modes)
	{
		this.out_dir = output_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();

		this.data = projections;
		this.num_top_evals = num_kpca_modes;

		verbose = true;

		data_reduced = data.getMatrix(0, data.getRowDimension() - 1, 0, NUMBER_PCs);
		data_reduced_transposed = data_reduced.transpose();

		rcd = new Row_Center_Data(data_reduced_transposed);
		data_reduced_transposed_centered = rcd.get_row_centered_data();

		pca = new PCA(data_reduced_transposed_centered);
		ds = new Descriptive_Stats();

		Matrix sigmaMatrix = rcd.get_variable_variances();
		double[] sigmaArray = sigmaMatrix.getRowPackedCopy();
		double varianceMean = ds.get_mean(sigmaArray);
		sigma = varianceMean;

		this.ROWS = NUMBER_PCs;
		this.KERNEL_FRAMES = data_reduced.getColumnDimension();
		this.num_top_evals = ROWS;

		rm = RoundingMode.HALF_UP;
		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);
	}

	public JEDi_Do_Kernel_PCA(Matrix normedPCs, String out) // Called by JEDi_Driver
	{
		this.data = normedPCs; // This should be the normed projections matrix from a PCA analysis (in columns)
		this.out_dir = out;
		this.verbose = Input_Parameters.verbose;

		create_Directory(out_dir);

		this.sigma = Input_Parameters.KPCA_SIGMA;
		this.slope = Input_Parameters.KPCA_SLOPE;
		this.MI_Kernel_Resolution = Input_Parameters.KDE_RESOLUTION;
		this.MI_Kernel_Cell = Input_Parameters.KDE_CELL;
		this.MI_Kernel_Margin = Input_Parameters.KDE_MARGIN;
		this.FLOOR = Input_Parameters.FLOOR;
		this.shrinkage = Input_Parameters.KERNEL_SHRINKAGE;
		this.NUMBER_PCs = Input_Parameters.NUMBER_PCs_INPUT;
		this.MAX_FRAMES = Input_Parameters.MAX_KERNEL_FRAMES;
		this.multiplier = Input_Parameters.MULTIPLIER;

		this.do_linear_kernel = Input_Parameters.Linear;
		this.do_Degree_2_Poly = Input_Parameters.Degree_2_Poly;
		this.do_Degree_3_Poly = Input_Parameters.Degree_3_Poly;
		this.do_Degree_4_Poly = Input_Parameters.Degree_4_Poly;
		this.do_XY_Poly = Input_Parameters.XY_Poly;
		this.do_Euclidean = Input_Parameters.Euclidean;
		this.do_Mahalanobis = Input_Parameters.Mahalanobis;
		this.do_Gaussian = Input_Parameters.Gaussian;
		this.do_Sigmoid = Input_Parameters.Sigmoid;
		this.do_Log = Input_Parameters.Log;
		this.do_Circular = Input_Parameters.Circular;
		this.do_Cauchy = Input_Parameters.Cauchy;
		this.do_MI = Input_Parameters.MI;
		this.do_MI_KDE = Input_Parameters.MI_KDE;

		rm = RoundingMode.HALF_UP;
		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);

		data_reduced = data.getMatrix(0, data.getRowDimension() - 1, 0, NUMBER_PCs - 1); // This limits the kernel PCA to the top PCs
		data_reduced_transposed = data_reduced.transpose(); // Transposing sets the rows to be the generalized variables

		rcd = new Row_Center_Data(data_reduced_transposed);
		data_reduced_transposed_centered = rcd.get_row_centered_data();

		this.ROWS = NUMBER_PCs;
		this.COLS = data_reduced_transposed_centered.getColumnDimension();

		if (COLS > MAX_FRAMES)
			{
				int ds = (COLS / MAX_FRAMES);
				double rm = (COLS % MAX_FRAMES);
				if (rm > 0) ds += 1;
				set_DOWNSAMPLE(ds);
				KERNEL_FRAMES = do_Down_Sample();
			}
		else
			this.KERNEL_FRAMES = COLS;

		this.num_top_evals = ROWS;

		pca = new PCA(data_reduced_transposed_centered);
		S = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		S_inv = S.inverse();
		K = pca.get_K_matrix(); // For Linear Kernel

		// MI = new Mutual_Information(out_dir, data_reduced_transposed_centered);
		ds = new Descriptive_Stats();

		Matrix u = rcd.get_variable_means();
		means = u.getRowPackedCopy();

		Matrix sigmaMatrix = rcd.get_variable_variances();
		double[] sigmaArray = sigmaMatrix.getRowPackedCopy();
		double varianceMean = ds.get_mean(sigmaArray);

		if (Input_Parameters.KPCA_SIGMA == 0) sigma = (varianceMean * multiplier);
		if (Input_Parameters.KPCA_SLOPE == 0) slope = ((1d / NUMBER_PCs) * multiplier);
	}

	public JEDi_Do_Kernel_PCA(Matrix PCs, String dir, int num_KPCs, String desc, int max_frames) // Called by KPCA_Driver
	{
		this.directory = dir;
		this.description = desc;
		this.data = PCs; // This should be the normed projections matrix from a PCA analysis (in columns)
		this.NUMBER_PCs = num_KPCs;
		this.MAX_FRAMES = max_frames;
		this.verbose = true;

		rm = RoundingMode.HALF_UP;
		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);

		this.out_dir = directory + "JEDi_KPCA_RESULTS_" + description + File.separatorChar;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();

		data_reduced = data.getMatrix(0, data.getRowDimension() - 1, 0, NUMBER_PCs - 1);
		data_reduced_transposed = data_reduced.transpose(); // Transposing sets the rows to be the generalized variables

		rcd = new Row_Center_Data(data_reduced_transposed);
		data_reduced_transposed_centered = rcd.get_row_centered_data();

		this.ROWS = NUMBER_PCs;
		this.COLS = data_reduced_transposed_centered.getColumnDimension();

		if (COLS > MAX_FRAMES)
			{
				int ds = (COLS / MAX_FRAMES);
				double rm = (COLS % MAX_FRAMES);
				if (rm > 0) ds += 1;
				set_DOWNSAMPLE(ds);
				KERNEL_FRAMES = do_Down_Sample();
			}
		else
			this.KERNEL_FRAMES = COLS;

		this.num_top_evals = ROWS;

		pca = new PCA(data_reduced_transposed_centered);
		S = pca.get_Covariance_Estimator_Using_Optimized_Shrinkage();
		S_inv = S.inverse();
		K = pca.get_K_matrix(); // For Linear Kernel

		MI = new Mutual_Information(out_dir, data_reduced_transposed_centered);
		ds = new Descriptive_Stats();

		Matrix u = rcd.get_variable_means();
		means = u.getRowPackedCopy();
	}

	/* ****************************************** METHODS ****************************************************************************** */
	public Matrix get_data_matrix(String path)
	{
		try
			{
				input = new File(path);
				reader = new BufferedReader(new FileReader(input));
				data = Matrix.read(reader);
				reader.close();

				data_reduced = data.getMatrix(0, data.getRowDimension() - 1, 0, NUMBER_PCs - 1);
				data_reduced_transposed = data_reduced.transpose();

				Row_Center_Data rcd = new Row_Center_Data(data_reduced_transposed);
				data_reduced_transposed_centered = rcd.get_row_centered_data();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return data;
	}

	private int do_Down_Sample()
	{
		int size = (COLS / DOWNSAMPLE);
		int[] column_indices = new int[size];
		for (int i = 0; i < size; i++)
			{
				column_indices[i] = i * DOWNSAMPLE;
			}
		Matrix Reduced_Matrix = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, column_indices);
		data_reduced_transposed_centered = Reduced_Matrix;
		return data_reduced_transposed_centered.getColumnDimension();
	}

	public Matrix get_Gaussian_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);
		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);
						Matrix diff = data_point_X.minus(data_point_Y);
						Matrix diff_sq = diff.transpose().times(diff);
						double dp = diff_sq.get(0, 0);
						double arg = (dp / (2 * sigma * sigma));
						double value = Math.exp((-1) * arg);
						kernel.set(col_index_1, col_index_2, value);
						kernel.set(col_index_2, col_index_1, value);
					}
			}
		return kernel;
	}

	public Matrix get_Log_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);
		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);
						Matrix diff = data_point_X.minus(data_point_Y);
						Matrix diff_sq = diff.transpose().times(diff);
						double dp = diff_sq.get(0, 0);
						double arg = (dp + 1d);
						double value = ((-1) * Math.log(arg));
						kernel.set(col_index_1, col_index_2, value);
						kernel.set(col_index_2, col_index_1, value);
					}
			}
		return kernel;
	}

	public Matrix get_Circular_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);
		final double C = (2 / Math.PI);
		double value = 0;
		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);
						Matrix diff = data_point_X.minus(data_point_Y);
						double norm = diff.normF();
						if (norm < sigma)
							{
								double norm_squared = (norm * norm);
								double arg1 = Math.acos((-1d) * (norm / sigma));
								double arg2 = (norm / sigma) * Math.sqrt(1d - (norm_squared / (sigma * sigma)));
								value = (C * (arg1 - arg2));
							}
						kernel.set(col_index_1, col_index_2, value);
						kernel.set(col_index_2, col_index_1, value);
						value = 0;
					}
			}
		return kernel;
	}

	public Matrix get_Cauchy_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);
		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);
						Matrix diff = data_point_X.minus(data_point_Y);
						double norm = diff.normF();
						double norm_squared = (norm * norm);
						double arg = (1d + (norm_squared / (sigma * sigma)));
						double value = (1d / arg);

						kernel.set(col_index_1, col_index_2, value);
						kernel.set(col_index_2, col_index_1, value);
					}
			}
		return kernel;
	}

	public Matrix get_Sigmoid_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);
		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);
						Matrix prod = data_point_X.transpose().times(data_point_Y);
						double dp = prod.get(0, 0);
						double value = Math.tanh(slope * dp);
						kernel.set(col_index_1, col_index_2, value);
						kernel.set(col_index_2, col_index_1, value);
					}
			}
		return kernel;
	}

	public Matrix get_Euclidean_Distance_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						Matrix diff = data_point_X.minus(data_point_Y);
						double dist = diff.normF();
						double similarity = (1.0d - dist);

						kernel.set(col_index_1, col_index_2, similarity);
						kernel.set(col_index_2, col_index_1, similarity);
					}
			}
		return kernel;
	}

	public Matrix get_Mahalanobis_Distance_Kernel_Matrix()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						Matrix diff = data_point_X.minus(data_point_Y);
						Matrix Dxy = diff.transpose().times(S_inv).times(diff);

						double dp = Dxy.get(0, 0);
						double dist = Math.sqrt(dp);
						double similarity = (1.0d - dist);

						kernel.set(col_index_1, col_index_2, similarity);
						kernel.set(col_index_2, col_index_1, similarity);
					}
			}
		return kernel;
	}

	public Matrix get_Polynomial_Kernel_Matrix_Degree_2()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						double dp = 0;

						for (int i = 0; i < ROWS; i++)
							{
								double X = data_point_X.get(i, 0);
								double Y = data_point_Y.get(i, 0);

								dp += (X * X) + (SQ_RT_2 * X * Y) + (Y * Y);
							}
						kernel.set(col_index_1, col_index_2, dp);
						kernel.set(col_index_2, col_index_1, dp);
					}
			}
		return kernel;
	}

	public Matrix get_Polynomial_Kernel_Matrix_Degree_3()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						double dp = 0;

						for (int i = 0; i < ROWS; i++)
							{
								double X = data_point_X.get(i, 0);
								double Y = data_point_Y.get(i, 0);

								dp += (X * X * X) + (SQ_RT_3 * X * X * Y) + (SQ_RT_3 * X * Y * Y) + (Y * Y * Y);
							}
						kernel.set(col_index_1, col_index_2, dp);
						kernel.set(col_index_2, col_index_1, dp);
					}
			}
		return kernel;
	}

	public Matrix get_Polynomial_Kernel_Matrix_Degree_4()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						double dp = 0;

						for (int i = 0; i < ROWS; i++)
							{
								double X = data_point_X.get(i, 0);
								double Y = data_point_Y.get(i, 0);

								dp += (X * X * X * X) + (2 * X * X * X * Y) + (SQ_RT_6 * X * X * Y * Y) + (2 * X * Y * Y * Y) + (Y * Y * Y * Y);
							}
						kernel.set(col_index_1, col_index_2, dp);
						kernel.set(col_index_2, col_index_1, dp);
					}
			}
		return kernel;
	}

	public Matrix get_Polynomial_Kernel_Matrix_XY()
	{
		Matrix kernel = new Matrix(KERNEL_FRAMES, KERNEL_FRAMES);

		for (int col_index_1 = 0; col_index_1 < KERNEL_FRAMES; col_index_1++)
			{
				for (int col_index_2 = col_index_1; col_index_2 < KERNEL_FRAMES; col_index_2++)
					{
						Matrix data_point_X = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_1, col_index_1);
						Matrix data_point_Y = data_reduced_transposed_centered.getMatrix(0, ROWS - 1, col_index_2, col_index_2);

						double dp = 0;

						for (int i = 0; i < ROWS; i++)
							{
								double X = data_point_X.get(i, 0);
								double Y = data_point_Y.get(i, 0);

								for (int j = 0; j < ROWS; j++)
									{
										double value_X = data_point_X.get(j, 0);
										double value_Y = data_point_Y.get(j, 0);
										double X_comp = X * value_X;
										double Y_comp = Y * value_Y;
										dp += X_comp * Y_comp;
									}
							}
						kernel.set(col_index_1, col_index_2, dp);
						kernel.set(col_index_2, col_index_1, dp);
					}
			}
		return kernel;
	}

	public Matrix get_MI_Kernel_Matrix_Version1()
	{
		Matrix kernel = MI.get_MI_Kernel_Version1();
		return kernel;
	}

	public Matrix get_MI_Kernel_Matrix_Version2()
	{
		Matrix kernel = MI.get_MI_Kernel_Version2();
		return kernel;
	}

	public Matrix get_MI_Kernel_Matrix_Version3()
	{
		Matrix kernel = MI.get_MI_Kernel_Version3();
		return kernel;
	}

	public Matrix get_MI_Kernel_Matrix_KDE()
	{
		MI.set_kernel_resolution(MI_Kernel_Resolution);
		Matrix kernel = MI.get_MI_Kernel_2D_KDE();
		return kernel;
	}

	public Matrix get_MI_Kernel_Matrix_KDE_Version_2()
	{
		return MI.get_MI_Kernel_2D_KDE_Version_2();
	}

	// ******************************************************************************************************** //

	public Matrix get_Normalized_Kernel_Eigenvectors(Matrix top_evects, ArrayList<Double> top_evals)
	{
		int rows = top_evects.getRowDimension();
		int cols = top_evects.getColumnDimension();
		Matrix scaled_evects = new Matrix(rows, cols);
		for (int i = 0; i < top_evals.size(); i++)
			{
				double val = top_evals.get(i);
				if (val < FLOOR)
					{
						val = FLOOR;
					}
				Matrix col = top_evects.getMatrix(0, rows - 1, i, i);
				double scale = (1d / Math.sqrt(val));
				Matrix normed_col = col.times(scale);
				scaled_evects.setMatrix(0, rows - 1, i, i, normed_col);
			}
		return scaled_evects;
	}

	public Matrix get_Normalized_Kernel_Eigenvectors(Matrix evects, double[] evals)
	{
		num_non_zero_evals = 0;
		non_zero_evals = new ArrayList<Double>();
		for (double d : evals)
			{
				if (d > FLOOR)
					{
						num_non_zero_evals++;
						non_zero_evals.add(d);
					}
			}
		Matrix evects_reduced = evects.getMatrix(0, KERNEL_FRAMES - 1, (KERNEL_FRAMES - num_non_zero_evals), KERNEL_FRAMES - 1);
		Matrix scaled_evects = new Matrix(KERNEL_FRAMES, num_non_zero_evals);
		for (int i = 0; i < num_non_zero_evals; i++)
			{
				double val = non_zero_evals.get(i);
				Matrix col = evects_reduced.getMatrix(0, KERNEL_FRAMES - 1, i, i);
				double scale = (1d / Math.sqrt(val));
				Matrix normed_col = col.times(scale);
				scaled_evects.setMatrix(0, KERNEL_FRAMES - 1, i, i, normed_col);
			}
		return scaled_evects;
	}

	public Matrix get_Top_Normalized_Kernel_Eigenvectors(Matrix evects, double[] evals)
	{
		int num_of_evals = evals.length;
		ArrayList<Double> top_evals = new ArrayList<Double>();
		for (int i = 0; i < num_top_evals; i++)
			{
				double d = evals[num_of_evals - num_top_evals + i];
				top_evals.add(d);
			}
		Matrix top_eigenvectors = evects.getMatrix(0, KERNEL_FRAMES - 1, (KERNEL_FRAMES - num_top_evals), KERNEL_FRAMES - 1);
		Matrix scaled_top_evects = new Matrix(KERNEL_FRAMES, num_top_evals);
		for (int i = 0; i < num_top_evals; i++)
			{
				int col_index = num_top_evals - i - 1;
				double val = top_evals.get(i);
				if (val < FLOOR)
					{
						val = FLOOR;
					}
				Matrix col = top_eigenvectors.getMatrix(0, KERNEL_FRAMES - 1, col_index, col_index);
				double scale = (1d / Math.sqrt(val));
				Matrix normed_col = col.times(scale);
				scaled_top_evects.setMatrix(0, KERNEL_FRAMES - 1, i, i, normed_col);
			}
		return scaled_top_evects;
	}

	public ArrayList<Double> get_Top_Eigenvalues(double[] evals, String type)
	{
		ArrayList<Double> eigenvalues = new ArrayList<Double>();
		ArrayList<Double> Top_eigenvalues = new ArrayList<Double>();
		for (double k : evals)
			{
				eigenvalues.add(k);
			}
		Collections.sort(eigenvalues, Collections.reverseOrder());

		for (int i = 0; i < num_top_evals; i++)
			{
				Top_eigenvalues.add(eigenvalues.get(i));
			}
		if (verbose) List_IO.write_Double_List(eigenvalues, out_dir + type + "_Eigenvalues.txt", 12);
		// List_IO.write_Double_List(Top_eigenvalues, out_dir + type + "_Top_" + num_top_evals + "_eigenvalues.txt", 12);
		return Top_eigenvalues;
	}

	public Matrix get_Top_Eigenvectors(Matrix evects, String type)
	{
		int rows = evects.getRowDimension();
		int cols = evects.getColumnDimension();
		Matrix top_eigenvectors = new Matrix(rows, num_top_evals);

		for (int i = 0; i < num_top_evals; i++)
			{
				Matrix eigenvector = evects.getMatrix(0, rows - 1, cols - i - 1, cols - i - 1);
				top_eigenvectors.setMatrix(0, rows - 1, i, i, eigenvector);
			}
		if (verbose) Matrix_IO.write_BZ2_Matrix(top_eigenvectors, out_dir + type + "_Top_Eigenvectors.txt.bz2", 20, 12);

		return top_eigenvectors;
	}

	public Matrix get_kPCs_using_Data(Matrix k, Matrix scaled_evects)
	{
		int rows = scaled_evects.getRowDimension();
		int num_of_evects = scaled_evects.getColumnDimension();
		int num_of_frames = k.getColumnDimension();
		Matrix kPCs = new Matrix(rows, num_of_evects);

		for (int i = 0; i < num_of_evects; i++)
			{
				Matrix alpha = scaled_evects.getMatrix(0, rows - 1, i, i);
				for (int j = 0; j < num_of_frames; j++)
					{
						Matrix col = k.getMatrix(0, rows - 1, j, j);
						double dp = col.transpose().times(alpha).get(0, 0);
						kPCs.set(j, i, dp);
					}
			}
		return kPCs;
	}

	public Matrix get_Top_kPCs_using_Data(Matrix k, Matrix scaled_top_evects)
	{
		Matrix top_kPCs = new Matrix(KERNEL_FRAMES, num_top_evals);
		for (int i = 0; i < num_top_evals; i++)
			{
				Matrix alpha = scaled_top_evects.getMatrix(0, KERNEL_FRAMES - 1, i, i);
				for (int j = 0; j < KERNEL_FRAMES; j++)
					{
						Matrix col = k.getMatrix(0, KERNEL_FRAMES - 1, j, j);
						double dp = col.transpose().times(alpha).get(0, 0);
						top_kPCs.set(j, i, dp);
					}
			}
		return top_kPCs;
	}

	/* *********************************************************************************************************** */
	public void do_MI_1()
	{
		Matrix kernel = get_MI_Kernel_Matrix_Version1();
		// System.out.println("MI 1 kernel matrix");
		// kernel.print(9, 6);
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "MI_Kernel_Version_1_Non_KDE");
		Top_Evals = get_Top_Eigenvalues(evals, "MI_Kernel_Version_1");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_Top_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "MI_Kernel_PCs_Version_1_Non_KDE.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "MI_Kernel_PCs_Version_1_Non_KDE", kpcs);
	}

	public void do_MI_2()
	{
		Matrix kernel = get_MI_Kernel_Matrix_Version2();
		// System.out.println("MI 2 kernel matrix");
		// kernel.print(9, 6);
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "MI_Kernel_Version_2_Non_KDE");
		Top_Evals = get_Top_Eigenvalues(evals, "MI_Kernel_Version_2_Non_KDE");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_Top_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "MI_Kernel_PCs_Version_2_Non_KDE.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "MI_Kernel_PCs_Version_2_Non_KDE", kpcs);
	}

	public void do_MI_3()
	{
		Matrix kernel = get_MI_Kernel_Matrix_Version3();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "MI_Kernel_Version_3_Non_KDE");
		Top_Evals = get_Top_Eigenvalues(evals, "MI_Kernel_Version_3_Non_KDE");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_Top_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "MI_Kernel_PCs_Version_3_Non_KDE.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "MI_Kernel_PCs_Version_3_Non_KDE_No_Marginal", kpcs);
	}

	public void do_MI_KDE()
	{
		if (!noSolution)
			{
				Matrix kernel = get_MI_Kernel_Matrix_KDE();
				// System.out.println("MI KDE kernel matrix");
				// kernel.print(9, 6);
				Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
				evd = new Jama.EigenvalueDecomposition(K);
				evals = evd.getRealEigenvalues();
				Matrix evects = evd.getV();
				ArrayList<Double> Top_Evals = new ArrayList<Double>();
				Matrix Top_Evects = get_Top_Eigenvectors(evects, "MI_Kernel_KDE");
				Top_Evals = get_Top_Eigenvalues(evals, "MI_Kernel_KDE");
				Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
				Matrix kpcs = get_Top_kPCs_using_Data(kernel, alpha);
				double[] bandwiths = MI.getBandwidths();
				double cellsize = MI.getCellSize();

				String name = "MI_Kernel_PCs_KDE_BW_" + nf3.format(bandwiths[0]) + "_" + nf3.format(bandwiths[3]) + "_CS_" + nf3.format(cellsize);
				if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + name + ".txt.bz2", 20, 12);
				PC_Plot.createChart4Series(out_dir, name, kpcs);
			}
	}

	public void do_MI_KDE_2()
	{
		Matrix kernel = get_MI_Kernel_Matrix_KDE_Version_2();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);

		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();

		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "MI_Kernel_KDE_Version_2");
		Top_Evals = get_Top_Eigenvalues(evals, "MI_Kernel_KDE_Version_2");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_Top_kPCs_using_Data(kernel, alpha);

		double[] bandwiths = MI.getBandwidths();
		double cellsize = MI.getCellSize();

		String name = "MI_Kernel_PCs_KDE_Version_2_BW_" + nf3.format(bandwiths[0]) + "_" + nf3.format(bandwiths[3]) + "_CS_" + nf3.format(cellsize);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + name + ".txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, name, kpcs);
	}

	public void do_Linear_Kernel()
	{
		Matrix kernel = pca.Shrink_COV_Matrix(K, shrinkage);
		evd = new Jama.EigenvalueDecomposition(kernel);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Linear_Kernel");
		Top_Evals = get_Top_Eigenvalues(evals, "Linear_Kernel");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Linear_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Linear_Kernel_PCs", kpcs);
	}

	public void do_Degree_2_Poly()
	{
		Matrix kernel = get_Polynomial_Kernel_Matrix_Degree_2();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Poly2");
		Top_Evals = get_Top_Eigenvalues(evals, "Poly2");

		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Polynomial_Degree_2_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Polynomial_Degree_2_Kernel_PCs", kpcs);
	}

	public void do_Degree_3_Poly()
	{
		Matrix kernel = get_Polynomial_Kernel_Matrix_Degree_3();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Poly3");
		Top_Evals = get_Top_Eigenvalues(evals, "Poly3");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Polynomial_Degree_3_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Polynomial_Degree_3_Kernel_PCs", kpcs);
	}

	public void do_Degree_4_Poly()
	{
		Matrix kernel = get_Polynomial_Kernel_Matrix_Degree_4();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Poly4");
		Top_Evals = get_Top_Eigenvalues(evals, "Poly4");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Polynomial_Degree_4_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Polynomial_Degree_4_Kernel_PCs", kpcs);
	}

	public void do_Poly_XY()
	{
		Matrix kernel = get_Polynomial_Kernel_Matrix_XY();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);

		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "PolyXY");
		Top_Evals = get_Top_Eigenvalues(evals, "PolyXY");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Polynomial_XY_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Polynomial_XY_Kernel_PCs", kpcs);
	}

	public void do_Euclidean()
	{
		Matrix kernel = get_Euclidean_Distance_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Euclidean");
		Top_Evals = get_Top_Eigenvalues(evals, "Euclidean");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Euclidean_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Euclidean_Kernel_PCs", kpcs);
	}

	public void do_Mahalanobis()
	{
		Matrix kernel = get_Mahalanobis_Distance_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();

		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Mahalanobis");
		Top_Evals = get_Top_Eigenvalues(evals, "Mahalanobis");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Mahalanobis_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Mahalanobis_Kernel_PCs", kpcs);
	}

	public void do_Gaussian()
	{
		Matrix kernel = get_Gaussian_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Gaussian");
		Top_Evals = get_Top_Eigenvalues(evals, "Gaussian");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Gaussian_Kernel_PCs_Sigma_" + nf3.format(sigma) + ".txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Gaussian_Kernel_PCs_Sigma_" + nf3.format(sigma), kpcs);
	}

	public void do_Sigmoid()
	{
		Matrix kernel = get_Sigmoid_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Sigmoid");
		Top_Evals = get_Top_Eigenvalues(evals, "Sigmoid");

		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Sigmoid_Kernel_PCs_Slope_" + nf3.format(slope) + ".txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Sigmoid_Kernel_PCs_Slope_" + nf3.format(slope), kpcs);
	}

	public void do_Log()
	{
		Matrix kernel = get_Log_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Log");
		Top_Evals = get_Top_Eigenvalues(evals, "Log");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Log_Kernel_PCs.txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Log_Kernel_PCs", kpcs);
	}

	public void do_Circular()
	{
		Matrix kernel = get_Circular_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();

		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Circular");
		Top_Evals = get_Top_Eigenvalues(evals, "Circular");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Circular_Kernel_PCs_Sigma_" + nf3.format(sigma) + ".txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Circular_Kernel_PCs_Sigma_" + nf3.format(sigma), kpcs);
	}

	public void do_Cauchy()
	{
		Matrix kernel = get_Cauchy_Kernel_Matrix();
		Matrix K = pca.Shrink_COV_Matrix(kernel, shrinkage);
		evd = new Jama.EigenvalueDecomposition(K);
		evals = evd.getRealEigenvalues();
		Matrix evects = evd.getV();
		ArrayList<Double> Top_Evals = new ArrayList<Double>();
		Matrix Top_Evects = get_Top_Eigenvectors(evects, "Cauchy");
		Top_Evals = get_Top_Eigenvalues(evals, "Cauchy");
		Matrix alpha = get_Normalized_Kernel_Eigenvectors(Top_Evects, Top_Evals);
		Matrix kpcs = get_kPCs_using_Data(kernel, alpha);

		if (verbose) Matrix_IO.write_BZ2_Matrix(kpcs, out_dir + "Cauchy_Kernel_PCs_Sigma_" + nf3.format(sigma) + ".txt.bz2", 20, 12);
		PC_Plot.createChart4Series(out_dir, "Cauchy_Kernel_PCs_Sigma_" + nf3.format(sigma), kpcs);
	}

	/* *************************** SETTERS ***************************************************** */

	public void set_DOWNSAMPLE(int downsample)
	{
		this.DOWNSAMPLE = downsample;
	}

	public void set_Sigma(double s)
	{
		this.sigma = s;
	}

	public void set_Slope(double m)
	{
		slope = m;
	}

	public void set_MI_kernelResolution(double resolution)
	{
		this.MI_Kernel_Resolution = resolution;
	}

	public void set_Num_top_evals(int num_evals)
	{
		this.num_top_evals = num_evals;
	}

	public void set_Out_Dir(String dir)
	{
		this.out_dir = dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}

	public void setMI_Kernel_Cell(double cell)
	{
		this.MI_Kernel_Cell = cell;
	}

	public void setShrinkage(double shrinkage)
	{
		this.shrinkage = shrinkage;
	}

	public void setMAX_FRAMES(int max)
	{
		this.MAX_FRAMES = max;
	}

	public void setMultiplier(int multiplier)
	{
		this.multiplier = multiplier;
	}

	public void setDOWNSAMPLE(int ds)
	{
		DOWNSAMPLE = ds;
	}

	public void setFLOOR(double floor)
	{
		this.FLOOR = floor;
	}

	/* ********************************************************************************************************** */

	public void setDo_MI(boolean do_MI)
	{
		this.do_MI = do_MI;
	}

	public void setDo_MI_KDE(boolean do_MI_KDE)
	{
		this.do_MI_KDE = do_MI_KDE;
	}

	public void setDo_linear_kernel(boolean do_linear_kernel)
	{
		this.do_linear_kernel = do_linear_kernel;
	}

	public void setDo_Degree_2_Poly(boolean do_Degree_2_Poly)
	{
		this.do_Degree_2_Poly = do_Degree_2_Poly;
	}

	public void setDo_Degree_3_Poly(boolean do_Degree_3_Poly)
	{
		this.do_Degree_3_Poly = do_Degree_3_Poly;
	}

	public void setDo_Degree_4_Poly(boolean do_Degree_4_Poly)
	{
		this.do_Degree_4_Poly = do_Degree_4_Poly;
	}

	public void setDo_XY_Poly(boolean do_XY_Poly)
	{
		this.do_XY_Poly = do_XY_Poly;
	}

	public void setDo_Euclidean(boolean do_Euclidean)
	{
		this.do_Euclidean = do_Euclidean;
	}

	public void setDo_Mahalanobis(boolean do_Mahalanobis)
	{
		this.do_Mahalanobis = do_Mahalanobis;
	}

	public void setDo_Gaussian(boolean do_Gaussian)
	{
		this.do_Gaussian = do_Gaussian;
	}

	public void setDo_Sigmoid(boolean do_Sigmoid)
	{
		this.do_Sigmoid = do_Sigmoid;
	}

	public void setDo_Log(boolean do_Log)
	{
		this.do_Log = do_Log;
	}

	public void setDo_Circular(boolean do_Circular)
	{
		this.do_Circular = do_Circular;
	}

	public void setDo_Cauchy(boolean do_Cauchy)
	{
		this.do_Cauchy = do_Cauchy;
	}

	private void create_Directory(String dir)
	{
		exist = new File(dir).exists();
		if (!exist) success = (new File(dir)).mkdirs();
	}

	/* *********************************** DRIVER ************************************************* */

	public void kPCA_Driver()
	{
		// System.out.println("\t\tPerforming Kernel PCA on the PCA Reduced Centered Data using the Top " + NUMBER_PCs + " modes.");
		// System.out.println("\t\t\tThe dimension of the Data is " + ROWS + " by " + KERNEL_FRAMES);
		// System.out.println("\t\t\tThe dimension of the Covariance matrix is " + ROWS + " by " + ROWS);
		// System.out.println("\t\t\tThe dimension of the Kernel matrix is " + KERNEL_FRAMES + " by " + KERNEL_FRAMES);
		// System.out.println("\t\t\tThe output directory is " + out_dir);

		if (do_linear_kernel) do_Linear_Kernel();
		if (do_Degree_2_Poly) do_Degree_2_Poly();
		if (do_Degree_3_Poly) do_Degree_3_Poly();
		if (do_Degree_4_Poly) do_Degree_4_Poly();
		if (do_XY_Poly) do_Poly_XY();
		if (do_Euclidean) do_Euclidean();
		if (do_Mahalanobis) do_Mahalanobis();
		if (do_Gaussian) do_Gaussian();
		if (do_Sigmoid) do_Sigmoid();
		if (do_Log) do_Log();
		if (do_Circular) do_Circular();
		if (do_Cauchy) do_Cauchy();

		if (do_MI || do_MI_KDE)
			{
				MI = new Mutual_Information(out_dir, data_reduced_transposed_centered);
				MI.get_Marginal_PROBS();
				MI.get_Joint_Probs_Table();
				noSolution = MI.is_No_Solution();
			}

		if (do_MI)
			{
				do_MI_1();
				if (!noSolution) do_MI_2();
			}

		if (do_MI_KDE)
			{
				if (!noSolution) do_MI_KDE();
				if (!noSolution) do_MI_KDE_2();
			}
	}
}
