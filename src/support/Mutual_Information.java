package support;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import Jama.Matrix;
import supportIO.Input_Parameters;
import supportKDE.KDE_2D;
import supportStats.EstimatePDF;

public class Mutual_Information
{
	public boolean no_solution;
	String title;
	final String out_dir;
	final int ROWS, COLS;
	int number_of_bins;
	double CO_EFF, kernelResolution, cellSize;
	final double MACH_EPS = Input_Parameters.FLOOR;
	final double PI2 = (2d * (Math.PI));
	final double shrinkage = Input_Parameters.KERNEL_SHRINKAGE;
	int number_of_points, KDE_offset;
	double[] X, Y, bandwidths;
	final Matrix marginalProbs, prob_table, jointProbsXX, jointProbsYY, jointProbsXY, PROBS, probabilities_X, probabilities_Y;
	final Matrix data, data_reduced, MI_KERNEL, MI_KERNEL2, cov, inv_cov;
	Matrix mutual_information_rows, mutual_information_cols;
	ArrayList<Double> X_values, MINS, MAXES, MEANS, STD_DEVS;
	final PCA pca;
	final Descriptive_Stats ds;

	/* ************************************* CONSTRUCTOR ************************************************************************************** */

	public Mutual_Information(String dir, Matrix m)
	{
		this.data = m;
		this.out_dir = dir;

		this.ROWS = 2;
		this.COLS = data.getColumnDimension();
		this.data_reduced = data.getMatrix(0, 1, 0, COLS - 1); // Gets the first 2 rows

		this.X_values = new ArrayList<Double>();
		this.MINS = new ArrayList<Double>();
		this.MAXES = new ArrayList<Double>();

		this.MI_KERNEL = new Matrix(COLS, COLS);
		this.MI_KERNEL2 = new Matrix(COLS, COLS);
		this.PROBS = new Matrix(ROWS, COLS);
		this.prob_table = new Matrix(COLS, COLS);
		this.probabilities_X = new Matrix(1, COLS);
		this.probabilities_Y = new Matrix(1, COLS);

		this.marginalProbs = new Matrix(2, COLS);
		this.jointProbsXX = new Matrix(COLS, COLS);
		this.jointProbsYY = new Matrix(COLS, COLS);
		this.jointProbsXY = new Matrix(COLS, COLS);

		this.ds = new Descriptive_Stats();
		this.pca = new PCA(data_reduced);

		Matrix Q = pca.get_Covariance_Matrix_by_Definition();
		cov = pca.Shrink_COV_Matrix(Q, shrinkage);
		inv_cov = cov.inverse();

		double det = cov.det();
		double pow = Math.pow(PI2, 2);
		double norm = Math.sqrt(pow * det);
		CO_EFF = Math.pow(norm, -1);
		if (CO_EFF > 10 | CO_EFF < MACH_EPS | CO_EFF == Double.NaN) CO_EFF = .5d;
	}

	/* *************************************** METHODS **************************************************************************************** */

	public Matrix get_MI_Kernel_2D_KDE() // FOR USE WITH 2 PCs ONLY: For this case, the MI double sum has 4 terms
	{
		if (!no_solution) get_Joint_Probs_KDE();

		int i, j;
		double pX1, pY1, pX2, pY2, pXX, pYY, pXY, pYX;
		double MI = 0;
		double arg1, arg2, arg3, arg4;
		double term1, term2, term3, term4;
		double sum = 0;
		for (i = 0; i < COLS; i++)
			{
				pX1 = marginalProbs.get(0, i);
				pY1 = marginalProbs.get(1, i);

				for (j = i; j < COLS; j++)
					{
						pX2 = marginalProbs.get(0, j);
						pY2 = marginalProbs.get(1, j);

						pXX = jointProbsXX.get(i, j);
						pYY = jointProbsYY.get(i, j);
						pXY = jointProbsXY.get(i, j);
						pYX = jointProbsXY.get(j, i);

						arg1 = (Math.log(pXX) - Math.log(pX1) - Math.log(pX2));
						arg2 = (Math.log(pXY) - Math.log(pX1) - Math.log(pY2));
						arg3 = (Math.log(pYX) - Math.log(pY1) - Math.log(pX2));
						arg4 = (Math.log(pYY) - Math.log(pY1) - Math.log(pY2));

						term1 = pXX * (arg1);
						term2 = pXY * (arg2);
						term3 = pYX * (arg3);
						term4 = pYY * (arg4);

						sum = (term1 + term2 + term3 + term4);
						MI = sum;
						double MI2 = (1d - sum);

						MI_KERNEL.set(i, j, MI);
						MI_KERNEL.set(j, i, MI);
						if (i == j) MI_KERNEL.set(j, i, 1d);

						MI_KERNEL2.set(i, j, MI2);
						MI_KERNEL2.set(j, i, MI2);
						if (i == j) MI_KERNEL2.set(j, i, 1d);
					}
			}
		return MI_KERNEL;
	}

	public Matrix get_MI_Kernel_2D_KDE_Version_2()
	{
		return MI_KERNEL2;
	}

	public Matrix get_MI_Kernel_Version1() // FOR USE WITH 2 PCs ONLY
	{
		for (int i = 0; i < COLS; i++)
			{
				Matrix X = prob_table.getMatrix(0, COLS - 1, i, i);
				Matrix Y = prob_table.getMatrix(i, i, 0, COLS - 1);
				double[] x_vals = X.getColumnPackedCopy();
				double[] y_vals = Y.getRowPackedCopy();
				double sum_X = 0;
				double sum_Y = 0;
				for (double x : x_vals)
					{
						double y = y_vals[i];
						sum_Y += y;
						sum_X += x;
					}
				probabilities_X.set(0, i, sum_X / COLS);
				probabilities_Y.set(0, i, sum_Y / COLS);
			}
		double MI, Probability_x, Probability_y, Probability_xy;
		for (int i = 0; i < COLS; i++)
			{
				Probability_y = probabilities_Y.get(0, i);
				for (int j = 0; j < COLS; j++)
					{
						Probability_x = probabilities_X.get(0, j);
						Probability_xy = prob_table.get(i, j);
						if (Probability_x < MACH_EPS) Probability_x = MACH_EPS;
						if (Probability_y < MACH_EPS) Probability_y = MACH_EPS;
						if (Probability_xy < MACH_EPS) Probability_xy = MACH_EPS;
						double arg = (Math.log(Probability_xy) - Math.log(Probability_x) - Math.log(Probability_y));
						MI = (Probability_xy * arg);
						MI_KERNEL.set(i, j, MI);
						MI_KERNEL.set(j, i, MI);
						if (i == j) MI_KERNEL.set(i, j, 1d);
					}
			}
		return MI_KERNEL;
	}

	public Matrix get_MI_Kernel_Version2() // FOR USE WITH 2 PCs ONLY
	{
		/* Get the Joint Probabilities of Frame i and Frame j: */
		/* The joint probabilities calculated using p(Frame i) * p(Frame j) */
		/* Where p(Frame i) and p(Frame j) are calculated from the ensemble covariance matrix */

		int i, j;
		double pX, pY, pXY;
		double MI = 0;
		double arg;
		double term;

		for (i = 0; i < COLS; i++)
			{
				pX = marginalProbs.get(0, i);
				pY = marginalProbs.get(1, i);

				for (j = i; j < COLS; j++)
					{
						pXY = prob_table.get(i, j);
						arg = (Math.log(pXY) - Math.log(pX) - Math.log(pY));
						term = pXY * (arg);
						// MI = (1d - term);
						MI = (term);

						MI_KERNEL.set(i, j, MI);
						MI_KERNEL.set(j, i, MI);
						if (i == j) MI_KERNEL.set(j, i, 1d);
					}
			}
		return MI_KERNEL;
	}

	public Matrix get_MI_Kernel_Version3() // FOR USE WITH 2 PCs ONLY: MI Approximated by joint probabilities.
	{
		/* No Marginal Probabilities in this version: */
		/* Get the Joint Probabilities of Frame i and Frame j: */
		/* The joint probabilities calculated using p(Frame i) * p(Frame j) */
		/* Where p(Frame i) and p(Frame j) are calculated from the ensemble covariance matrix */
		int i, j;
		double pXY, MI, arg, term;
		for (i = 0; i < COLS; i++)
			{
				for (j = i; j < COLS; j++)
					{
						pXY = prob_table.get(i, j);
						arg = (Math.log(pXY));
						term = pXY * (arg);
						MI = (term);
						// MI = (1d - term);

						MI_KERNEL.set(i, j, MI);
						MI_KERNEL.set(j, i, MI);
						// if (i == j) MI_KERNEL.set(j, i, 1d);
					}
			}
		return MI_KERNEL;
	}

	public Matrix get_PROBS()
	{
		for (int i = 0; i < ROWS; i++)
			{
				Matrix row = data.getMatrix(i, i, 0, COLS - 1);
				double[] row_sample = row.getRowPackedCopy();
				double mu = ds.get_mean(row_sample);
				double sigma = ds.get_standard_deviation(row_sample);
				double[] probs = ds.get_probabilities(row_sample, mu, sigma);
				Matrix row_probs = new Matrix(probs, 1);
				PROBS.setMatrix(i, i, 0, COLS - 1, row_probs);
			}
		// System.out.println("PROBS");
		// PROBS.print(9, 6);
		return PROBS;
	}

	public Matrix get_Joint_Probs_Table()
	{
		for (int i = 0; i < COLS; i++)
			{
				Matrix X = data_reduced.getMatrix(0, ROWS - 1, i, i); // Gets Frame i
				Matrix argument_X = X.transpose().times(inv_cov).times(X);
				double arg_X = (-0.500) * argument_X.get(0, 0);
				double prob_X = (CO_EFF * Math.exp(arg_X) / COLS); // Probability of Frame i

				for (int j = i; j < COLS; j++) // Gets Frame j
					{
						Matrix Y = data_reduced.getMatrix(0, ROWS - 1, j, j);
						Matrix argument_Y = Y.transpose().times(inv_cov).times(Y);
						double arg_Y = (-0.500) * argument_Y.get(0, 0);
						double prob_Y = (CO_EFF * Math.exp(arg_Y) / COLS); // Probability of Frame j
						double prob_XY = ((prob_X * prob_Y)); // Probability of Frame i and Frame j
						prob_table.set(i, j, prob_XY);
						prob_table.set(j, i, prob_XY);
					}
			}
		// System.out.println("jointProbsTable");
		// prob_table.print(9, 6);
		return prob_table;
	}

	public Matrix get_Marginal_PROBS()
	{
		/* Get the Marginal Probabilities of Variables X and Y */
		/* Each ROW is a variable: 'marginalProbs' */
		/* The marginal probabilities calculated using PDF(x) * dx */
		for (int i = 0; i < 2; i++)
			{
				Matrix X = data_reduced.getMatrix(i, i, 0, COLS - 1);
				ArrayList<Double> x = new ArrayList<Double>();
				for (int index = 0; index < COLS; index++)
					{
						x.add(X.get(0, index));
					}
				title = "MI_Input_PC_" + (i + 1);

				EstimatePDF est = new EstimatePDF(out_dir, title, x);
				ArrayList<Double> mpx = est.estimate_ALD();

				if (mpx.size() == 0 || mpx.equals(null))
					{
						System.err.println("Warning: PDF Estimator was not able to compute a probability distribution for the data.");
						System.err.println("Some MI kernels will not be available...");
						no_solution = true;
						break;
					}

				double[] marginalX = new double[COLS];
				for (int index = 0; index < COLS; index++)
					{
						double prob = mpx.get(index);
						marginalX[index] = prob;
					}

				Matrix pX = new Matrix(marginalX, 1);
				marginalProbs.setMatrix(i, i, 0, COLS - 1, pX);
			}
		// System.out.println("marginalProbs");
		// marginalProbs.print(9, 6);

		return marginalProbs;
	}

	public void get_Joint_Probs_KDE()
	{
		/* Get the Joint Probabilities of Variables 'X' and 'X' --> 'jointProbsXX' */
		/* The joint probabilities calculated using PDF(XX) * dxdx */
		/* The dx and dy values depend on the kernel bandwidths */
		for (int i = 0; i < 2 - 1; i++)
			{
				Matrix rowX = data_reduced.getMatrix(i, i, 0, COLS - 1);
				X = rowX.getRowPackedCopy();
				double[] kde_arrayX = new double[2 * COLS];
				for (int k = 0; k < COLS; k++)
					{
						kde_arrayX[k + k] = X[k];
						kde_arrayX[k + k + 1] = X[k];
					}
				KDE_offset = 0;
				number_of_points = COLS;

				KDE_2D kde = new KDE_2D(kde_arrayX, KDE_offset, number_of_points, null, null, null);

				Arrays.sort(X);
				double minX = X[0];
				double maxX = X[COLS - 1];
				kde.computeBoundsReal(minX, maxX, minX, maxX);
				kde.compute();

				double dxdy = ((1d / kernelResolution) * (1d / kernelResolution));

				for (int kX = 0; kX < COLS; kX++)
					{
						for (int kY = kX; kY < COLS; kY++)
							{
								double prob = kde.getProb(X[kX], X[kY]);
								jointProbsXX.set(kX, kY, prob * dxdy);
								jointProbsXX.set(kY, kX, prob * dxdy);
							}
					}
				/* Get the Joint Probabilities of Variables 'Y' and 'Y' --> 'jointProbsYY' */
				/* The joint probabilities calculated using PDF(YY) * dydy */
				for (int j = i + 1; j < 2; j++)
					{
						Matrix rowY = data_reduced.getMatrix(j, j, 0, COLS - 1);
						Y = rowY.getRowPackedCopy();
						double[] kde_arrayY = new double[2 * COLS];
						for (int k = 0; k < COLS; k++)
							{
								kde_arrayY[k + k] = Y[k];
								kde_arrayY[k + k + 1] = Y[k];
							}
						KDE_offset = 0;
						number_of_points = COLS;

						Arrays.sort(Y);
						double minY = Y[0];
						double maxY = Y[COLS - 1];

						kde = new KDE_2D(kde_arrayY, KDE_offset, number_of_points, null, null, null);
						kde.computeBoundsReal(minY, maxY, minY, maxY);
						kde.compute();

						for (int kX = 0; kX < COLS; kX++)
							{
								for (int kY = kX; kY < COLS; kY++)
									{
										double prob = kde.getProb(Y[kX], Y[kY]);
										jointProbsYY.set(kX, kY, prob * dxdy);
										jointProbsYY.set(kY, kX, prob * dxdy);
									}
							}
						/* Get the Joint Probabilities of Variables 'X' and 'Y' --> 'jointProbsXY' */
						/* The joint probabilities calculated using PDF(XY) * dxdy */
						double[] kde_arrayXY = new double[2 * COLS];

						for (int k = 0; k < COLS; k++)
							{
								kde_arrayXY[k + k] = X[k];
								kde_arrayXY[k + k + 1] = Y[k];
							}
						KDE_offset = 0;
						number_of_points = COLS;

						kde = new KDE_2D(kde_arrayXY, KDE_offset, number_of_points, null, null, null);
						kde.computeBoundsReal(minX, maxX, minY, maxY);
						kde.compute();

						bandwidths = kde.get_band_Widths();
						cellSize = kde.get_cell_Size();

						for (int kX = 0; kX < COLS; kX++)
							{
								for (int kY = kX; kY < COLS; kY++)
									{
										double prob = kde.getProb(X[kX], Y[kY]);
										jointProbsXY.set(kX, kY, prob * dxdy);
										jointProbsXY.set(kY, kX, prob * dxdy);
									}
							}
					}
			}
		// System.out.println("joint probs");
		// jointProbsXY.print(20, 12);
	}

	public Matrix get_Mutual_Information_Rows()
	{
		number_of_bins = COLS;
		if (COLS > 100) number_of_bins = 100;

		mutual_information_rows = new Matrix(ROWS, ROWS);
		for (int i = 0; i < ROWS; i++)
			{
				for (int j = 0; j < COLS; j++)
					{
						double d = data.get(i, j);
						if (Math.abs(d) < MACH_EPS) // if there are zeros in the data matrix, set them to machine epsilon
							{
								if (d > 0)
									{
										d = MACH_EPS;
										data.set(i, j, d);
									}
								if (d < 0)
									{
										d = -MACH_EPS;
										data.set(i, j, d);
									}
							}
						X_values.add(d);
					}
				Collections.sort(X_values); // sort the rows and extract the min and max value
				MINS.add(X_values.get(0));
				MAXES.add(X_values.get(X_values.size() - 1));
				// System.out.println("X_MAX from sort: " + X_values.get(X_values.size() - 1) + "\n");
				// System.out.println("X_MIN from sort: " + X_values.get(0) + "\n");
				X_values.clear();
			}

		for (int i = 0; i < ROWS; i++) // GET X LOOP: iterate over all rows of variable samples
			{
				Matrix X_data = data.getMatrix(i, i, 0, COLS - 1);
				double X_MAX = MAXES.get(i);
				double X_MIN = MINS.get(i);
				// System.out.println("X_MAX: " + X_MAX + "\n");
				// System.out.println("X_MIN: " + X_MIN + "\n");
				double X_bin_size = ((X_MAX - X_MIN) / (number_of_bins - 1)); // determine the bin size for X
				double[] X_freq = new double[number_of_bins]; // define the frequency bins for X

				for (int k = 0; k < COLS; k++) // iterate over all the samples of X
					{
						double X_element = X_data.get(0, k);
						double X_diff = X_element - X_MIN;
						double X_distance_from_MIN = (X_diff / (X_bin_size));
						int X_num_of_bins = (int) (X_distance_from_MIN);
						if (X_distance_from_MIN - X_num_of_bins >= .5) X_num_of_bins++;
						if (X_diff == 0) X_num_of_bins = 0;
						X_freq[X_num_of_bins] = X_freq[X_num_of_bins] + 1;
					}
				double[] X_probabilities = new double[number_of_bins];

				for (int z = 0; z < X_freq.length; z++) // turn the frequencies into probabilities
					{
						X_probabilities[z] = (X_freq[z]) / COLS;
					}

				Matrix X = new Matrix(X_probabilities, 1);
				// X.print(5, 5);

				for (int j = i; j < ROWS; j++) // GET Y LOOP: iterate over all rows of variable samples
					{
						Matrix Y_data = data.getMatrix(j, j, 0, COLS - 1);
						double Y_MAX = MAXES.get(j);
						double Y_MIN = MINS.get(j);
						// System.out.println("Y_MAX: " + Y_MAX + "\n");
						// System.out.println("Y_MIN: " + Y_MIN + "\n");
						double Y_bin_size = ((Y_MAX - Y_MIN) / (number_of_bins - 1)); // determine the bin size for Y
						double[] Y_freq = new double[number_of_bins]; // define the frequency bins for Y

						for (int k = 0; k < COLS; k++) // iterate over all the samples of Y
							{
								double Y_element = Y_data.get(0, k);
								double Y_diff = Y_element - Y_MIN;
								double Y_distance_from_MIN = (Y_diff / (Y_bin_size));
								int Y_num_of_bins = (int) (Y_distance_from_MIN);
								if (Y_distance_from_MIN - Y_num_of_bins >= .5) Y_num_of_bins++;
								if (Y_diff == 0) Y_num_of_bins = 0;
								Y_freq[Y_num_of_bins] = Y_freq[Y_num_of_bins] + 1;
							}

						double[] Y_probabilities = new double[number_of_bins];

						for (int z = 0; z < X_freq.length; z++) // turn the frequencies into probabilities
							{
								Y_probabilities[z] = (Y_freq[z]) / COLS;
							}

						Matrix Y = new Matrix(Y_probabilities, 1);
						// Y.print(5, 5);
						Matrix joint_prob_table = (X.transpose()).times(Y); // create the joint probability table
						// joint_prob_table.print(5, 5);
						double Probability_x, Probability_y, Probability_xy = MACH_EPS;
						double MI = 0;
						for (int p = 0; p < number_of_bins; p++) // get the values for the outer sum
							{
								Probability_y = Y.get(0, p);
								for (int q = 0; q < number_of_bins; q++) // get the values for the inner sum
									{
										Probability_x = X.get(0, q);
										Probability_xy = joint_prob_table.get(p, q);
										if (Probability_x < MACH_EPS) Probability_x = MACH_EPS;
										if (Probability_y < MACH_EPS) Probability_y = MACH_EPS;
										if (Probability_xy < MACH_EPS) Probability_xy = MACH_EPS;
										double log_argument = ((Probability_xy) / (Probability_x * Probability_y));
										double log_part = Math.log(log_argument);
										double argument = (Probability_xy * log_part);
										MI = (MI + argument); // calculate the double sum
										if (MI < 0.000001) MI = 0;
									}
							}
						mutual_information_rows.set(i, j, MI);
						mutual_information_rows.set(j, i, MI);
						MI = MACH_EPS;
					}
			}
		return mutual_information_rows;
	}

	public Matrix get_Mutual_Information_Cols()
	{

		number_of_bins = 5;

		mutual_information_cols = new Matrix(COLS, COLS);
		for (int i = 0; i < COLS; i++)
			{
				for (int j = 0; j < ROWS; j++)
					{
						double d = data.get(j, i);

						if (Math.abs(d) < MACH_EPS) // if there are zeros in the data matrix, set them to machine epsilon
							{
								if (d > 0)
									{
										d = MACH_EPS;
										data.set(j, i, d);
									}
								if (d < 0)
									{
										d = -MACH_EPS;
										data.set(j, i, d);
									}
							}
						X_values.add(d);
					}
				Collections.sort(X_values); // sort the rows and extract the min and max value
				MINS.add(X_values.get(0));
				MAXES.add(X_values.get(X_values.size() - 1));
				// System.out.println("X_MAX from sort: " + X_values.get(X_values.size() - 1) + "\n");
				// System.out.println("X_MIN from sort: " + X_values.get(0) + "\n");
				X_values.clear();
			}

		for (int i = 0; i < COLS; i++) // GET X LOOP: iterate over all COLS of variable samples
			{
				Matrix X_data = data.getMatrix(0, (ROWS - 1), i, i);
				double X_MAX = MAXES.get(i);
				double X_MIN = MINS.get(i);

				double X_bin_size = ((X_MAX - X_MIN) / (number_of_bins - 1)); // determine the bin size for X
				double[] X_freq = new double[number_of_bins]; // define the frequency bins for X

				// System.out.println("COL: " + i + "\n");
				// System.out.println("X_MAX: " + X_MAX + "\n");
				// System.out.println("X_MIN: " + X_MIN + "\n");
				// System.out.println("X_bin_size: " + X_bin_size + "\n");
				// System.out.println("number_of_bins: " + number_of_bins + "\n");

				for (int k = 0; k < ROWS; k++) // iterate over all the samples of X
					{
						double X_element = X_data.get(k, 0);
						double X_diff = X_element - X_MIN;
						double X_distance_from_MIN = (X_diff / (X_bin_size));
						int X_index_of_bin = (int) (X_distance_from_MIN);

						// System.out.println("ROW: " + k + "\n");
						// System.out.println("X_element: " + X_element + "\n");
						// System.out.println("X_diff: " + X_diff + "\n");
						// System.out.println("X_distance_from_MIN: " + X_distance_from_MIN + "\n");
						// System.out.println("X_index_of_bin: " + X_index_of_bin + "\n");

						if (X_distance_from_MIN - X_index_of_bin >= .5) X_index_of_bin++;
						if (X_diff == 0) X_index_of_bin = 0;
						X_freq[X_index_of_bin] = X_freq[X_index_of_bin] + 1;
					}
				double[] X_probabilities = new double[number_of_bins];

				for (int z = 0; z < X_freq.length; z++) // turn the frequencies into probabilities
					{
						X_probabilities[z] = (X_freq[z]) / ROWS;
					}

				Matrix X = new Matrix(X_probabilities, 1);
				// X.print(5, 5);

				for (int j = i; j < COLS; j++) // GET Y LOOP: iterate over all rows of variable samples
					{
						Matrix Y_data = data.getMatrix(0, (ROWS - 1), j, j);
						double Y_MAX = MAXES.get(j);
						double Y_MIN = MINS.get(j);
						// System.out.println("Y_MAX: " + Y_MAX + "\n");
						// System.out.println("Y_MIN: " + Y_MIN + "\n");
						double Y_bin_size = ((Y_MAX - Y_MIN) / (number_of_bins - 1)); // determine the bin size for Y
						double[] Y_freq = new double[number_of_bins]; // define the frequency bins for Y

						for (int k = 0; k < ROWS; k++) // iterate over all the samples of Y
							{
								double Y_element = Y_data.get(k, 0);
								double Y_diff = Y_element - Y_MIN;
								double Y_distance_from_MIN = (Y_diff / (Y_bin_size));
								int Y_num_of_bins = (int) (Y_distance_from_MIN);
								if (Y_distance_from_MIN - Y_num_of_bins >= .5) Y_num_of_bins++;
								if (Y_diff == 0) Y_num_of_bins = 0;
								Y_freq[Y_num_of_bins] = Y_freq[Y_num_of_bins] + 1;
							}

						double[] Y_probabilities = new double[number_of_bins];

						for (int z = 0; z < X_freq.length; z++) // turn the frequencies into probabilities
							{
								Y_probabilities[z] = (Y_freq[z]) / ROWS;
							}

						Matrix Y = new Matrix(Y_probabilities, 1);
						// Y.print(5, 5);
						Matrix joint_prob_table = (X.transpose()).times(Y); // create the joint probability table
						// prob_table.print(5, 5);
						double Probability_x, Probability_y, Probability_xy = MACH_EPS;
						double MI = 0;
						for (int p = 0; p < number_of_bins; p++) // get the values for the outer sum
							{
								Probability_y = Y.get(0, p);
								for (int q = 0; q < number_of_bins; q++) // get the values for the inner sum
									{
										Probability_x = X.get(0, q);
										Probability_xy = joint_prob_table.get(p, q);
										if (Probability_x < MACH_EPS) Probability_x = MACH_EPS;
										if (Probability_y < MACH_EPS) Probability_y = MACH_EPS;
										if (Probability_xy < MACH_EPS) Probability_xy = MACH_EPS;
										double log_argument = ((Probability_xy) / (Probability_x * Probability_y));
										double log_part = Math.log(log_argument);
										double argument = (Probability_xy * log_part);
										MI = (MI + argument); // calculate the double sum
										// if (MI < 0.000001)
										// MI = 0;
									}
							}
						// System.out.println("MI: " + MI + "\n");
						mutual_information_cols.set(j, i, MI);
						mutual_information_cols.set(i, j, MI);
						MI = MACH_EPS;
					}
			}
		// mutual_information_cols.print(6, 3);
		return mutual_information_cols;
	}

	public void set_kernel_resolution(double resolution)
	{
		kernelResolution = resolution;
	}

	public void set_CO_EFF(double co_eff)
	{
		CO_EFF = co_eff;
	}

	public double[] getBandwidths()
	{
		return bandwidths;
	}

	public double getCellSize()
	{
		return cellSize;
	}

	public boolean is_No_Solution()
	{
		return no_solution;
	}
}
