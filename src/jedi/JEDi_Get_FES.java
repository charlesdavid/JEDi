package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import Jama.Matrix;
import supportIO.DateUtils;
import supportIO.Input_Parameters;
import supportIO.Matrix_IO;
import supportKDE.KDE_2D;
import supportPlot.FES_Plot;

public class JEDi_Get_FES
{
	boolean exist, success;
	String out_dir, path;
	final int ROWS, COLS;
	final int OP1, OP2, number_of_points, OP_offset, KDE_offset;
	final double kernel_resolution = Input_Parameters.KDE_RESOLUTION;
	final double kernel_cell = Input_Parameters.KDE_CELL;
	final double margin = Input_Parameters.KDE_MARGIN;
	double cellsize;
	double[] KDE_bounds, KDE_bandwidths;
	final Matrix normed_projections, order_param1, order_param2;
	Matrix FE;
	File LOG;
	PrintWriter LOG_writer;
	NumberFormat nf;
	KDE_2D KDE;
	FES_Plot plot1;
	FES_Plot plot2;

	public JEDi_Get_FES(Matrix DVPs, int op1, int op2, int num_points, int offset, double size)
	{
		this.normed_projections = DVPs;
		this.OP1 = op1;
		this.OP2 = op2;
		this.number_of_points = num_points;
		this.OP_offset = offset;
		this.cellsize = size;
		this.KDE_offset = 0;

		ROWS = normed_projections.getRowDimension();
		COLS = normed_projections.getColumnDimension();
		int row_index1 = OP_offset;
		int row_index2 = OP_offset + number_of_points - 1;
		if (row_index2 > ROWS - 1) row_index2 = ROWS - 1;
		order_param1 = normed_projections.getMatrix(row_index1, row_index2, OP1, OP1);
		order_param2 = normed_projections.getMatrix(row_index1, row_index2, OP2, OP2);

		nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setMinimumFractionDigits(3);
		nf.setRoundingMode(RoundingMode.HALF_UP);
	}

	public void get_FES()
	{
		/* Get the DVPs to use as order parameters in the deltaG free energy calculations */
		if (OP_offset + number_of_points > ROWS) System.err.println("ERROR! The OP offset plus the number of points exceeds the length of the PCs!");
		double[] order_parameter_1 = order_param1.getColumnPackedCopy();
		List<Double> order_parameter_1_sorted = new ArrayList<Double>();
		for (double d : order_parameter_1)
			{
				order_parameter_1_sorted.add(d);
			}
		Collections.sort(order_parameter_1_sorted);
		double[] order_parameter_2 = order_param2.getColumnPackedCopy();
		List<Double> order_parameter_2_sorted = new ArrayList<Double>();
		for (double d : order_parameter_2)
			{
				order_parameter_2_sorted.add(d);
			}
		Collections.sort(order_parameter_2_sorted);

		/* 2D KDE using Gaussian Functions */
		int length = order_parameter_1.length;
		double[] kde_array = new double[2 * length];
		for (int i = 0; i < length; i++)
			{
				kde_array[i + i] = order_parameter_1[i];
				kde_array[i + i + 1] = order_parameter_2[i];
			}
		KDE = new KDE_2D(kde_array, KDE_offset, number_of_points, null, null, null);

		int size = order_parameter_1_sorted.size();
		double min1 = order_parameter_1_sorted.get(0);
		double max1 = order_parameter_1_sorted.get(size - 1);
		double min2 = order_parameter_2_sorted.get(0);
		double max2 = order_parameter_2_sorted.get(size - 1);

		KDE.computeBoundsReal(min1, max1, min2, max2);
		KDE.compute();
		KDE_bandwidths = KDE.get_band_Widths();
		cellsize = KDE.get_cell_Size();
		KDE_bounds = KDE.getBounds();

		double[] probabilities = new double[length];
		double[] probabilities_sorted = new double[length];
		for (int i = 0; i < length; i++)
			{
				double prob = KDE.getProb(order_parameter_1[i], order_parameter_2[i]);
				probabilities[i] = prob;
				probabilities_sorted[i] = prob;
			}
		Arrays.sort(probabilities_sorted);
		final double prob_max = probabilities_sorted[length - 1];
		final double ln_prob_max = Math.log(prob_max);
		final double KBT = (-0.600); // Units are in kcal/mol, T = 300K
		FE = new Matrix(length, 3);
		for (int i = 0; i < length; i++)
			{
				double prob = probabilities[i];
				if (prob < 1.00E-16) prob = 1.00E-16;
				double ln_prob = Math.log(prob);
				double delta_G = KBT * (ln_prob - ln_prob_max);
				FE.set(i, 0, order_parameter_1[i]);
				FE.set(i, 1, order_parameter_2[i]);
				FE.set(i, 2, delta_G);
			}
		if (Input_Parameters.verbose) Matrix_IO.write_Matrix(FE, path, 20, 12);
		plot1 = new FES_Plot(out_dir, "Free_Energy_Landscape_Above_Left", "Derived from Normed Projections of PC Modes 1 & 2", "PC " + (OP1 + 1), "FreeEnergy", "PC " + (OP2 + 1),
				FE, "LEFT");
		plot2 = new FES_Plot(out_dir, "Free_Energy_Landscape_Above_Right", "Derived from Normed Projections PC Modes 1 & 2", "PC " + (OP1 + 1), "FreeEnergy", "PC " + (OP2 + 1), FE,
				"RIGHT");
	}

	public void write_FES_Log()
	{
		LOG = new File(out_dir + "FES_Log.txt");
		try
			{
				LOG_writer = new PrintWriter(new BufferedWriter(new FileWriter(LOG)));
			}
		catch (IOException e)
			{
				System.err.println("Could not create the Job Log file.");
				e.printStackTrace();
			}
		LOG_writer.write("Total Number of Conformations (rows) used to create the DVPs: " + ROWS + "\n");
		LOG_writer.write("Total Number of Modes (cols) used to create the DVPs: " + COLS + "\n");
		LOG_writer.write("Order Parameter 1 (Col# 1): " + (OP1 + 1) + "\n");
		LOG_writer.write("Order Parameter 2 (Col# 2): " + (OP2 + 1) + "\n");
		LOG_writer.write("Number of points to extract from the OPs to use for KDE: " + number_of_points + "\n");
		LOG_writer.write("Offset in selecting points from the order parameters: " + OP_offset + "\n");
		LOG_writer.write("The bounds for the KDE are: OP1(" + nf.format(KDE_bounds[0]) + "," + nf.format(KDE_bounds[2]) + "); OP2(" + nf.format(KDE_bounds[1]) + ","
				+ nf.format(KDE_bounds[3]) + ")" + "\n");
		LOG_writer.write("The kernel bandwidths are: " + nf.format(KDE_bandwidths[0]) + "   " + nf.format(KDE_bandwidths[3]) + "\n\n");
		LOG_writer.write("The kernel cell-size is: " + nf.format(cellsize) + "\n\n");
		LOG_writer.write(String.format("%-12s%-12s%-12s%-12s", " ", "OP1", "OP2", "FE"));
		LOG_writer.flush();
		FE.print(LOG_writer, 12, 3);
		LOG_writer.flush();
		LOG_writer.write("\nAnalysis completed: " + DateUtils.now());
		LOG_writer.close();
	}

	public double[] get_KDE_bandwidths()
	{
		return KDE_bandwidths;
	}

	public Matrix get_FE()
	{
		return FE;
	}

	public void set_Out_dir(String dir)
	{
		this.out_dir = dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
		path = out_dir + "FES_Matrix_" + (OP1 + 1) + "_" + (OP2 + 1) + ".txt";
	}
}
