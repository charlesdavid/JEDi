package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;

import Jama.Matrix;
import support.Subspace_Analysis2;
import supportIO.DateUtils;
import supportIO.Input_Parameters;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportPlot.RMSIP_Plot;

/**
 * JED class JED_Get_Subspace_Analysis: Performs Subspace Analysis on sets of eigenvectors with varying degree of outputs.
 * 
 * Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 *
 */

public class JEDi_Get_Subspace_Analysis
{
	boolean exist, success, doRandom, verbose;
	final int ROWS, COLS;
	String out, out_dir, date = DateUtils.now();
	final String type_subset, type_comparison;
	double RMSIP, max_angle;
	Matrix projections, CO_matrix_1_2, CO_matrix_2_1, PAs, avg_PAs, avg_RMSIPs, avg_COs, RMSIP_std_devs;
	final Matrix data1, data2;
	ArrayList<Double> cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products, vectorial_sum_of_angles, RMSIPs, Z_Scores;
	File Job_log;
	PrintWriter job_log_writer;
	final NumberFormat nf0, nf2, nf3, nf6;
	RoundingMode rm = RoundingMode.HALF_UP;

	/**
	 * Constructor for performing Subspace Analysis
	 * 
	 * Note: The matrices of eigenvectors MUST have the same dimensions for both rows and columns.
	 *
	 * @param subset     The subset
	 * @param comparison The kind of comparison: For example, Direct versus Hierarchical, COV versus CORR
	 * @param evects1    The first set of eigenvectors for the subspace comparison
	 * @param evects2    The second set of eigenvectors for the subspace comparison
	 */

	public JEDi_Get_Subspace_Analysis(String subset, String comparison, Matrix evects1, Matrix evects2)
	{
		this.verbose = Input_Parameters.verbose;
		this.type_subset = subset;
		this.type_comparison = comparison;
		this.data1 = evects1;
		this.data2 = evects2;
		this.ROWS = data1.getRowDimension();
		this.COLS = data1.getColumnDimension();

		if (data1.getRowDimension() != data2.getRowDimension())
			{
				System.err.println("Vector space dimensions do not match:");
				System.err.println("VS1 dim = " + data1.getRowDimension() + " VS2 dim = " + data2.getRowDimension());
			}

		this.nf0 = NumberFormat.getInstance();
		this.nf0.setRoundingMode(rm);
		this.nf0.setMaximumFractionDigits(0);
		this.nf0.setMinimumFractionDigits(0);

		this.nf2 = NumberFormat.getInstance();
		this.nf2.setRoundingMode(rm);
		this.nf2.setMaximumFractionDigits(2);
		this.nf2.setMinimumFractionDigits(2);

		this.nf3 = NumberFormat.getInstance();
		this.nf3.setRoundingMode(rm);
		this.nf3.setMaximumFractionDigits(3);
		this.nf3.setMinimumFractionDigits(3);

		this.nf6 = NumberFormat.getInstance();
		this.nf6.setRoundingMode(rm);
		this.nf6.setMaximumFractionDigits(6);
		this.nf6.setMinimumFractionDigits(6);
	}

	/**
	 * Calls the Subspace Analysis class and runs the SSA method: This compares the two subspaces directly using the metrics of RMSIP, PA, CO, CP, and VAS.
	 */
	public void get_SSA()
	{
		Subspace_Analysis2 ssa = new Subspace_Analysis2(data1, data2);
		ssa.get_SSA();

		RMSIP = ssa.getRMSIP();
		max_angle = ssa.getMax_angle();
		projections = ssa.getProjections();
		CO_matrix_1_2 = ssa.getCO_matrix_1_2();
		CO_matrix_2_1 = ssa.getCO_matrix_2_1();
		cumulative_overlaps_1_2 = ssa.getCumulative_overlaps_1_2();
		cumulative_overlaps_2_1 = ssa.getCumulative_overlaps_2_1();
		principle_angles_svd = ssa.getPrinciple_angles_svd();
		cosine_products = ssa.getCosine_products();

		double sum_PA_squared = 0;
		double vector_sum_of_angles = 0;
		vectorial_sum_of_angles = new ArrayList<>();
		for (double d : principle_angles_svd)
			{
				sum_PA_squared += (d * d);
				vector_sum_of_angles = Math.sqrt(sum_PA_squared);
				vectorial_sum_of_angles.add(vector_sum_of_angles);
			}

		out = out_dir + type_comparison + File.separatorChar;
		exist = new File(out).exists();
		if (!exist) success = (new File(out)).mkdirs();

		if (verbose)
			{
				Matrix_IO.write_Matrix(CO_matrix_1_2, out + "CO_1_2_dim_" + COLS + ".txt", 20, 12);
				Matrix_IO.write_Matrix(CO_matrix_2_1, out + "CO_2_1_dim_" + COLS + ".txt", 20, 12);
				Matrix_IO.write_Matrix(projections, out + "Projections_dim_" + COLS + ".txt", 20, 12);
				List_IO.write_Double_List(principle_angles_svd, out + "PAs_dim_" + COLS + ".txt", 12);
				List_IO.write_Double_List(cosine_products, out + "Cosine_Products_dim_" + COLS + ".txt", 12);
				List_IO.write_Double_List(vectorial_sum_of_angles, out + "Vector_Sums_of_Angles_dim_" + COLS + ".txt", 6);
			}
		write_SSA_Log();
	}

	private void write_SSA_Log()
	{
		Job_log = new File(out + "JED_SSA_Log_dim_" + COLS + ".txt");
		try
			{
				job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
			}
		catch (IOException e)
			{
				System.err.println("Problem writing the log file...Check permissions... " + out + "JED_SSA_Log_dim_" + COLS + ".txt");
				e.printStackTrace();
			}
		job_log_writer.write("Type of analysis: " + type_subset + "\n");
		job_log_writer.write("Comparison: " + type_comparison + "\n");
		job_log_writer.write("First Set of Eigenvectors:\n");
		job_log_writer.write("Rows: " + ROWS + "\n");
		job_log_writer.write("Cols: " + COLS + "\n");
		job_log_writer.write("Second Set of Eigenvectors:\n");
		job_log_writer.write("Rows: " + ROWS + "\n");
		job_log_writer.write("Cols: " + COLS + "\n");
		job_log_writer.write("\nOutput Directory: " + out + "\n");
		if (verbose)
			{
				job_log_writer.write("Projections file written to: " + "Projections_dim_" + COLS + ".txt\n");
				job_log_writer.write("Cumulative overlaps 1 --> 2 file written to: " + "CO_1_2_dim_" + COLS + ".txt\n");
				job_log_writer.write("Cumulative overlaps 2 --> 1 file written to: " + "CO_2_1_dim_" + COLS + ".txt\n");
				job_log_writer.write("Principle Angles file written to: " + "PAs_dim_" + COLS + ".txt" + "\n");
				job_log_writer.write("Cosine Products file written to: " + "Cosine_Products_dim_" + COLS + ".txt\n");
				job_log_writer.write("Vectorial sums of angles file written to: " + "Vector_Sums_of_Angles_dim_" + COLS + ".txt\n");
			}
		job_log_writer.write("\nThe Inner Products of each vector in subspace 1 with each vector in subspace 2 are:\n");
		projections.print(job_log_writer, 9, 3);
		job_log_writer.write("\nThe cumulative overlaps CO_" + COLS + " for each vector in subspace 1 with all the vectors in subspace 2 are:\n\n");
		int j = 1;
		for (double d : cumulative_overlaps_1_2)
			{
				job_log_writer.write(String.format("%-8s%-12s%-12s%n", "Vector ", j, nf3.format(d)));
				j++;
			}
		job_log_writer.write("\nThe cumulative overlaps CO_" + COLS + " for each vector in subspace 2 with all the vectors in subspace 1 are:\n/n");
		j = 1;
		for (double d : cumulative_overlaps_2_1)
			{
				job_log_writer.write(String.format("%-8s%-12s%-12s%n", "Vector ", j, nf3.format(d)));
				j++;
			}
		job_log_writer.write("\nThe RMSIP score is " + nf6.format(RMSIP) + "\n");
		job_log_writer.write("\nThe principle angles (in degrees) are: " + "\n\n");
		int i = 1;
		for (double p : principle_angles_svd)
			{
				job_log_writer.write(String.format("%-6s%-12s%-12s%n", "PA", i, nf0.format(p)));
				i++;
			}
		job_log_writer.write("\nThe cosine products (in degrees) are: " + "\n\n");
		i = 1;
		for (double c : cosine_products)
			{
				job_log_writer.write(String.format("%-6s%-12s%-12s%n", "CP", i, nf0.format(Math.acos(c) * 180 / Math.PI)));
				i++;
			}
		job_log_writer.write("\nThe vectorial sums of angles (in degrees) are: " + "\n/n");
		i = 1;
		for (double p : vectorial_sum_of_angles)
			{
				job_log_writer.write(String.format("%-6s%-12s%-12s%n", "VS", i, nf0.format(p)));
				i++;
			}
		job_log_writer.write("\nMaximum possible angle between two subspaces of this dimension is " + nf0.format(max_angle) + " degrees\n\n");
		job_log_writer.write("Analysis completed: " + date);
		job_log_writer.flush();
		job_log_writer.close();
	}

	/**
	 * Calls the Subspace Analysis class and runs the FSSA Iterated method: This iteratively compares the two subspaces using the metrics of RMSIP, PA, CO, CP, and VAS. Comparisons
	 * are made for SS dims from 1 to the entered dim.
	 */
	public void get_FSSA_Iterated()
	{
		Subspace_Analysis2 fSSAi = new Subspace_Analysis2(data1, data2);
		fSSAi.get_fast_SSA_iterated();

		RMSIPs = fSSAi.getRMSIPs();
		PAs = fSSAi.getPAs();
		avg_RMSIPs = fSSAi.getAvg_RMSIP_Score_Matrix();
		avg_PAs = fSSAi.getAvg_PA_Matrix();
		avg_COs = fSSAi.getAvg_CO_Score_Matrix();
		RMSIP_std_devs = fSSAi.getRMSIP_Std_Dev_Matrix();

		Z_Scores = new ArrayList<Double>();
		Matrix rmsipPlot2 = new Matrix(RMSIPs.size(), 4);

		int i = 0;
		for (double d : RMSIPs)
			{
				double meanRMSIP = avg_RMSIPs.get(0, i);
				double sd = RMSIP_std_devs.get(0, i);
				double zscore = ((d - meanRMSIP) / RMSIP_std_devs.get(0, i));
				Z_Scores.add(zscore);

				rmsipPlot2.set(i, 0, d);
				rmsipPlot2.set(i, 1, meanRMSIP);
				rmsipPlot2.set(i, 2, sd);
				rmsipPlot2.set(i, 3, zscore);

				i++;
			}

		out = out_dir + type_comparison + File.separatorChar;
		exist = new File(out).exists();
		if (!exist) success = (new File(out)).mkdirs();

		if (verbose)
			{
				String name = "Avg_Random_RMSIPs.txt";
				Matrix_IO.write_Matrix(avg_RMSIPs.transpose(), out, name, 20, 12);
				name = "Avg_Random_PAs.txt";
				Matrix_IO.write_Matrix(avg_PAs.transpose(), out, name, 6, 0);
				name = "Avg_Random_COs.txt";
				Matrix_IO.write_Matrix(avg_COs.transpose(), out, name, 12, 6);
				name = "Random_RMSIP_Std_Devs.txt";
				Matrix_IO.write_Matrix(RMSIP_std_devs.transpose(), out, name, 20, 12);
				name = "Iterated_PAs.txt";
				Matrix_IO.write_Matrix(PAs, out, name, 9, 3);
				name = "Iterated_RMSIPs.txt";
				List_IO.write_Double_List(RMSIPs, out + name, 12);
			}

		RMSIP_Plot.create_RMSIP_Chart_ErrorBars(out, "Iterated Subspace Analysis with RMSIP Scores", "Subspace Dimension", "RMSIP", rmsipPlot2);

		write_FSSA_Iterated_Log();
		write_Random_FSSA_Log();
	}

	private void write_FSSA_Iterated_Log()
	{
		Job_log = new File(out + "JED_FSSA_Iterated_Log.txt");
		try
			{
				job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
			}
		catch (IOException e)
			{
				System.err.println("Problem writing the log file... Check permissions...  " + out + "JED_FSSA_Iterated_Log.txt");
				e.printStackTrace();
			}
		job_log_writer.write("Type of analysis: " + type_subset + "\n");
		job_log_writer.write("Comparison: " + type_comparison + "\n");
		job_log_writer.write("First Set of Eigenvectors:\n");
		job_log_writer.write("Rows: " + ROWS + "\n");
		job_log_writer.write("Cols: " + COLS + "\n");
		job_log_writer.write("Second Set of Eigenvectors:\n");
		job_log_writer.write("Rows: " + ROWS + "\n");
		job_log_writer.write("Cols: " + COLS + "\n");
		job_log_writer.write("Output Directory: " + out + "\n");
		if (verbose) job_log_writer.write("Principle Angle Spectra file written to: Iterated_PAs.txt\n");
		if (verbose) job_log_writer.write("RMSIPs file written to: Iterated_RMSIPs.txt\n\n");
		job_log_writer.write(String.format("%-12s%-12s%-12s%n%n", "SS_DIM", "RMSIP", "Z-SCORE"));
		int j = 1;
		for (double d : RMSIPs)
			{
				double zscore = Z_Scores.get(j - 1);
				job_log_writer.write(String.format("%-12s%-12s%-12s%n", j, nf6.format(d), nf2.format(zscore)));
				j++;
			}
		job_log_writer.write("\n" + "The PA spectra for the range of subspaces are:\n");
		PAs.print(job_log_writer, 3, 0);
		job_log_writer.write("\nAnalysis completed: " + date);
		job_log_writer.flush();
		job_log_writer.close();
	}

	/**
	 * Calls the Subspace Analysis class and runs the Random FSSA method: This compares two random subspaces having the same dim as the entered ones, using the metrics of RMSIP,
	 * PA, CO.
	 */
	public void get_Random_FSSA()
	{
		Subspace_Analysis2 rSSA = new Subspace_Analysis2(data1, data2);
		rSSA.get_random_SSA();

		avg_RMSIPs = rSSA.getAvg_RMSIP_Score_Matrix();
		avg_PAs = rSSA.getAvg_PA_Matrix();
		avg_COs = rSSA.getAvg_CO_Score_Matrix();
		RMSIP_std_devs = rSSA.getRMSIP_Std_Dev_Matrix();

		out = out_dir + type_comparison + File.separatorChar;
		exist = new File(out).exists();
		if (!exist) success = (new File(out)).mkdirs();

		if (verbose)
			{
				String name = "Avg_Random_RMSIPs.txt";
				Matrix_IO.write_Matrix(avg_RMSIPs.transpose(), out, name, 20, 12);
				name = "Avg_Random_PAs.txt";
				Matrix_IO.write_Matrix(avg_PAs.transpose(), out, name, 6, 0);
				name = "Avg_Random_COs.txt";
				Matrix_IO.write_Matrix(avg_COs.transpose(), out, name, 12, 6);
				name = "Random_RMSIP_Std_Devs.txt";
				Matrix_IO.write_Matrix(RMSIP_std_devs.transpose(), out, name, 20, 12);
			}
		write_Random_FSSA_Log();
	}

	private void write_Random_FSSA_Log()
	{
		Job_log = new File(out + "JED_Random_FSSA_Iterated_Log.txt");
		try
			{
				job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
			}
		catch (IOException e)
			{
				System.err.println("Problem writing the log file... Check permissions... " + out + "JED_Random_SSA_Log.txt");
				e.printStackTrace();
			}
		job_log_writer.write("Random analysis to establish baselines for the given vector space and dimension of subspaces.\n");
		job_log_writer.write("Output Directory: " + out + "\n");
		if (verbose) job_log_writer.write("Average RMSIPs file written to: Avg_Random_RMSIPs.txt\n");
		if (verbose) job_log_writer.write("Average PAs file written to: Avg_Random_PAs.txt\n");
		if (verbose) job_log_writer.write("Average COs file written to: Avg_Random_COs.txt\n");
		if (verbose) job_log_writer.write("RMSIP Std Devs file written to: Random_RMSIP_Std_Devs.txt\n");
		nf3.setMaximumFractionDigits(6);
		nf3.setMinimumFractionDigits(6);
		double[] rmsips = avg_RMSIPs.getRowPackedCopy();
		double[] std_devs = RMSIP_std_devs.getRowPackedCopy();
		int j = 1;
		job_log_writer.write("\nThe dimension of the vector space is  " + ROWS + "\n\n");
		job_log_writer.write(String.format("%-12s%-20s%-206s%n%n", "SS_DIM", "Mean RMSIP", "Standard_Deviation"));
		for (double r : rmsips)
			{
				double s = std_devs[j - 1];
				job_log_writer.write(String.format("%-12s%-20s%-20s%n", j, nf6.format(r), nf6.format(s)));
				j++;
			}
		job_log_writer.write("\n" + "The avg random PA spectra for the range of subspaces are: " + "\n");
		avg_PAs.transpose().print(job_log_writer, 3, 0);

		job_log_writer.write("\n" + "The avg random CO scores for the range of subspaces are: " + "\n");
		avg_COs.transpose().print(job_log_writer, 6, 3);

		job_log_writer.write("Analysis completed: " + date);
		job_log_writer.flush();
		job_log_writer.close();
	}
	// *************************************Setters & Getters**********************************************************

	public void setOut_dir(String out_dir)
	{
		this.out_dir = out_dir;
		exist = new File(out_dir).exists();
		if (!exist) success = (new File(out_dir)).mkdirs();
	}
}
