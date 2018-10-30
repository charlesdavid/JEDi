package jedi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;

import Jama.Matrix;

/**
 * JED class JED_Get_Subspace_Analysis: Performs Subspace Analysis on sets of eigenvectors with varying degree of outputs. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 *
 */

public class JEDi_Get_Subspace_Analysis
	{

		int num_of_jobs, ROWS, COLS, ROWS2, COLS2;
		String directory, out_dir, description, type, comparison, date = DateUtils.now();
		double RMSIP, max_angle;
		Matrix data1, data2, projections, projections_abs, CO_matrix_1_2, CO_matrix_2_1, PAs, avg_PAs, avg_RMSIPs, avg_COs, RMSIP_std_devs;
		ArrayList<Double> cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products, vectorial_sum_of_angles, RMSIPs;
		File Job_log;
		PrintWriter job_log_writer;
		NumberFormat nf;
		boolean exist, success;

		/**
		 * Constructor for performing Subspace Analysis Note: The matrices of eigenvectors must have the same dimensions for both rows and columns.
		 *
		 * @param dir             The working directory
		 * @param desc            The job description
		 * @param models_compared The type of PCA: COV or CORR (Q or R)
		 * @param evects1         The first set of eigenvectors for the subspace comparison
		 * @param evects2         The second set of eigenvectors for the subspace comparison
		 */

		public JEDi_Get_Subspace_Analysis(String dir, String desc, String models_compared, String pca_model, Matrix evects1, Matrix evects2)
			{

				this.directory = dir;
				this.description = desc;
				this.type = models_compared;
				this.comparison = pca_model;
				this.data1 = evects1;
				this.data2 = evects2;
				this.ROWS = data1.getRowDimension();
				this.COLS = data1.getColumnDimension();

				this.out_dir = directory + "JEDi_RESULTS_" + description + File.separatorChar + "SSA" + File.separatorChar + type + File.separatorChar + comparison
						+ File.separatorChar;
				exist = new File(out_dir).exists();
				if (!exist)
					success = (new File(out_dir)).mkdirs();

				nf = NumberFormat.getInstance();
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);

			}

		/**
		 * Calls the Subspace Analysis class and runs the SSA method: This compares the two subspaces directly using the metrics of RMSIP, PA, CO, CP, and VAS.
		 */
		public void get_SSA_JED()
			{
				Subspace_Analysis ssa = new Subspace_Analysis(data1, data2);
				ssa.get_SSA();

				RMSIP = ssa.getRMSIP();
				max_angle = ssa.getMax_angle();
				projections = ssa.getProjections();
				projections_abs = ssa.getProjections_abs();
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
				Matrix_IO.write_Matrix(CO_matrix_1_2, out_dir + "CO_1_2_dim_" + COLS + ".txt", 6, 3);
				Matrix_IO.write_Matrix(CO_matrix_2_1, out_dir + "CO_2_1_dim_" + COLS + ".txt", 6, 3);
				Matrix_IO.write_Matrix(projections, out_dir + "Projections_dim_" + COLS + ".txt", 6, 3);
				List_IO.write_Double_List(principle_angles_svd, out_dir + "PAs_dim_" + COLS + ".txt", 3);
				List_IO.write_Double_List(cosine_products, out_dir + "Cosine_Products_dim_" + COLS + ".txt", 3);
				List_IO.write_Double_List(vectorial_sum_of_angles, out_dir + "Vector_Sums_of_Angles_dim_" + COLS + ".txt", 3);
				write_SSA_Log();
			}

		private void write_SSA_Log()
			{
				Job_log = new File(out_dir + "JED_SSA_Log_dim_" + COLS + ".txt");
				try
					{
						job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
					} catch (IOException e)
					{
						System.err.println("Problem writing the log file: " + out_dir + "JED_SSA_Log_dim_" + COLS + ".txt");
						e.printStackTrace();
					}
				job_log_writer.write("Top_COV_Eigenvectors:\n");
				job_log_writer.write("Rows: " + ROWS + "\n");
				job_log_writer.write("Cols: " + COLS + "\n");
				job_log_writer.write("Top_CORR_Eigenvectors:\n");
				job_log_writer.write("Rows: " + ROWS + "\n");
				job_log_writer.write("Cols: " + COLS + "\n");
				job_log_writer.write("\nOutput Directory: " + out_dir + "\n");
				job_log_writer.write("Projections file written to: " + "Projections_dim_" + COLS + ".txt\n");
				job_log_writer.write("Cumulative overlaps 1 --> 2 file written to: " + "CO_1_2_dim_" + COLS + ".txt\n");
				job_log_writer.write("Cumulative overlaps 2 --> 1 file written to: " + "CO_2_1_dim_" + COLS + ".txt\n");
				job_log_writer.write("Principle Angles file written to: " + "PAs_dim_" + COLS + ".txt" + "\n");
				job_log_writer.write("Cosine Products file written to: " + "Cosine_Products_dim_" + COLS + ".txt\n");
				job_log_writer.write("Vectorial sums of angles file written to: " + "Vector_Sums_of_Angles_dim_" + COLS + ".txt\n");
				job_log_writer.write("\nThe Inner Products of each vector in subspace 1 with each vector in subspace 2 are:\n");
				projections.print(job_log_writer, 9, 3);
				job_log_writer.write("\nThe cumulative overlaps CO_" + COLS + " for each vector in subspace 1 with all the vectors in subspace 2 are:\n");
				int j = 1;
				for (double d : cumulative_overlaps_1_2)
					{
						job_log_writer.write(String.format("%-8s%-12s%-12s%n", "Vector ", j, nf.format(d)));
						j++;
					}
				job_log_writer.write("\nThe cumulative overlaps CO_" + COLS + " for each vector in subspace 2 with all the vectors in subspace 1 are:\n");
				j = 1;
				for (double d : cumulative_overlaps_2_1)
					{
						job_log_writer.write(String.format("%-8s%-12s%-12s%n", "Vector ", j, nf.format(d)));
						j++;
					}
				job_log_writer.write("\nThe RMSIP score is " + nf.format(RMSIP) + "\n");
				job_log_writer.write("\nThe principle angles (in degrees) are: " + "\n");
				int i = 1;
				nf.setMaximumFractionDigits(0);
				nf.setMinimumFractionDigits(0);
				for (double p : principle_angles_svd)
					{
						job_log_writer.write(String.format("%-6s%-12s%-12s%n", "PA", i, nf.format(p)));
						i++;
					}
				job_log_writer.write("\nThe cosine products (in degrees) are: " + "\n");
				i = 1;
				for (double c : cosine_products)
					{
						job_log_writer.write(String.format("%-6s%-12s%-12s%n", "CP", i, nf.format(Math.acos(c) * 180 / Math.PI)));
						i++;
					}
				job_log_writer.write("\nThe vectorial sums of angles (in degrees) are: " + "\n");
				i = 1;
				for (double p : vectorial_sum_of_angles)
					{
						job_log_writer.write(String.format("%-6s%-12s%-12s%n", "VS", i, nf.format(p)));
						i++;
					}
				job_log_writer.write("\nMaximum possible angle between two subspaces of this dimension is " + nf.format(max_angle) + " degrees\n\n");
				job_log_writer.write("Analysis completed: " + date);
				job_log_writer.flush();
				job_log_writer.close();
			}

		/**
		 * Calls the Subspace Analysis class and runs the FSSA Iterated method: This iteratively compares the two subspaces using the metrics of RMSIP, PA, CO, CP, and VAS.
		 * Comparisons are made for SS dims from 1 to the entered dim.
		 */
		public void get_FSSA_Iterated_JED()
			{
				Subspace_Analysis fssa = new Subspace_Analysis(data1, data2);
				fssa.get_fast_SSA_iterated();

				RMSIPs = fssa.getRMSIPs();
				PAs = fssa.getPAs();

				Matrix_IO.write_Matrix(PAs, out_dir + "Iterated_PAs.txt", 6, 3);
				List_IO.write_Double_List(RMSIPs, out_dir + "Iterated_RMSIPs.txt", 3);
				write_FSSA_Iterated_Log();
			}

		private void write_FSSA_Iterated_Log()
			{
				Job_log = new File(out_dir + "JED_FSSA_Iterated_Log.txt");
				try
					{
						job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
					} catch (IOException e)
					{
						System.err.println("Problem writing the log file: " + out_dir + "JED_FSSA_Iterated_Log.txt");
						e.printStackTrace();
					}
				job_log_writer.write("Output Directory: " + out_dir + "\n");
				job_log_writer.write("Principle Angle Spectra file written to: Iterated_PAs.txt\n");
				job_log_writer.write("RMSIPs file written to: Iterated_RMSIPs.txt\n");
				job_log_writer.write("\n\nRMSIPs:\n");
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
				int j = 1;
				for (double d : RMSIPs)
					{
						job_log_writer.write(String.format("%-6s%-12s%-12s%n", "Dim", j, nf.format(d)));
						j++;
					}
				job_log_writer.write("\n" + "The PA spectra for the range of subspaces are:\n");
				PAs.print(job_log_writer, 3, 0);
				job_log_writer.write("\nAnalysis completed: " + date);
				job_log_writer.flush();
				job_log_writer.close();
			}

		/**
		 * Calls the Subspace Analysis class and runs the Random FSSA method: This compares two random subspaces having the same dim as the entered ones, using the metrics of
		 * RMSIP, PA, CO.
		 */
		public void get_Random_FSSA_JED()
			{
				Subspace_Analysis rssa = new Subspace_Analysis(data1, data2);
				rssa.get_random_SSA();

				avg_RMSIPs = rssa.getAvg_RMSIP_Score_Matrix();
				avg_PAs = rssa.getAvg_PA_Matrix();
				avg_COs = rssa.getAvg_CO_Score_Matrix();
				RMSIP_std_devs = rssa.getRMSIP_Std_Dev_Matrix();

				String name = "Avg_Random_RMSIPs.txt";
				Matrix_IO.write_Matrix(avg_RMSIPs.transpose(), out_dir, name, 6, 3);
				name = "Avg_Random_PAs.txt";
				Matrix_IO.write_Matrix(avg_PAs.transpose(), out_dir, name, 6, 0);
				name = "Avg_Random_COs.txt";
				Matrix_IO.write_Matrix(avg_COs.transpose(), out_dir, name, 6, 3);
				name = "Random_RMSIP_Std_Devs.txt";
				Matrix_IO.write_Matrix(RMSIP_std_devs.transpose(), out_dir, name, 6, 3);
				write_Random_FSSA_Log();
			}

		private void write_Random_FSSA_Log()
			{
				Job_log = new File(out_dir + "JED_Random_FSSA_Iterated_Log.txt");
				try
					{
						job_log_writer = new PrintWriter(new BufferedWriter(new FileWriter(Job_log)));
					} catch (IOException e)
					{
						System.err.println("Problem writing the log file: " + out_dir + "JED_Random_SSA_Log.txt");
						e.printStackTrace();
					}
				job_log_writer.write("Output Directory: " + out_dir + "\n");
				job_log_writer.write("Average RMSIPs file written to: Avg_Random_RMSIPs.txt\n");
				job_log_writer.write("Average PAs file written to: Avg_Random_PAs.txt\n");
				job_log_writer.write("Average COs file written to: Avg_Random_COs.txt\n");
				job_log_writer.write("RMSIP Std Devs file written to: Random_RMSIP_Std_Devs.txt\n");
				nf.setMaximumFractionDigits(3);
				nf.setMinimumFractionDigits(3);
				double[] rmsips = avg_RMSIPs.getRowPackedCopy();
				double[] std_devs = RMSIP_std_devs.getRowPackedCopy();
				int j = 1;
				job_log_writer.write("\nThe dimension of the vector space is  " + ROWS + "\n\n");
				job_log_writer.write(String.format("%-12s%-16s%-16s%n%n", "SS_DIM", "RMSIP", "Std_Dev"));
				for (double r : rmsips)
					{
						double s = std_devs[j - 1];
						job_log_writer.write(String.format("%-12s%-16s%-16s%n", j, nf.format(r), nf.format(s)));
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
			}
	}
