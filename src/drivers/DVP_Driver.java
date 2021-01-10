package drivers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;

import Jama.Matrix;
import support.Projector;
import supportIO.List_IO;
import supportIO.Matrix_IO;
import supportPlot.PC_Plot;

public class DVP_Driver
{

	static int num_of_jobs, number_of_modes = 2, ROWS1, COLS1, ROWS2, COLS2;
	static double RMSIP;
	static String directory1, directory2, directory3, name1, name2, name3, out_dir, description, batch_description, date = supportIO.DateUtils.now();
	static Matrix delta_vectors, eigenvectors;
	static Matrix projections, normed_projections, weighted_projections, weighted_normed_projections;
	static ArrayList<Double> eigenvalues;
	static File Job_log, Batch_log;
	static BufferedReader reader;
	static BufferedWriter Batch_log_Writer;
	static PrintWriter Job_log_writer;
	static NumberFormat nf;

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception
	{

		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setMinimumFractionDigits(1);
		reader = new BufferedReader(new FileReader("DVP.txt"));
		num_of_jobs = Integer.parseInt(reader.readLine());
		out_dir = reader.readLine();
		batch_description = reader.readLine();
		Batch_log = new File(out_dir + "data_log_" + batch_description + ".txt");
		Batch_log_Writer = new BufferedWriter(new FileWriter(Batch_log));

		for (int index = 0; index < (num_of_jobs); index++)
		{
			description = reader.readLine();

			String nextline = reader.readLine();
			StringTokenizer stoken = new StringTokenizer(nextline);
			directory1 = stoken.nextToken();
			name1 = stoken.nextToken(); // DELTA VECTORS

			nextline = reader.readLine();
			stoken = new StringTokenizer(nextline);
			directory2 = stoken.nextToken();
			name2 = stoken.nextToken(); // EIGENVECTORS

			nextline = reader.readLine();
			stoken = new StringTokenizer(nextline);
			directory3 = stoken.nextToken();
			name3 = stoken.nextToken(); // EIGENVALUES

			delta_vectors = Matrix_IO.read_Matrix(directory1, name1);
			ROWS1 = delta_vectors.getRowDimension();
			COLS1 = delta_vectors.getColumnDimension();

			eigenvectors = Matrix_IO.read_Matrix(directory2, name2);
			ROWS2 = eigenvectors.getRowDimension();
			COLS2 = eigenvectors.getColumnDimension();

			eigenvalues = (ArrayList<Double>) List_IO.read_List_From_File(directory3 + name3, "Double");

			System.out.println("Matrix 1: " + directory1 + name1 + "\n");
			System.out.println("Rows: " + ROWS1 + "\n");
			System.out.println("Cols: " + COLS1 + "\n");
			System.out.println("Matrix 2: " + directory2 + name2 + "\n");
			System.out.println("Rows: " + ROWS2 + "\n");
			System.out.println("Cols: " + COLS2 + "\n");

			projections = new Matrix(COLS1, COLS2);
			normed_projections = new Matrix(COLS1, COLS2);
			weighted_projections = new Matrix(COLS1, COLS2);
			weighted_normed_projections = new Matrix(COLS1, COLS2);

			for (int mode = 0; mode < COLS2; mode++)
			{
				double val = eigenvalues.get(mode);
				if (val < 0) val = 0;
				double weight = Math.sqrt(val); // Weight has units of Angstroms

				Matrix data1 = eigenvectors.getMatrix(0, ROWS1 - 1, mode, mode);
				Matrix vector1 = Projector.get_Normed_arrayF(data1);

				for (int conf = 0; conf < COLS1; conf++)
				{
					Matrix data2 = delta_vectors.getMatrix(0, ROWS1 - 1, conf, conf);
					Matrix vector2 = Projector.get_Normed_arrayF(data2);

					double dp = Projector.get_InnerProduct(data1, data2);
					double normed_dp = Projector.get_InnerProduct(vector1, vector2);
					double w_dp = (weight * dp);
					double weighted_normed_dp = (weight * normed_dp);

					projections.set(conf, mode, dp);
					normed_projections.set(conf, mode, normed_dp);
					weighted_projections.set(conf, mode, w_dp);
					weighted_normed_projections.set(conf, mode, weighted_normed_dp);
				}
			}
			String path = out_dir + COLS2 + "_DVPs.txt";
			Matrix_IO.write_Matrix_adaptive_spacing(projections, path);
			path = out_dir + COLS2 + "_normed_DVPs.txt";
			Matrix_IO.write_Matrix_adaptive_spacing(normed_projections, path);
			path = out_dir + COLS2 + "_weighted_DVPs.txt";
			Matrix_IO.write_Matrix_adaptive_spacing(weighted_projections, path);
			path = out_dir + COLS2 + "_weighted_normed_DVPs.txt";
			Matrix_IO.write_Matrix_adaptive_spacing(weighted_normed_projections, path);

			PC_Plot.createChart(out_dir, description + "_Projections", projections);
			PC_Plot.createChart(out_dir, description + "_Normed_Projections", normed_projections);
			PC_Plot.createChart(out_dir, description + "_Weighted_Projections", weighted_projections);
			PC_Plot.createChart(out_dir, description + "_Weighted_Normed_Projections", weighted_normed_projections);

			Job_log = new File(out_dir + "DVP_Job_log_" + description + "_dim_" + COLS1 + ".txt");
			Job_log_writer = new PrintWriter(new FileWriter(Job_log));
			Job_log_writer.write("Matrix 1: " + directory1 + name1 + "\n");
			Job_log_writer.write("Rows: " + delta_vectors.getRowDimension() + "\n");
			Job_log_writer.write("Cols: " + delta_vectors.getColumnDimension() + "\n");
			Job_log_writer.write("Matrix 2: " + directory2 + name2 + "\n");
			Job_log_writer.write("Rows: " + eigenvectors.getRowDimension() + "\n");
			Job_log_writer.write("Cols: " + eigenvectors.getColumnDimension() + "\n" + "\n");
			Batch_log_Writer.write("Job " + index + description);
			Job_log_writer.write("Analysis completed: " + date);
			Job_log_writer.flush();
			Job_log_writer.close();
			Batch_log_Writer.flush();
		}
		reader.close();
		Batch_log_Writer.close();
	}
}
