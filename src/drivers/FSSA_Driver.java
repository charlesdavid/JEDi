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
import support.Subspace_Analysis;
import supportIO.Matrix_IO;

public class FSSA_Driver
{

	static int num_of_jobs, ROWS1, COLS1, ROWS2, COLS2;
	static double RMSIP;
	static String directory1, directory2, name1, name2, out_dir, description, batch_description, date = supportIO.DateUtils.now();
	static Matrix matrix1, matrix2;
	static ArrayList<Double> cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products, vectorial_sum_of_angles, RMSIPs;
	static File file_1, file_2, Job_log, Batch_log;
	static BufferedReader reader;
	static BufferedWriter Batch_log_Writer;
	static PrintWriter Job_log_writer;
	static NumberFormat nf;

	public static void main(String[] args) throws Exception
	{

		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setMinimumFractionDigits(1);
		reader = new BufferedReader(new FileReader("SSA.txt"));
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
				name1 = stoken.nextToken();
				nextline = reader.readLine();
				stoken = new StringTokenizer(nextline);
				directory2 = stoken.nextToken();
				name2 = stoken.nextToken();
				file_1 = new File(directory1 + name1);
				file_2 = new File(directory2 + name2);
				matrix1 = Matrix_IO.read_Matrix(directory1, name1);
				ROWS1 = matrix1.getRowDimension();
				COLS1 = matrix1.getColumnDimension();
				matrix2 = Matrix_IO.read_Matrix(directory2, name2);
				ROWS2 = matrix2.getRowDimension();
				COLS2 = matrix2.getColumnDimension();

				cumulative_overlaps_1_2 = new ArrayList<Double>();
				cumulative_overlaps_2_1 = new ArrayList<Double>();
				principle_angles_svd = new ArrayList<Double>();
				cosine_products = new ArrayList<Double>();
				vectorial_sum_of_angles = new ArrayList<Double>();

				System.out.println("Matrix 1: " + directory1 + name1 + "\n");
				System.out.println("Rows: " + ROWS1 + "\n");
				System.out.println("Cols: " + COLS1 + "\n");
				System.out.println("Matrix 2: " + directory2 + name2 + "\n");
				System.out.println("Rows: " + ROWS2 + "\n");
				System.out.println("Cols: " + COLS2 + "\n");

				Subspace_Analysis ssa = new Subspace_Analysis(matrix1, matrix2);
				ssa.get_fast_SSA();
				RMSIP = ssa.getRMSIP();
				principle_angles_svd = ssa.getPrinciple_angles_svd();

				Job_log = new File(out_dir + "FSSA_Job_log_" + description + "_dim_" + matrix1.getColumnDimension() + ".txt");
				Job_log_writer = new PrintWriter(new FileWriter(Job_log));
				Job_log_writer.write("Matrix 1: " + directory1 + name1 + "\n");
				Job_log_writer.write("Rows: " + matrix1.getRowDimension() + "\n");
				Job_log_writer.write("Cols: " + matrix1.getColumnDimension() + "\n");
				Job_log_writer.write("Matrix 2: " + directory2 + name2 + "\n");
				Job_log_writer.write("Rows: " + matrix2.getRowDimension() + "\n");
				Job_log_writer.write("Cols: " + matrix2.getColumnDimension() + "\n" + "\n");
				Job_log_writer.write("Principle Angles file written to: " + "\n" + out_dir + description + "_PA_dim_" + matrix1.getColumnDimension() + ".txt" + "\n" + "\n");
				Batch_log_Writer.write(description + "\t" + "SS_dim " + "\t" + COLS1 + "\t" + "VS_dim" + "\t" + ROWS1 + "\t" + "RMSIP" + "\t" + nf.format(RMSIP) + "\t"
						+ "Principle_Angles" + "\t");
				Job_log_writer.write("\n" + "The RMSIP for the two subspaces is " + nf.format(RMSIP) + "\n");
				Job_log_writer.write("\n" + "The principle angles between these spaces (in degrees) are: " + "\n");
				int i = 1;
				nf.setMaximumIntegerDigits(3);
				nf.setMaximumFractionDigits(2);
				for (double p : principle_angles_svd)
					{
						Job_log_writer.write("Angle " + i + "\t" + nf.format(p) + "\n");
						Batch_log_Writer.write(nf.format(p) + "\t");
						i++;
					}
				Job_log_writer.write("Analysis completed: " + date);
				Job_log_writer.close();
				Batch_log_Writer.flush();
				principle_angles_svd.clear();
			}
		reader.close();
		Batch_log_Writer.close();
	}

}
