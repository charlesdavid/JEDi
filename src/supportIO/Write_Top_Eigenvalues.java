package supportIO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;

public class Write_Top_Eigenvalues
{
	static NumberFormat nf6, nf12;
	static RoundingMode rm;

	public static void write_evals(String path, List<Double> evals, double trace)
	{
		rm = RoundingMode.HALF_UP;

		nf6 = NumberFormat.getInstance();
		nf6.setRoundingMode(rm);
		nf6.setMaximumFractionDigits(6);
		nf6.setMinimumFractionDigits(6);

		nf12 = NumberFormat.getInstance();
		nf12.setRoundingMode(rm);
		nf12.setMaximumFractionDigits(12);
		nf12.setMinimumFractionDigits(12);

		try
			{
				File out = new File(path);
				BufferedWriter writer = new BufferedWriter(new FileWriter(out));
				writer.write(String.format("%-30s%-30s%-30s%n", "Eigenvalue", "% Variance", "Cumulative Variance"));

				double cumulative_variance = 0;

				for (int i = 0; i < evals.size(); i++)
					{
						double val = evals.get(i);
						double normed_val = (val / trace);
						cumulative_variance += normed_val;
						writer.write(String.format("%-30s%-30s%-30s%n", nf12.format(val), nf6.format(normed_val), nf6.format(cumulative_variance)));
					}
				writer.close();
			}
		catch (IOException io)
			{
				System.err.println("Could not write to the file: " + path);
				io.printStackTrace();
			}
	}
}
