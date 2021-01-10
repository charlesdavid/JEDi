package supportIO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.List;

public class Write_DP_Stats
{
	static NumberFormat nf3, nf12;
	static RoundingMode rm;
	static boolean exist, success;

	public static void write_stats(String path, double[] means, double[] std_devs, List<String> IDs1, List<String> IDs2, List<Integer> Atoms1, List<Integer> Atoms2)
	{
		rm = RoundingMode.HALF_UP;
		nf3 = NumberFormat.getInstance();
		nf3.setRoundingMode(rm);
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		int number_of_pairs = IDs1.size();
		try
			{
				File stats = new File(path);
				BufferedWriter writer = new BufferedWriter(new FileWriter(stats));
				writer.write("MEANs and STANDARD DEVIATIONS for the Atom Pair Distances: " + "\n");
				writer.write(String.format("%-20s%-20s%-20s%-20s%n", "Atom1", "Atom2", "Mean", "Std_Dev"));
				for (int i = 0; i < number_of_pairs; i++)
					{
							{
								writer.write(String.format("%-20s%-20s%-20s%-20s%n", IDs1.get(i) + Atoms1.get(i), IDs2.get(i) + Atoms2.get(i), nf3.format(means[i]),
										nf3.format(std_devs[i])));
							}
					}
				writer.close();
			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not write the file: " + path);
				System.err.println("Please check file and/or directory permissions.");
				io.printStackTrace();
			}
	}
}
