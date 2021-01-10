package supportIO;

import java.io.File;
import java.util.Hashtable;

public class Input_Parameters_ALIGN
{
	public static boolean exist, success, PROCESS_OUTLIERS, COMPRESS;
	public static int REFERENCE_COL;
	public static double MAD_Score, Z_Score;
	public static String OUT_DIR, DESCRIPTION, PDB_COORDS, COMPRESS_METHOD;
	public static Hashtable<String, String> parameters;

	// ****************************************** CONSTRUCTOR *********************************************************************************** //

	public Input_Parameters_ALIGN(Hashtable<String, String> parameters_read)
	{
		parameters = parameters_read;

		OUT_DIR = parameters.get("OUT_DIR");

		if (!(OUT_DIR.endsWith(File.separator)))
			{
				System.err.println("Expected the directory to end with " + File.separator + ", but got: " + OUT_DIR);
				System.err.println("Attempting to fix...");
				OUT_DIR = OUT_DIR + File.separator;
			}

		boolean dir = new File(OUT_DIR).isDirectory();
		if (!dir)
			{
				System.err.println("Sorry, but the entered directory is not recognized as a proper directory: " + OUT_DIR);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}

		exist = new File(OUT_DIR).exists();
		if (!exist)
			{
				success = (new File(OUT_DIR)).mkdirs();
				if (!success)
					{
						System.err.println("Sorry, unable to create the output directory. Please check synatx and file permissions: " + OUT_DIR);
						System.exit(0);
					}
			}

		PDB_COORDS = parameters.get("PDB_COORDS");
		exist = new File(PDB_COORDS).exists();
		if (!exist)

			{
				System.err.println("Sorry, the entered PDB Coordinates File can not be found: " + PDB_COORDS);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}

		DESCRIPTION = parameters.get("DESCRIPTION");
		REFERENCE_COL = Integer.valueOf(parameters.get("REFERENCE_COL"));
		MAD_Score = Double.valueOf(parameters.get("MAD_Score"));
		Z_Score = Double.valueOf(parameters.get("Z_Score"));

		if (MAD_Score > 0 && Z_Score > 0)
			{
				System.err.println("Sorry, but you can NOT use both MAD scores and Z_Scores simultaneously.");
				System.err.println("Setting MAD_Score cutoff to zero and and using Z_Score Cutoff.");
				MAD_Score = 0;
			}

		if (parameters.get("PROCESS_OUTLIERS").equals("true")) PROCESS_OUTLIERS = true;
		if (parameters.get("COMPRESS").equals("true")) COMPRESS = true;

		if (COMPRESS) COMPRESS_METHOD = parameters.get("COMPRESS_METHOD"); // Either 'bz2' or 'gz' ONLY!
	}
}