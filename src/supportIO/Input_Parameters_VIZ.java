package supportIO;

import java.io.File;
import java.util.Hashtable;

public class Input_Parameters_VIZ
{
	public static boolean exist, success, doESSENTIAL, doINDIVIDUAL;
	public static int MODE_START, MODE_END, numberOfModes, numberOfFrames, numberOfCycles, numberOfFramesEssential, numberOfCyclesEssential;
	public static double FLOOR, mode_amplitude, threshold_low, threshold_high;
	public static String OUT_DIR, DESCRIPTION, REFERENCE_PDB, EVALS, EVECTS, modelPCA;
	public static Hashtable<String, String> parameters;

	// ****************************************** CONSTRUCTOR *********************************************************************************** //

	public Input_Parameters_VIZ(Hashtable<String, String> parameters_read)
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

		DESCRIPTION = parameters.get("DESCRIPTION");
		REFERENCE_PDB = parameters.get("REFERENCE_PDB");

		exist = new File(REFERENCE_PDB).exists();
		if (!exist)

			{
				System.err.println("Sorry, the entered PDB Reference File can not be found: " + REFERENCE_PDB);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}

		EVALS = parameters.get("EIGENVALUES");
		exist = new File(EVALS).exists();

		if (!exist)

			{
				System.err.println("Sorry, the entered PDB Reference File can not be found: " + EVALS);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}

		EVECTS = parameters.get("EIGENVECTORS");
		exist = new File(EVECTS).exists();

		if (!exist)

			{
				System.err.println("Sorry, the entered PDB Reference File can not be found: " + EVECTS);
				System.err.println("Terminating program execution.");
				System.exit(0);
			}

		FLOOR = Double.valueOf(parameters.get("FLOOR"));
		mode_amplitude = Double.valueOf(parameters.get("mode_amplitude"));
		threshold_low = Double.valueOf(parameters.get("threshold_low"));
		threshold_high = Double.valueOf(parameters.get("threshold_high"));

		numberOfFrames = Integer.valueOf(parameters.get("numberOfFrames"));
		numberOfCycles = Integer.valueOf(parameters.get("numberOfCycles"));
		numberOfFramesEssential = Integer.valueOf(parameters.get("numberOfFramesEssential"));
		numberOfCyclesEssential = Integer.valueOf(parameters.get("numberOfCyclesEssential"));

		MODE_START = Integer.valueOf(parameters.get("MODE_START"));
		MODE_END = Integer.valueOf(parameters.get("MODE_END"));
		numberOfModes = (MODE_END - MODE_START + 1);
		modelPCA = parameters.get("modelPCA");

		doESSENTIAL = Boolean.valueOf(parameters.get("doESSENTIAL"));
		doINDIVIDUAL = Boolean.valueOf(parameters.get("doINDIVIDUAL"));
	}
}