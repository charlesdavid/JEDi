package supportIO;

import java.io.File;
import java.util.Hashtable;

public class Input_Parameters
{
	public static boolean exist, success;
	public static boolean doPREPROCESS, doREAD_ARCHIVE, doParityCheck, doOutputCoordinates, verbose;
	public static boolean doAA, doBB, doHA, doCA;
	public static boolean doDownSample, doFrameSelect, doREDUCE, doFES, doEssentialViz, doModeViz, doKPCA, doOutlierProcessing, doSPARSIFY, do_StatThresholds;
	public static boolean doRESIDUE_INDIVIDUAL, doRESIDUE_PAIRS, doHIERARCHICAL_AA, doHIERARCHICAL_HA, doATOM_LIST, doDIST_PAIRS, doCORR, doPCORR;
	public static boolean do_cartesian_AA, do_cartesian_BB, do_cartesian_HA, do_cartesian_CA, do_hierarchical_AA, do_hierarchical_HA, do_dist_pairs, do_residue_individual,
			do_residue_pairs, do_atom_list;
	public static boolean Linear, Degree_2_Poly, Degree_3_Poly, Degree_4_Poly, XY_Poly, Euclidean, Mahalanobis, Gaussian, Sigmoid, Log, Circular, Cauchy, MI, MI_KDE;


	public static int MODES_RESIDUE_INDIVIDUAL, MODES_EIGEN_RESIDUE_PAIRS, MODES_EIGEN_RESIDUE_AA, MODES_EIGEN_RESIDUE_HA, MODES_ALL_ATOM, MODES_HEAVY_ATOM, MODES_BACKBONE,
			MODES_ALPHA_CARBON, MODES_HIERARCHICAL_AA, MODES_HIERARCHICAL_HA, MODES_ATOMS_LIST, MODES_DISTANCE_PAIRS;

	public static int DOWNSAMPLE, FRAME_START, FRAME_END, MULTIPLIER;
	public static int MODES_VIZ, MAX_KERNEL_FRAMES, NUMBER_PCs_INPUT, numberModeCycles, numberModeFrames, numberModeComponents, numberEssentialCycles, numberEssentialFrames;

	public static double Z_SCORE_CUTOFF, MAD_SCORE_CUTOFF, KDE_RESOLUTION, KDE_CELL, KDE_MARGIN, KERNEL_SHRINKAGE, KPCA_SIGMA, KPCA_SLOPE;
	public static double VIZ_MODE_SCALE_FACTOR, LOG_FLOOR, FLOOR, NOISE_LEVEL, THRESHOLD_COV, THRESHOLD_CORR, THRESHOLD_PCORR, THRESHOLD_RP_DIFF;
	public static double VARIANCE_THRESHOLD, SKEW_THRESHOLD, KURTOSIS_THRESHOLD;

	public static String DIRECTORY, DESCRIPTION, REFERENCE_PDB, READ_PDBS_FILTER_STRING, ORIGINAL_PDB_COORDS, OUT_DIR, ARCHIVE_NAME;

	public static String RESIDUE_LIST_ALL_ATOM, RESIDUE_LIST_HEAVY_ATOM, RESIDUE_LIST_BACKBONE, RESIDUE_LIST_ALPHA_CARBON, RESIDUE_LIST_INDIVIDUAL, RESIDUE_LIST_PAIRS,
			RESIDUE_LIST_HIERARCHICAL_AA, RESIDUE_LIST_HIERARCHICAL_HA, ATOMS_LIST, ATOM_PAIRS_LIST;

	public static Hashtable<String, String> parameters;

	// ****************************************** CONSTRUCTOR *********************************************************************************** //

	public Input_Parameters(Hashtable<String, String> parameters_read)
	{
		parameters = parameters_read;

		DIRECTORY = parameters.get("DIRECTORY");
		DESCRIPTION = parameters.get("DESCRIPTION");
		REFERENCE_PDB = parameters.get("REFERENCE_PDB");

		if (parameters.get("doPREPROCESS").equals("true")) doPREPROCESS = true;

		if (doPREPROCESS)
			{
				if (parameters.get("doParityCheck").equals("true")) doParityCheck = true;

				if (parameters.get("doREAD_ARCHIVE").equals("true"))
					{
						doREAD_ARCHIVE = true;
						ARCHIVE_NAME = parameters.get("ARCHIVE_NAME");
					}

				do_cartesian_AA = false;
				do_cartesian_BB = false;
				do_cartesian_HA = false;
				do_cartesian_CA = false;
				do_residue_individual = false;
				do_residue_pairs = false;
				do_hierarchical_AA = false;
				do_hierarchical_HA = false;
				do_atom_list = false;
				do_dist_pairs = false;
				doEssentialViz = false;
				doModeViz = false;

				READ_PDBS_FILTER_STRING = parameters.get("READ_PDBS_FILTER_STRING");

				if (parameters.get("do_StatThresholds").equals("true")) do_StatThresholds = true;
				VARIANCE_THRESHOLD = Double.valueOf(parameters.get("VARIANCE_THRESHOLD"));
				SKEW_THRESHOLD = Double.valueOf(parameters.get("SKEW_THRESHOLD"));
				KURTOSIS_THRESHOLD = Double.valueOf(parameters.get("KURTOSIS_THRESHOLD"));

				if (parameters.get("doOutputCoordinates").equals("true")) doOutputCoordinates = true;
				if (parameters.get("verbose").equals("true")) verbose = true;

				if (!(DIRECTORY.endsWith(File.separator)))
					{
						System.err.println("Expected the directory to end with " + File.separator + ", but got: " + DIRECTORY);
						System.err.println("Attempting to fix...");
						DIRECTORY = DIRECTORY + File.separator;
					}
				boolean dir = new File(DIRECTORY).isDirectory();
				if (!dir)
					{
						System.err.println("Sorry, but the entered directory is not recognized as a proper directory: " + DIRECTORY);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				boolean check = new File(DIRECTORY + REFERENCE_PDB).exists();
				if (!check)

					{
						System.err.println("Sorry, the entered PDB Reference File can not be found: " + DIRECTORY + REFERENCE_PDB);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
			}

		if (!doPREPROCESS)
			{
				if (parameters.get("doAA").equals("true")) doAA = true;
				if (parameters.get("doBB").equals("true")) doBB = true;
				if (parameters.get("doCA").equals("true")) doCA = true;
				if (parameters.get("doHA").equals("true")) doHA = true;

				if (parameters.get("doRESIDUE_INDIVIDUAL").equals("true")) doRESIDUE_INDIVIDUAL = true;
				if (parameters.get("doRESIDUE_PAIRS").equals("true")) doRESIDUE_PAIRS = true;
				if (parameters.get("doHIERARCHICAL_AA").equals("true")) doHIERARCHICAL_AA = true;
				if (parameters.get("doHIERARCHICAL_HA").equals("true")) doHIERARCHICAL_HA = true;
				if (parameters.get("doATOM_LIST").equals("true")) doATOM_LIST = true;
				if (parameters.get("doDIST_PAIRS").equals("true")) doDIST_PAIRS = true;

				if (parameters.get("doCORR").equals("true")) doCORR = true;
				if (parameters.get("doPCORR").equals("true")) doPCORR = true;

				if (parameters.get("do_StatThresholds").equals("true")) do_StatThresholds = true;
				VARIANCE_THRESHOLD = Double.valueOf(parameters.get("VARIANCE_THRESHOLD"));
				SKEW_THRESHOLD = Double.valueOf(parameters.get("SKEW_THRESHOLD"));
				KURTOSIS_THRESHOLD = Double.valueOf(parameters.get("KURTOSIS_THRESHOLD"));

				if (parameters.get("doSPARSIFY").equals("true")) doSPARSIFY = true;
				THRESHOLD_CORR = Double.valueOf(parameters.get("THRESHOLD_CORR"));
				THRESHOLD_PCORR = Double.valueOf(parameters.get("THRESHOLD_PCORR"));
				THRESHOLD_RP_DIFF = Double.valueOf(parameters.get("THRESHOLD_RP_DIFF"));

				if (parameters.get("doREDUCE").equals("true")) doREDUCE = true;
				if (parameters.get("doFES").equals("true")) doFES = true;
				if (parameters.get("doKPCA").equals("true")) doKPCA = true;

				MAX_KERNEL_FRAMES = Integer.valueOf(parameters.get("MAX_KERNEL_FRAMES"));
				NUMBER_PCs_INPUT = Integer.valueOf(parameters.get("NUMBER_PCs_INPUT"));
				KERNEL_SHRINKAGE = Double.valueOf(parameters.get("KERNEL_SHRINKAGE"));
				KDE_RESOLUTION = Double.valueOf(parameters.get("KDE_RESOLUTION"));
				KDE_CELL = Double.valueOf(parameters.get("KDE_CELL"));
				KDE_MARGIN = Double.valueOf(parameters.get("KDE_MARGIN"));
				KPCA_SIGMA = Double.valueOf(parameters.get("KPCA_SIGMA"));
				KPCA_SLOPE = Double.valueOf(parameters.get("KPCA_SLOPE"));
				MULTIPLIER = Integer.valueOf(parameters.get("MULTIPLIER"));

				if (parameters.get("Linear").equals("true")) Linear = true;
				if (parameters.get("Degree_2_Poly").equals("true")) Degree_2_Poly = true;
				if (parameters.get("Degree_3_Poly").equals("true")) Degree_3_Poly = true;
				if (parameters.get("Degree_4_Poly").equals("true")) Degree_4_Poly = true;
				if (parameters.get("XY_Poly").equals("true")) XY_Poly = true;
				if (parameters.get("Euclidean").equals("true")) Euclidean = true;
				if (parameters.get("Mahalanobis").equals("true")) Mahalanobis = true;
				if (parameters.get("Gaussian").equals("true")) Gaussian = true;
				if (parameters.get("Sigmoid").equals("true")) Sigmoid = true;
				if (parameters.get("Log").equals("true")) Log = true;
				if (parameters.get("Circular").equals("true")) Circular = true;
				if (parameters.get("Cauchy").equals("true")) Cauchy = true;
				if (parameters.get("MI").equals("true")) MI = true;
				if (parameters.get("MI_KDE").equals("true")) MI_KDE = true;

				if (parameters.get("doOutlierProcessing").equals("true")) doOutlierProcessing = true;

				Z_SCORE_CUTOFF = Double.valueOf(parameters.get("Z_SCORE_CUTOFF"));
				MAD_SCORE_CUTOFF = Double.valueOf(parameters.get("MAD_SCORE_CUTOFF"));

				if (parameters.get("doModeViz").equals("true")) doModeViz = true;

				numberModeCycles = Integer.valueOf(parameters.get("numberModeCycles"));
				numberModeFrames = Integer.valueOf(parameters.get("numberModeFrames"));

				if (parameters.get("doEssentialViz").equals("true")) doEssentialViz = true;

				numberEssentialCycles = Integer.valueOf(parameters.get("numberEssentialCycles"));
				numberEssentialFrames = Integer.valueOf(parameters.get("numberEssentialFrames"));
				numberModeComponents = Integer.valueOf(parameters.get("numberModeComponents"));

				MODES_VIZ = Integer.valueOf(parameters.get("MODES_VIZ"));
				VIZ_MODE_SCALE_FACTOR = Double.valueOf(parameters.get("VIZ_MODE_SCALE_FACTOR"));
				LOG_FLOOR = Double.valueOf(parameters.get("LOG_FLOOR"));

				if (parameters.get("doDownSample").equals("true")) doDownSample = true;

				DOWNSAMPLE = Integer.valueOf(parameters.get("DOWNSAMPLE"));

				if (parameters.get("doFrameSelect").equals("true")) doFrameSelect = true;

				FRAME_START = Integer.valueOf(parameters.get("FRAME_START"));
				FRAME_END = Integer.valueOf(parameters.get("FRAME_END"));

				RESIDUE_LIST_ALL_ATOM = parameters.get("RESIDUE_LIST_ALL_ATOM");
				RESIDUE_LIST_HEAVY_ATOM = parameters.get("RESIDUE_LIST_HEAVY_ATOM");
				RESIDUE_LIST_BACKBONE = parameters.get("RESIDUE_LIST_BACKBONE");
				RESIDUE_LIST_ALPHA_CARBON = parameters.get("RESIDUE_LIST_ALPHA_CARBON");

				RESIDUE_LIST_INDIVIDUAL = parameters.get("RESIDUE_LIST_INDIVIDUAL");
				RESIDUE_LIST_PAIRS = parameters.get("RESIDUE_LIST_PAIRS");
				RESIDUE_LIST_HIERARCHICAL_AA = parameters.get("RESIDUE_LIST_HIERARCHICAL_AA");
				RESIDUE_LIST_HIERARCHICAL_HA = parameters.get("RESIDUE_LIST_HIERARCHICAL_HA");
				ATOMS_LIST = parameters.get("ATOMS_LIST");
				ATOM_PAIRS_LIST = parameters.get("ATOM_PAIRS_LIST");

				ORIGINAL_PDB_COORDS = parameters.get("ORIGINAL_PDB_COORDS");

				FLOOR = Double.valueOf(parameters.get("FLOOR"));
				NOISE_LEVEL = Double.valueOf(parameters.get("NOISE_LEVEL"));
				THRESHOLD_COV = Double.valueOf(parameters.get("THRESHOLD_COV"));

				MODES_RESIDUE_INDIVIDUAL = Integer.valueOf(parameters.get("MODES_RESIDUE_INDIVIDUAL"));
				MODES_EIGEN_RESIDUE_PAIRS = Integer.valueOf(parameters.get("MODES_EIGEN_RESIDUE_PAIRS"));

				MODES_EIGEN_RESIDUE_AA = Integer.valueOf(parameters.get("MODES_EIGEN_RESIDUE_AA"));
				MODES_HIERARCHICAL_AA = Integer.valueOf(parameters.get("MODES_HIERARCHICAL_AA"));

				MODES_EIGEN_RESIDUE_HA = Integer.valueOf(parameters.get("MODES_EIGEN_RESIDUE_HA"));
				MODES_HIERARCHICAL_HA = Integer.valueOf(parameters.get("MODES_HIERARCHICAL_HA"));

				MODES_ALL_ATOM = Integer.valueOf(parameters.get("MODES_ALL_ATOM"));
				MODES_HEAVY_ATOM = Integer.valueOf(parameters.get("MODES_HEAVY_ATOM"));
				MODES_BACKBONE = Integer.valueOf(parameters.get("MODES_BACKBONE"));
				MODES_ALPHA_CARBON = Integer.valueOf(parameters.get("MODES_ALPHA_CARBON"));

				MODES_ATOMS_LIST = Integer.valueOf(parameters.get("MODES_ATOMS_LIST"));
				MODES_DISTANCE_PAIRS = Integer.valueOf(parameters.get("MODES_DISTANCE_PAIRS"));

				if (parameters.get("doOutputCoordinates").equals("true")) doOutputCoordinates = true;
				if (parameters.get("verbose").equals("true")) verbose = true;

				/* -------------------------------------------------------------------------------------------------------------------- */

				if (doAA) do_cartesian_AA = true;
				if (doHA) do_cartesian_HA = true;
				if (doBB) do_cartesian_BB = true;
				if (doCA) do_cartesian_CA = true;

				if (doRESIDUE_INDIVIDUAL) do_residue_individual = true;
				if (doRESIDUE_PAIRS) do_residue_pairs = true;

				if (doHIERARCHICAL_AA) do_hierarchical_AA = true;
				if (doHIERARCHICAL_HA) do_hierarchical_HA = true;

				if (doATOM_LIST) do_atom_list = true;
				if (doDIST_PAIRS) do_dist_pairs = true;

				if (!(DIRECTORY.endsWith(File.separator)))
					{
						System.err.println("Expected the directory to end with " + File.separator + ", but got: " + DIRECTORY);
						System.err.println("Attempting to fix...");
						DIRECTORY = DIRECTORY + File.separator;
					}
				boolean dir = new File(DIRECTORY).isDirectory();
				if (!dir)
					{
						System.err.println("Sorry, but the entered directory is not recognized as a proper directory: " + DIRECTORY);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				boolean check = new File(DIRECTORY + REFERENCE_PDB).exists();
				if (!check)

					{
						System.err.println("Sorry, the entered PDB Reference File can not be found: " + DIRECTORY + REFERENCE_PDB);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				check = new File(DIRECTORY + ORIGINAL_PDB_COORDS).exists();
				if (!check)
					{
						System.err.println("Sorry, the entered All-Atom Coordinates Matrix File can not be found. Please check: " + DIRECTORY + ORIGINAL_PDB_COORDS);
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
				if (MODES_EIGEN_RESIDUE_AA > 21)
					{
						MODES_EIGEN_RESIDUE_AA = 21;
						System.err.println("WARNING: Number of modes RESIDUE_AA reset to maximum allowed: 21 (Limit is glycine)");
					}
				if (MODES_EIGEN_RESIDUE_HA > 12)
					{
						MODES_EIGEN_RESIDUE_HA = 12;
						System.err.println("WARNING: Number of modes RESIDUE_HA reset to maximum allowed: 12 (Limit is glycine)");
					}
				if (doOutlierProcessing)
					{
						if (MAD_SCORE_CUTOFF > 0 & Z_SCORE_CUTOFF > 0)
							{
								System.err.println("WARNING: Only one type of outlier processing is allowed per analysis.");
								System.err.println("\tSetting MAD_SCORE_CUTOFF to zero.");
								MAD_SCORE_CUTOFF = 0;
							}
					}
				if (do_residue_pairs);
					{
						check = new File(DIRECTORY + RESIDUE_LIST_PAIRS).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_PAIRS: Sorry, the entered Residue List File can not be found. Please check: " + DIRECTORY + RESIDUE_LIST_PAIRS);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (do_residue_individual)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_INDIVIDUAL).exists();
						if (!check)
							{
								System.err
										.println("RESIDUE_LIST_LOCAL: Sorry, the entered Residue List File can not be found. Please check: " + DIRECTORY + RESIDUE_LIST_INDIVIDUAL);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (do_hierarchical_AA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_HIERARCHICAL_AA).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_HIERARCHICAL_AA: Sorry, the entered Residue List File can not be found. Please check: " + DIRECTORY
										+ RESIDUE_LIST_HIERARCHICAL_AA);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (do_hierarchical_HA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_HIERARCHICAL_HA).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_HIERARCHICAL_HA: Sorry, the entered Residue List File can not be found. Please check: " + DIRECTORY
										+ RESIDUE_LIST_HIERARCHICAL_HA);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (doAA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_ALL_ATOM).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_ALL_ATOM: Sorry, the entered All-Atom Residue List File can not be found. Please check: " + DIRECTORY
										+ RESIDUE_LIST_ALL_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (doBB)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_BACKBONE).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_BACKBONE: Sorry, the entered Backbone Residue List File can not be found. Please check: " + DIRECTORY
										+ RESIDUE_LIST_BACKBONE);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (doHA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_HEAVY_ATOM).exists();
						if (!check)
							{
								System.err.println("RESIDUE_LIST_HEAVY_ATOM: Sorry, the entered Heavy-Atom Residue List File can not be found. Please check: " + DIRECTORY
										+ RESIDUE_LIST_HEAVY_ATOM);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (doCA)
					{
						check = new File(DIRECTORY + RESIDUE_LIST_ALPHA_CARBON).exists();
						if (!check)
							{
								System.err.println(
										"RESIDUE_LIST_ALPHA_CARBON: Sorry, the entered Residue List File can not be found. Please check: " + DIRECTORY + RESIDUE_LIST_ALPHA_CARBON);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (do_atom_list)
					{
						check = new File(DIRECTORY + ATOMS_LIST).exists();
						if (!check)
							{
								System.err.println("ATOMS_LIST: Sorry, the entered Atom List File can not be found. Please check: " + DIRECTORY + ATOMS_LIST);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				if (do_dist_pairs)
					{
						check = new File(DIRECTORY + ATOM_PAIRS_LIST).exists();
						if (!check)
							{
								System.err.println("ATOM_PAIRS_LIST: Sorry, the entered Atom-Pairs List File can not be found. Please check: " + DIRECTORY + ATOM_PAIRS_LIST);
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
			}
		OUT_DIR = DIRECTORY + "JEDi_RESULTS_" + DESCRIPTION + File.separatorChar;
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
	}
}