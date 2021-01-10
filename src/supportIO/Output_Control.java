package supportIO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import Jama.Matrix;
import support.Atom;

public class Output_Control
{
	final String input_file = "JEDi_Parameters.txt", delim = "=", Q = "COV", R = "CORR", P = "PCORR", AA = "All_Atom", HA = "Heavy_Atom", BB = "BackBone", CA = "Alpha_Carbon",
			cAL = "Atom_List_PCA", cAA = "All_Atom_PCA", cBB = "Backbone_PCA", cAC = "Alpha_Carbon_PCA", cHA = "Heavy_Atom_PCA", Res_AA = "Residue_All_Atom_PCA",
			Res_PAIRS = "Residue_Pair_Analysis", HAA = "Hierarchical_All_Atom_PCA", HHA = "Hierarchical_Heavy_Atom_PCA", DP = "Distance_Pair_PCA";
	final File log;
	final NumberFormat nf0, nf3, nf6, df;
	final RoundingMode rm;
	BufferedWriter log_writer;

	// ********************************************* CONSTRUCTOR *********************************************************************************** //

	public Output_Control()
	{
		super();

		rm = RoundingMode.HALF_UP;
		df = new DecimalFormat("0.###E0");
		df.setRoundingMode(rm);

		nf0 = NumberFormat.getInstance();
		nf0.setRoundingMode(rm);
		nf0.setMaximumFractionDigits(0);
		nf0.setMinimumFractionDigits(0);

		nf3 = NumberFormat.getInstance();
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);

		nf6 = NumberFormat.getInstance();
		nf6.setMaximumFractionDigits(6);
		nf6.setMinimumFractionDigits(6);
		nf6.setRoundingMode(rm);

		log = new File(Input_Parameters.OUT_DIR + "JEDi_LOG_" + DateUtils.now2() + ".txt");

		try
			{
				log_writer = new BufferedWriter(new FileWriter(log));
			}
		catch (IOException e)
			{
				System.err.println("Could not create the JEDi LOG file. Please check directory permissions: " + Input_Parameters.OUT_DIR);
				e.printStackTrace();
			}
	}

	// *********************************************** METHODS ************************************************************************************* //

	public void initialize_JED_Log()
	{
		try
			{
				log_writer.write("JEDi: Java Essential Dynamics Inspector - Multi-Threaded" + "\n");
				log_writer.write("Release 1, December, 2020" + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();

			}
		catch (IOException e)
			{
				System.err.println("Could not write the JEDi LOG file. Please check directory permissions: " + Input_Parameters.OUT_DIR);
				e.printStackTrace();
			}
	}

	public void write_Parameter_Log(List<String> lines)
	{
		try
			{
				if (Input_Parameters.doPREPROCESS)
					{
						log_writer.write("\nParameters for the PRE PROCESSING Run.\n");
						log_writer.write("\tDIRECTORY=" + Input_Parameters.DIRECTORY + "\n");
						log_writer.write("\tDESCRIPTION=" + Input_Parameters.DESCRIPTION + "\n");
						log_writer.write("\tREFERENCE_PDB=" + Input_Parameters.REFERENCE_PDB + "\n");
						log_writer.write("\tREAD_PDBS_FILTER_STRING=" + Input_Parameters.READ_PDBS_FILTER_STRING + "\n");
						log_writer.write("----------------------------------------------------------------------------------------------\n");
						log_writer.flush();
					}
				else
					{

						int numOfLines = 0;
						log_writer.write("\nParameters for the analysis (as read from input file):\n");
						for (String line : lines)
							{
								log_writer.write("\t" + line + "\n");
								if (!(line.startsWith("#") || line.startsWith("-"))) numOfLines++;
							}
						log_writer.write("\nThe total number of KEY=VALUE pairs is " + numOfLines + "\n");
						log_writer.write("----------------------------------------------------------------------------------------------\n");
					}
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Reference_PDB_Log(int numRes, int numAtoms_AA, int numAtoms_HA)
	{
		try
			{
				log_writer.write("\nProcessed the Reference PDB file: " + Input_Parameters.REFERENCE_PDB + "\n");
				log_writer.write("\tThe number of residues found in the Reference PDB file = " + numRes + "\n");
				log_writer.write("\tThe number of atoms found in the Reference PDB file = " + numAtoms_AA + "\n");
				log_writer.write("\tThe number of Heavy atoms found in the Reference PDB file = " + numAtoms_HA + "\n");
				log_writer.write("\tCreated file: 'All_PDB_Residues_JEDi.txt'. This file contains all chainID-residue number pairs found in the Reference PDB file." + "\n");
				log_writer.write("\tNumber of atoms per residue list file created: 'numbers_Of_Atoms_in_Residues.txt' " + "\n");
				log_writer.write("\tNumber of heavy atoms per residue list file created: 'numbers_Of_Heavy_Atoms_in_Residues.txt' " + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();

			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Coordinate_File_Log(int rows, int cols)
	{
		try
			{
				log_writer.write("\nRead the All Atom PDB coordinates file: " + Input_Parameters.ORIGINAL_PDB_COORDS + "\n");
				log_writer.write("\tThe expected packing in this Matrix is: {X}{Y}{Z} stacking." + "\n");
				log_writer.write("\tThe dimension of the coordinates matrix is = " + rows + " by " + cols + "\n");
				log_writer.write("\tTotal number of atoms in matrix = " + (rows / 3) + "\n");
				log_writer.write("\tTotal number of conformations in matrix = " + cols + "\n");
				log_writer.write("\tThe matching reference structure is: " + Input_Parameters.REFERENCE_PDB + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();

			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Down_Sample_Log(int cols_orig, int cols_final)
	{
		try
			{
				log_writer.write("\nDown Sampled the coordinates file:\n");
				log_writer.write("\tThe down-sample-factor was: " + Input_Parameters.DOWNSAMPLE + "\n");
				log_writer.write("\tThe number of frames in the original matrix was: " + cols_orig + "\n");
				log_writer.write("\tThe number of frames in the reduced matrix is: " + cols_final + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Frame_Select_Log(int cols_orig)
	{
		try
			{
				log_writer.write("\nSelected frames from the coordinates file:\n");
				log_writer.write("\tFirst frame is: " + Input_Parameters.FRAME_START + " The last frame is: " + Input_Parameters.FRAME_END + "\n");
				log_writer.write("\tThe number of frames in the original matrix was: " + cols_orig + "\n");
				log_writer.write("\tThe number of frames in the reduced matrix is: " + (Input_Parameters.FRAME_END - Input_Parameters.FRAME_START + 1) + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_PreProcess_Log(int numRes, int numAtoms, int numConformations)
	{
		try
			{
				log_writer.write("\nPerformed a PRE PROCESSING Run.\n");
				log_writer.write("\tThe Working Directory:\n" + "\t" + Input_Parameters.DIRECTORY + "\n");
				if (Input_Parameters.doREAD_ARCHIVE) log_writer.write("\tRead PDB files from Archive:\n" + "\t" + Input_Parameters.ARCHIVE_NAME + "\n");
				else
					log_writer.write("\tRead PDB files from the Working Directory:\n");
				log_writer.write("\t\tUsed the PDB Filter String: " + Input_Parameters.READ_PDBS_FILTER_STRING + "\n");
				log_writer.write("\tThe names of the PDB files read were logged to the file: 'PDB_Read_Log.txt' \n");
				log_writer.write("\tThe number of PDB files read = " + numConformations + "\n");
				log_writer.write("\tThe number of residues found in the PDB files = " + numRes + "\n");
				log_writer.write("\tThe number of atoms found in the PDB files = " + numAtoms + "\n");
				log_writer.write("\tAll Atom Coordinates Matrix created: 'original_PDB_Coordinates_AA.txt'" + "\n");
				log_writer.write("\tThe packing in this Matrix is: {X}{Y}{Z} stacking." + "\n");
				log_writer.write("\tThe dimension of this matrix is: " + (numAtoms * 3) + "  by  " + numConformations + "\n");
				log_writer.write("\tPDB reference structure is: " + Input_Parameters.REFERENCE_PDB + "\n");
				log_writer.write("\tConformation RMSDs and Atom RMSFs were calculated and plotted." + "\n");
				log_writer.write("\tA PDB file with B-factors replaced by atomic RMSFs was created: \n");
				log_writer.write("\t\tss_" + numRes + "_" + numAtoms + "_RMSF_edited.pdb" + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Stat_Thresholding_Log(String type)
	{
		try
			{
				log_writer.write("\nThe " + type + " Variable Statistics were calculated and plotted including: \n");
				log_writer.write("\tThe Means \n");
				log_writer.write("\tThe Variances \n");
				log_writer.write("\tThe Skews \n");
				log_writer.write("\tThe Kurtosis \n");
				log_writer.write("\tAtom lists and Residue lists were generated for each atom and residue that exceeded the entered thresholds. \n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Subset_Log(String type, int numRes, int numAtoms, Matrix ref, Matrix coords, List<Atom> atoms)
	{
		try
			{
				log_writer.write("\nThe " + type + " Subset coordinates were obtained from: " + Input_Parameters.RESIDUE_LIST_ALL_ATOM + "\n");
				log_writer.write("\tThe number of residues in the subset = " + numRes + "\n");
				log_writer.write("\tThe number of atoms in the subset = " + numAtoms + "\n");
				log_writer.write("\tThe dimension of the Reference Coordinates matrix  = " + ref.getRowDimension() + " by " + ref.getColumnDimension() + "\n");
				log_writer.write("\tThe dimension of the Coordinates matrix  = " + coords.getRowDimension() + " by " + coords.getColumnDimension() + "\n");
				log_writer.write("\tThe coordinates were aligned to the reference coordinates using quaternion algebra. \n");
				log_writer.write("\tThe atomic RMSFs were calculated and used to create a PDB file with B-Factors replaced with the atomic RMSF. \n");
				log_writer.write("\tThe subset atoms are:\n");
				for (Atom atm : atoms)
					{
						log_writer.write("\t\t" + atm.toString() + "\n");
					}
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Subset_Log_AL(int numAtoms, Matrix ref, Matrix coords, List<Atom> atoms)
	{
		try
			{
				log_writer.write("\nThe Atom List Subset coordinates were obtained from: " + Input_Parameters.ATOMS_LIST + "\n");
				log_writer.write("\tThe number of atoms in the Atom List Subset = " + numAtoms + "\n");
				log_writer.write("\tThe dimension of the Reference Coordinates matrix  = " + ref.getRowDimension() + " by " + ref.getColumnDimension() + "\n");
				log_writer.write("\tThe dimension of the Coordinates matrix  = " + coords.getRowDimension() + " by " + coords.getColumnDimension() + "\n");
				log_writer.write("\tThe coordinates were aligned to the reference coordinates using quaternion algebra. \n");
				log_writer.write("\tThe atomic RMSFs were calculated and used to create a PDB file with B-Factors replaced with the atomic RMSF. \n");
				log_writer.write("\t\tss_" + numAtoms + "_AL_RMSF_edited.pdb" + "\n");
				log_writer.write("\tThe Atom List Subset Reference atoms are:\n");
				for (Atom atm : atoms)
					{
						log_writer.write("\t\t" + atm.toShortString() + "\n");
						System.out.println("\t\t" + atm.toShortString());
					}
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Subset_Log_DP(int numAtomPairs, List<Atom> atoms1, List<Atom> atoms2)
	{
		try
			{
				log_writer.write("\nThe Distance-Pair Subset coordinates were obtained from: " + Input_Parameters.ATOM_PAIRS_LIST + "\n");
				log_writer.write("\tThe number of atom pairs in the Distance Pairs Subset = " + numAtomPairs + "\n");
				log_writer.write("\tThe Reference Atoms Pairs are:\n");
				for (int i = 0; i < atoms1.size(); i++)
					{
						Atom a1 = atoms1.get(i);
						Atom a2 = atoms2.get(i);
						log_writer.write("\t\t" + a1.toShortString() + " <----->    " + a2.toShortString() + "\n");
						System.out.println("\t\t" + a1.toShortString() + " <----->    " + a2.toShortString());
					}
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_PCA_Log(String type, int numRes, int numAtoms, int numModes, int rank, double trace, double cond, double det, double kmo, double shrinkage)
	{
		try
			{
				log_writer.write("\nPerformed " + type + ". Computed Top " + numModes + " modes." + "\n\n");
				log_writer.write("\tNumber of Residues in the subset: " + numRes + "\n");
				log_writer.write("\tNumber of Atoms in the subset: " + numAtoms + "\n");
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						log_writer.write("\tThe coordinate values with MAD-scores GREATER THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Removal.\n");
						log_writer.write("\tThe coordinate values with MAD-scores LESS THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Selection.\n");
					}
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						log_writer.write(
								"\tThe coordinate values with Z-scores GREATER THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Removal.\n");
						log_writer.write(
								"\tThe coordinate values with Z-scores LESS THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Selection.\n");
					}
				else
					log_writer.write("\tNo outlier processing was done.\n");
				if (Input_Parameters.THRESHOLD_COV > 0)
					log_writer.write("\tThe covariance matrix was numerically stabililzed at threshold " + df.format(Input_Parameters.THRESHOLD_COV) + "\n");
				if (Input_Parameters.doSPARSIFY)
					{
						if (Input_Parameters.THRESHOLD_CORR > 0)
							log_writer.write("\tThe correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_CORR + "\n");
						if (Input_Parameters.THRESHOLD_PCORR > 0)
							log_writer.write("\tThe partial correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_PCORR + "\n");
					}
				log_writer.write("\tThe KMO Statistic for the variables in the subset is: " + nf6.format(kmo) + "\n");
				log_writer.write("\tShrinkage intensity used to shrink the Covariance Matrix = " + nf6.format(shrinkage) + "\n");
				log_writer.write("\tRank of the Covariance Matrix = " + (rank) + "\n");
				log_writer.write("\tCondition Number of the Covariance Matrix = " + df.format(cond) + "\n");
				log_writer.write("\tTrace of the Covariance Matrix = " + df.format(trace) + "\n");
				log_writer.write("\tDeterminant of the Covariance Matrix = " + df.format(det) + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Hierarchical_PCA_Log(String type, int numRes, int numAtoms, int numModesResidue, int numModesHierarchical, int rank, double trace, double cond, double det,
			double shrinkage)
	{
		try
			{
				log_writer.write("\nPerformed " + type + ". Computed Top " + numModesHierarchical + " modes." + "\n\n");
				log_writer.write("\tThe Eigen Residues and residue PCs were calculated by performing residue pca using a global alignment over the entire subset.\n");
				log_writer.write("\t\tEach residue was represented with the top " + numModesResidue + " PCs (DOFs).\n");
				log_writer.write("\tNumber of Residues in the subset: " + numRes + "\n");
				log_writer.write("\tNumber of Atoms in the subset: " + numAtoms + "\n");
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						log_writer.write("\tThe coordinate values with MAD-scores GREATER THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Removal.\n");
						log_writer.write("\tThe coordinate values with MAD-scores LESS THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Selection.\n");
					}
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						log_writer.write(
								"\tThe coordinate values with Z-scores GREATER THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Removal.\n");
						log_writer.write(
								"\tThe coordinate values with Z-scores LESS THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Selection.\n");
					}
				else
					log_writer.write("\tNo outlier processing was done.\n");
				if (Input_Parameters.THRESHOLD_COV > 0)
					log_writer.write("\tThe covariance matrix was numerically stabililzed at threshold " + df.format(Input_Parameters.THRESHOLD_COV) + "\n");
				if (Input_Parameters.doSPARSIFY)
					{
						if (Input_Parameters.THRESHOLD_CORR > 0)
							log_writer.write("\tThe correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_CORR + "\n");
						if (Input_Parameters.THRESHOLD_PCORR > 0)
							log_writer.write("\tThe partial correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_PCORR + "\n");
					}
				// log_writer.write("\tThe KMO Statistic for the variables in the subset is: " + nf6.format(kmo) + "\n");
				log_writer.write("\tShrinkage intensity used to shrink the Covariance Matrix = " + nf6.format(shrinkage) + "\n");
				log_writer.write("\tRank of the Covariance Matrix = " + (rank) + "\n");
				log_writer.write("\tCondition Number of the Covariance Matrix = " + df.format(cond) + "\n");
				log_writer.write("\tTrace of the Covariance Matrix = " + df.format(trace) + "\n");
				log_writer.write("\tDeterminant of the Covariance Matrix = " + df.format(det) + "\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_DP_PCA_Log(String type, int numPairs, int rank, double trace, double cond, double det, double kmo, double shrinkage, List<String> IDs1, List<Integer> List1,
			List<String> IDs2, List<Integer> List2, double[] means, double[] std_devs)
	{
		try
			{
				log_writer.write("\nPerformed Distance-Pair PCA, " + type + ". Computed Top " + Input_Parameters.MODES_DISTANCE_PAIRS + " modes.\n");
				log_writer.write("\tNumber of atom pairs: " + numPairs + "\n");
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						log_writer.write("\tThe coordinate values with MAD-scores GREATER THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Removal.\n");
						log_writer.write("\tThe coordinate values with MAD-scores LESS THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Selection.\n");
					}
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						log_writer.write(
								"\tThe coordinate values with Z-scores GREATER THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Removal.\n");
						log_writer.write(
								"\tThe coordinate values with Z-scores LESS THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Selection.\n");
					}
				else
					log_writer.write("\tNo outlier processing was done.\n");
				if (Input_Parameters.THRESHOLD_COV > 0)
					log_writer.write("\tThe covariance matrix was numerically stabililzed at threshold " + df.format(Input_Parameters.THRESHOLD_COV) + "\n");
				if (Input_Parameters.doSPARSIFY)
					{
						if (Input_Parameters.THRESHOLD_CORR > 0)
							log_writer.write("\tThe correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_CORR + "\n");
						if (Input_Parameters.THRESHOLD_PCORR > 0)
							log_writer.write("\tThe partial correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_PCORR + "\n");
					}
				log_writer.write("\tShrinkage intensity used to shrink the Covariance Matrix = " + nf6.format(shrinkage) + "\n");
				log_writer.write("\tThe KMO Statistic for the variables in the dataset is: " + nf6.format(kmo) + "\n");
				log_writer.write("\tRank of the Covariance Matrix = " + (rank) + "\n");
				log_writer.write("\tCondition Number of the Covariance Matrix = " + df.format(cond) + "\n");
				log_writer.write("\tTrace of the Covariance Matrix = " + df.format(trace) + "\n");
				log_writer.write("\tDeterminant of the Covariance Matrix = " + df.format(det) + "\n");
				log_writer.write("\n\tMEANs and STANDARD DEVIATIONs for the Atom Pair Distances: " + "\n\n");
				log_writer.write(String.format("%-5s%-20s%-20s%-20s%-20s%n", "", "Atom1", "Atom2", "Mean", "Std_Dev"));
				for (int i = 0; i < numPairs; i++)
					{
						log_writer.write(String.format("%-5s%-20s%-20s%-20s%-20s%n", "", IDs1.get(i) + List1.get(i), IDs2.get(i) + List2.get(i), nf3.format(means[i]),
								nf6.format(std_devs[i])));
					}
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();

			}
		catch (IOException e)
			{
				System.err.println("Could not write the JEDi_LOG file: " + Input_Parameters.OUT_DIR + "JEDi_LOG.txt");
				e.printStackTrace();
			}
	}

	public void write_Individual_Residue_PCA_Log()
	{
		try
			{
				log_writer.write("\nPerformed Individual Residue PCA using Local Alignment for each residue, Computed Top " + Input_Parameters.MODES_RESIDUE_INDIVIDUAL + " modes."
						+ "\n\n");
				log_writer.write("\tThe results for each residue was put in a separate directory in the JEDi output tree.\n");
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.MAD_SCORE_CUTOFF > 0)
					{
						log_writer.write("\tThe coordinate values with MAD-scores GREATER THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Removal.\n");
						log_writer.write("\tThe coordinate values with MAD-scores LESS THAN |" + Input_Parameters.MAD_SCORE_CUTOFF
								+ "| were set to their median value for Outlier Selection.\n");
					}
				if (Input_Parameters.doOutlierProcessing && Input_Parameters.Z_SCORE_CUTOFF > 0)
					{
						log_writer.write(
								"\tThe coordinate values with Z-scores GREATER THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Removal.\n");
						log_writer.write(
								"\tThe coordinate values with Z-scores LESS THAN |" + Input_Parameters.Z_SCORE_CUTOFF + "| were set to their mean value for Outlier Selection.\n");
					}
				else
					log_writer.write("\tNo outlier processing was done.\n");
				if (Input_Parameters.THRESHOLD_COV > 0)
					log_writer.write("\tThe covariance matrix was numerically stabililzed at threshold " + df.format(Input_Parameters.THRESHOLD_COV) + "\n");
				if (Input_Parameters.doSPARSIFY)
					{
						if (Input_Parameters.THRESHOLD_CORR > 0)
							log_writer.write("\tThe correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_CORR + "\n");
						if (Input_Parameters.THRESHOLD_PCORR > 0)
							log_writer.write("\tThe partial correlation matrix was sparsified using the threshold: : " + Input_Parameters.THRESHOLD_PCORR + "\n");
					}
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_Residue_Pair_Analysis_Log()
	{
		try
			{
				log_writer.write("\nPerformed RESIDUE PAIR ANALYSIS using Global Alignment for each Residue Pair. \n");
				log_writer.write("\tUsed the top " + Input_Parameters.MODES_EIGEN_RESIDUE_PAIRS + " residue PCs (DOFs) to compute the coupling score. \n");
				log_writer.write("\tNo coordinate outliers were adjusted and no thresholding was done.\n");
				log_writer.write("\tThe coupling scores matrix was written to file and plotted as a heatmap in the Residue Pair Analysis Directory in the JEDi output tree.\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_FES_Log()
	{
		try
			{
				log_writer.write("\nFree-Energy calculations were done using two order parameters to create a Free-Energy Surface: \n");
				log_writer.write("\tThe first and second PCs or DVPs (as normed-projections). \n");
				log_writer.write("\tA 3D ScatterPlot in PNG format was generated for each analysis. \n");
				log_writer.write("\tSee the FES Logs in the FES subdirectory for more details.\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_KPCA_Log()
	{
		try
			{
				log_writer.write("\nKernel PCA Analysis was performed on the PCA Reduced Data: \n");
				log_writer.write("\tUsing the top 2 DVPs or PCs (as weighted-projections) \n");
				log_writer.write("\tThe number of frames is capped at: " + Input_Parameters.numberEssentialFrames + "\n"); // Need to ass this parameter!
				log_writer.write("\tAll Q and K matrices were optimized using a shrinkage method, either optimal or fixed at 0.250 \n");
				log_writer.write("\tThe following kernels were applied: \n");
				if (Input_Parameters.Cauchy) log_writer.write("\t\tCauchy \n");
				if (Input_Parameters.Circular) log_writer.write("\t\tCircular \n");
				if (Input_Parameters.Gaussian) log_writer.write("\t\tGaussian (Radial Basis Function, Power = 2) \n");
				if (Input_Parameters.Linear) log_writer.write("\t\tLinear \n");
				if (Input_Parameters.Log) log_writer.write("\t\tLog, Power = 2 \n");
				if (Input_Parameters.MI) log_writer.write("\t\tMutual Information \n");
				if (Input_Parameters.MI_KDE) log_writer.write("\t\tMutual Information with KDE \n");
				if (Input_Parameters.Degree_2_Poly) log_writer.write("\t\tPolynomial Degree 2 \n");
				if (Input_Parameters.Degree_3_Poly) log_writer.write("\t\tPolynomial Degree 3 \n");
				if (Input_Parameters.Degree_4_Poly) log_writer.write("\t\tPolynomial Degree 4 \n");
				if (Input_Parameters.XY_Poly) log_writer.write("\t\tPolynomial XY \n");
				if (Input_Parameters.Sigmoid) log_writer.write("\t\tSigmoid (Hyperbolic Tangent Function) \n");
				if (Input_Parameters.Euclidean) log_writer.write("\t\tEuclidean Distance (Using the similarity metric) \n");
				if (Input_Parameters.Mahalanobis) log_writer.write("\t\tMahalanobis Distance (Using the similarity metric) \n");
				log_writer.write("\tXY ScatterPlots were generated in PNG format for each of the kernel's: Showing the top 2 kPCs. \n");
				log_writer.write("\tNOTE: NO OPTIMIZATION of the kernel parameters was done... \n");
				log_writer.write("\tThis analysis is only meant to provide a high level view of patterns in the data. \n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_DVP_Log()
	{
		try
			{
				log_writer.write("\nThe Displacement Vectors (DVs) were calculated using the reference PDB: " + Input_Parameters.REFERENCE_PDB + "\n");
				log_writer.write("\nThe Displacement Vector Projections (DVPs) from the selected PCA models were calculated in the following ways:" + "\n");
				log_writer.write("\tStandard dot product(dp), normed dp, weighted dp (by sqrt of eigenvalue), and weighted normed dp" + "\n");
				log_writer.write("\tThe weighted DVPs were plotted as 2-D Scatter Plots (PC Plots).\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_SSA_Log(boolean inter_ssa_aa, boolean inter_ssa_ha)
	{
		try
			{
				if (Input_Parameters.doCORR || Input_Parameters.doPCORR)
					{
						log_writer.write("\n'Intra-Subset' Subspace Analysis was done for each subset comparing the top vector spaces from the selected PCA models." + "\n");
						if (Input_Parameters.doSPARSIFY) log_writer.write("\tSparse matrices were compared to their non-sparse counterparts.\n");
						log_writer.write("\tNote that the Hierarchical Methods only use the Covariance PCA model.\n");
						log_writer.write("\tKey metrics include RMSIP and Principle Angle Spectra. \n");
						log_writer.write("\tThe iterated RMSIPs with comparisons to random were plotted for each analysis with Z-Scores for assessing significance. \n");
						log_writer.write("\tAdditional log files can be found in the /SSA directory tree." + "\n");
						log_writer.write("----------------------------------------------------------------------------------------------\n");
						log_writer.flush();
					}
				if (inter_ssa_aa)
					{
						log_writer.write("\n'Inter-Subset' Subspace Analysis was done for Hierarchical and Direct All Atom Subsets, COV Model." + "\n");
						log_writer.write("\tThe iterated RMSIPs with comparisons to random were plotted along with Z-Scores for assessing significance. \n");
						log_writer.write("\tAdditional log files can be found in the /SSA directory tree." + "\n");
						log_writer.write("----------------------------------------------------------------------------------------------\n");
						log_writer.flush();
					}
				if (inter_ssa_ha)
					{
						log_writer.write("\n'Inter-Subset' Subspace Analysis was done for Hierarchical and Direct Heavy Atom Subsets, COV Model." + "\n");
						log_writer.write("\tThe iterated RMSIPs with comparisons to random were plotted along with Z-Scores for assessing significance. \n");
						log_writer.write("\tAdditional log files can be found in the /SSA directory tree." + "\n");
						log_writer.write("----------------------------------------------------------------------------------------------\n");
						log_writer.flush();
					}
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public void write_VIZ_Log()
	{
		try
			{
				if (Input_Parameters.doEssentialViz)
					{
						log_writer.write("\nEssential Mode Visualization done for selected PCA models:\n");
						log_writer.write("\tTo construct the essential motion, a super-position of the top modes was created spanning multiple full periods of the first mode:\n");
						log_writer.write("\t\tNumber of modes in the superposition: " + Input_Parameters.numberModeComponents + "\n");
						log_writer.write("\t\tNumber of Cycles for the superposition: " + Input_Parameters.numberEssentialCycles + "\n");
						log_writer.write("\t\tNumber of Frames for the superposition: " + Input_Parameters.numberEssentialFrames + "\n");
						log_writer.write("\tA PyMol(TM) scripts was generated for each essential subspace, to play the mode structures as movies.\n");
					}
				if (Input_Parameters.doModeViz)
					{
						log_writer.write("\nPerforming Individual Mode Visualization on the Top Selected PCA modes:\n");
						log_writer.write("\tTo construct the PCA modes, atoms of each residue were perturbed along the mode eigenvector using a sine function.\n");
						log_writer.write("\t\tNumber of modes to visualize: " + Input_Parameters.MODES_VIZ + "\n");
						log_writer.write("\t\tNumber of Cycles for each mode: " + Input_Parameters.numberModeCycles + "\n");
						log_writer.write("\t\tNumber of Frames for the mode: " + Input_Parameters.numberModeFrames + "\n");
						log_writer.write("\tPyMol(TM) scripts were generated for each individual mode, to play the mode structures as movies.\n");
					}
				log_writer.write("\t\tMODE AMPLITUDE = " + nf3.format(Input_Parameters.VIZ_MODE_SCALE_FACTOR) + "\n");
				log_writer.write("\t\tColor FLOOR for the log scale coloring = " + nf3.format(Input_Parameters.LOG_FLOOR) + "\n");
				log_writer.write("\t\tNote: Hierarchical PCA Methods only use the Covariance PCA model.\n");
				log_writer.write("\t\tNote: The Distance Pair PCA Method is NOT visualized.\n");
				log_writer.write("----------------------------------------------------------------------------------------------\n");
				log_writer.flush();
			}
		catch (IOException e)
			{
				System.err.println("Could not write the JEDi_LOG file: " + Input_Parameters.OUT_DIR + "JEDi_LOG.txt");
				e.printStackTrace();
			}
	}

	public void close_JEDi_Log(String date, long time)
	{
		try
			{
				log_writer.write("\nJEDi Analysis Completed: " + date);
				log_writer.write("\n\tTotal Run Time = " + nf3.format(time / (60.000E9)) + " minutes.");
				log_writer.close();

			}
		catch (IOException e)
			{
				System.err.println("Could not close the JEDi Log file...");
				e.printStackTrace();
			}
	}
}

