package jedi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

import Jama.Matrix;
import support.Atom;
import support.Atom_ID_Pair;
import support.Residue_ID_Pair;
import supportIO.Input_Parameters;
import supportIO.List_IO;
import supportIO.PDB_Filter;
import supportIO.PDB_IO;
import supportIO.ZIP_IO;

public class JEDi_Do_PDB_Processing
{
	final boolean read_archive;
	int number_of_atoms, number_of_atoms_AA, number_of_atoms_BB, number_of_atoms_HA, number_of_residues, number_of_atom_pairs, number_of_residues_SS;
	final double FLOOR = 1.00E-6, delta_y = 99;
	final String filter_String, directory, pdb_ref_file;
	final String C = "C", CA = "CA", N = "N", O = "O", H = "H", P = "P", S = "S", M = "M";
	String out_dir, archive_name;

	List<Integer> atom_list, atom_list_original, atom_list1, atom_list1_original, atom_list2, atom_list2_original;
	List<Integer> number_of_Atoms_in_Residues, number_of_Heavy_Atoms_in_Residues;
	List<Integer> atom_numbers_read_original, backbone_atom_numbers_read_original, heavy_atom_numbers_read_original, alpha_carbon_numbers_read_original,
			residue_numbers_read_original;
	List<Integer> SS_atom_numbers_original, SS_atom_numbers, SS_backbone_atom_numbers_original, SS_backbone_atom_numbers, SS_heavy_atom_numbers_original, SS_heavy_atom_numbers,
			SS_alpha_carbon_numbers_original, SS_alpha_carbon_numbers;
	final List<Integer> atom_numbers_read, heavy_atom_numbers_read, backbone_atom_numbers_read, alpha_carbon_numbers_read, residue_numbers_read;

	List<String> chain_ids, chain_ids0, chain_ids1, chain_ids2, chain_IDs_read, PDB_file_names;
	final List<String> ref_PDB_file_lines;

	Matrix reference_PDB_coordinates_AA, reference_subset_PDB_coordinates_AA, original_PDB_coordinates_AA, subset_PDB_coordinates_AA;
	Matrix reference_PDB_coordinates_BB, reference_subset_PDB_coordinates_BB, subset_PDB_coordinates_BB;
	Matrix reference_PDB_coordinates_CA, reference_subset_PDB_coordinates_CA, subset_PDB_coordinates_CA;
	Matrix reference_PDB_coordinates_HA, reference_subset_PDB_coordinates_HA, subset_PDB_coordinates_HA;
	Matrix subset_Hierarchical_Coordinates;

	List<Atom> atoms_AA, atoms_BB, atoms_CA, atoms_HA;
	List<Atom> ref_atoms_AA, ref_atoms_BB, ref_atoms_CA, ref_atoms_HA;
	List<Atom> ref_subset_atoms;

	List<Residue_ID_Pair> residue_ID_pairs_read, residue_list, residue_list_orig;
	List<Atom_ID_Pair> atom_ID_pairs_read;

	List<Matrix> Residue_Coords_List, Residue_Coords_HA_List, Reference_Residue_Coords_List;
	List<List<Atom>> Residue_Atoms_List;

	StringTokenizer sToken;
	final PDB_IO pio;
	final ZIP_IO zio;

	/* ************************************************* CONSTRUCTOR **************************************************************************** */

	public JEDi_Do_PDB_Processing()
	{
		this.directory = Input_Parameters.DIRECTORY;
		this.pdb_ref_file = Input_Parameters.REFERENCE_PDB;
		this.filter_String = Input_Parameters.READ_PDBS_FILTER_STRING;

		this.read_archive = Input_Parameters.doREAD_ARCHIVE;
		this.zio = new ZIP_IO();

		if (read_archive) this.archive_name = Input_Parameters.ARCHIVE_NAME;

		if (pdb_ref_file.endsWith("bz2") || pdb_ref_file.endsWith("gz"))
			{
				ref_PDB_file_lines = zio.read_Lines_from_Compressed_File(directory + pdb_ref_file);
			}
		else if (pdb_ref_file.endsWith("tar.bz2") || pdb_ref_file.endsWith("tar.gz"))
			{
				ref_PDB_file_lines = zio.read_Lines_from_TAR_File(directory + pdb_ref_file);
			}
		else if (pdb_ref_file.endsWith("zip"))
			{
				ref_PDB_file_lines = zio.read_Lines_from_ZIP_File(directory + pdb_ref_file);
			}
		else
			ref_PDB_file_lines = List_IO.read_Lines_From_File(directory + pdb_ref_file);

		this.pio = new PDB_IO(ref_PDB_file_lines);

		atom_numbers_read = new ArrayList<Integer>(10000);
		backbone_atom_numbers_read = new ArrayList<Integer>(2500);
		heavy_atom_numbers_read = new ArrayList<Integer>(5000);
		alpha_carbon_numbers_read = new ArrayList<Integer>(1000);
		residue_numbers_read = new ArrayList<Integer>(1000);

		number_of_Atoms_in_Residues = new ArrayList<Integer>(1000);
		number_of_Heavy_Atoms_in_Residues = new ArrayList<Integer>(1000);
		chain_ids = new ArrayList<String>(1000); // list of all chain identifiers specified
	}

	/* ************************************************* READ PDB FILE METHODS **************************************************************************** */

	/**
	 * Reads the Reference PDB file and extracts all needed information for JEDi analysis.
	 * 
	 */
	public void parse_Reference_PDB()
	{
		ref_atoms_AA = pio.parse_PDB_Lines_List();
		ref_atoms_BB = pio.get_Backbone_Atoms();
		ref_atoms_CA = pio.get_Atoms_CA();
		ref_atoms_HA = pio.get_Heavy_Atoms();

		atom_numbers_read_original = pio.get_Atom_Numbers_Read();
		for (int number : atom_numbers_read_original) // Need to adjust in case of Non-Contiguous atom numbering and JAN to be consistent with coordinates matrix
			{
				int index = atom_numbers_read_original.indexOf(number);
				atom_numbers_read.add(index);
			}

		backbone_atom_numbers_read_original = pio.get_Backbone_Atom_Numbers_Read();
		for (int number : backbone_atom_numbers_read_original)
			{
				int index = backbone_atom_numbers_read_original.indexOf(number);
				backbone_atom_numbers_read.add(index);
			}

		heavy_atom_numbers_read_original = pio.get_Heavy_Atom_Numbers_Read();
		for (int number : heavy_atom_numbers_read_original)
			{
				int index = heavy_atom_numbers_read_original.indexOf(number);
				heavy_atom_numbers_read.add(index);
			}

		alpha_carbon_numbers_read_original = pio.get_Alpha_Carbon_Numbers_Read();
		for (int number : alpha_carbon_numbers_read_original)
			{
				int index = alpha_carbon_numbers_read_original.indexOf(number);
				alpha_carbon_numbers_read.add(index);
			}

		residue_numbers_read_original = pio.get_Residue_Numbers_Read();
		chain_IDs_read = pio.get_Chain_IDs_Read();

		residue_ID_pairs_read = pio.getResidue_ID_pairs();
		atom_ID_pairs_read = pio.getAtom_ID_pairs();

		for (Residue_ID_Pair pair : residue_ID_pairs_read)
			{
				int numAtomsAA = 0;
				int numAtomsHA = 0;

				// System.out.println("Checking: " + pair.toString());

				String id = pair.getChain_ID();
				int num = pair.getResidue_Number();
				for (Atom a : ref_atoms_AA)
					{
						if (a.chainID.equals(id) && a.res_number == num)
							{
								// System.out.println("Match to Residue ID Pair: " + a.chainID + "\t" + a.res_number);
								numAtomsAA++;
								// if (!(a.symbol.contains("H")) || a.symbol.contains("CH") || a.symbol.contains("NH"))
								if (a.symbol.startsWith(C) || a.symbol.startsWith(N) || a.symbol.startsWith(O) || a.symbol.startsWith(P) || a.symbol.startsWith(S)
										|| a.symbol.startsWith(M))
									{
										numAtomsHA++;
									}
							}
					}
				number_of_Atoms_in_Residues.add(numAtomsAA);
				number_of_Heavy_Atoms_in_Residues.add(numAtomsHA);
			}

		reference_PDB_coordinates_AA = get_Cartesian_PDB_Coords(ref_atoms_AA);
		reference_PDB_coordinates_BB = get_Cartesian_PDB_Coords(ref_atoms_BB);
		reference_PDB_coordinates_CA = get_Cartesian_PDB_Coords(ref_atoms_CA);
		reference_PDB_coordinates_HA = get_Cartesian_PDB_Coords(ref_atoms_HA);

		number_of_atoms_AA = ref_atoms_AA.size();
		number_of_atoms_BB = ref_atoms_BB.size();
		number_of_atoms_HA = ref_atoms_HA.size();

		number_of_residues = residue_ID_pairs_read.size();
	}

	/**
	 * This Method is used for building original_PDB_coordinates_AA while reading and parsing the PDB files in the working directory.
	 * 
	 * Processes all atoms in specified PDB file (frame)
	 * 
	 * @param directory
	 * @param pdb_file
	 */
	private List<Atom> parse_PDB_as_File(String directory, String pdb_file)
	{
		PDB_IO pio_AA = new PDB_IO(directory + pdb_file);
		return pio_AA.parse_PDB_File();
	}

	/**
	 * This Method is used for building original_PDB_coordinates_AA while parsing the lines of the PDB files in the working directory that have been read into a list of strings.
	 * 
	 * Processes all atoms in specified PDB file (frame)
	 * 
	 * @param directory
	 * @param pdb_file
	 */
	private List<Atom> parse_PDB_as_Lines_List(List<String> lines)
	{
		PDB_IO pio_AA = new PDB_IO(lines);
		return pio_AA.parse_PDB_Lines_List();
	}

	/* ************************************************* GET MATRIX OF COORDINATES METHODS **************************************************************************** */

	/**
	 * Returns the matrix of the all atom coordinates from all PDB files in specified directory (Subject to matching the PDB_FILTER_STRING)
	 * 
	 * Supported file types are '.pdb', '.pdb.bz2', '.pdb.gz', and '.pdb.zip' (A SINGLE compressed file, NOT a directory archive)
	 * 
	 * NOTE: These MUST be SINGLE FILES (A single entry zip or tar.gz or tar.bz2 archive is fine)
	 * 
	 * TO PROCESS A MULTI FILE ARCHIVE, the 'doARCHIVE=TRUE' flag must be set in the Input Parameters File.
	 * 
	 * @return original_PDB_coordinates_AA
	 */
	public Matrix get_Matrix_of_Original_PDB_Coords_from_PDB_Files()
	{
		File pdb_file = new File(directory);
		FilenameFilter filter = new PDB_Filter(filter_String);
		String[] ls = pdb_file.list(filter);
		Arrays.sort(ls);
		int number_of_conformations = ls.length;
		int ROWS_AA = (number_of_atoms_AA * 3);
		List<Atom> atoms = null;

		PDB_file_names = new ArrayList<>(ls.length);
		original_PDB_coordinates_AA = new Matrix(ROWS_AA, number_of_conformations);

		for (int i = 0; ls != null && i < ls.length;) for (i = 0; i < ls.length; i++)
			{
				System.out.println("PDB file: " + ls[i]);
				PDB_file_names.add(ls[i]);
				if (ls[i].endsWith("pdb")) atoms = parse_PDB_as_File(directory, ls[i]);
				if (ls[i].endsWith("bz2") || ls[i].endsWith("gz")) atoms = parse_PDB_as_Lines_List(zio.read_Lines_from_Compressed_File(directory + ls[i]));
				if (ls[i].endsWith("zip")) atoms = parse_PDB_as_Lines_List(zio.read_Lines_from_ZIP_File(directory + ls[i]));
				if (ls[i].endsWith("tar.bz2") || ls[i].endsWith("tar.gz")) atoms = parse_PDB_as_Lines_List(zio.read_Lines_from_TAR_File(directory + ls[i]));

				if (Input_Parameters.doParityCheck)
					{
						// Verify that each PDB file read is equivalent to the reference PDB file:
						// Here, minimum parity means same atom type, same residue type, and same residue number.
						// Atom numbers may be different... This is needed for pooling mutants with INDELs.
						for (int k = 0; k < atoms.size(); k++)
							{
								Atom refAtom = ref_atoms_AA.get(k);
								Atom test = atoms.get(k);
								boolean ok = refAtom.minimumParity(refAtom, test);
								if (!ok)
									{
										System.err.println("The PDB file " + ls[i] + " differes from reference at position " + (k + 1));
										System.err.println(refAtom.toVeryShortString());
										System.err.println(test.toVeryShortString());
									}
							}
					}

				Matrix cv_aa = get_Cartesian_PDB_Coords(atoms);
				original_PDB_coordinates_AA.setMatrix(0, ROWS_AA - 1, i, i, cv_aa);
			}
		return original_PDB_coordinates_AA;
	}

	/**
	 * Returns the matrix of the all atom coordinates from all PDB files in specified a TAR.BZ2 or TAR.GZ Archive
	 * 
	 * NOTE: The archive MUST be a SINGLE DIRECTORY containing all the PDB files to process (*.pdb)
	 * 
	 * NOTE: The archive MUST be compressed using BZIP2 or GZIP only. No other compression schemes are supported.
	 * 
	 * TO PROCESS A MULTI FILE ARCHIVE, the 'doARCHIVE=TRUE' flag must be set in the Input Parameters File.
	 * 
	 * @return original_PDB_coordinates_AA
	 */
	public Matrix get_Matrix_of_Original_PDB_Coords_from_TAR_Archive()
	{
		int ROWS_AA = (number_of_atoms_AA * 3);
		List<TarArchiveEntry> entries = zio.read_Entries_In_TAR_Archive(directory + archive_name);
		int frames = 0;
		for (TarArchiveEntry tae : entries)
			{
				if (tae.isFile()) frames++;
			}
		PDB_file_names = new ArrayList<>(frames);
		original_PDB_coordinates_AA = new Matrix(ROWS_AA, frames);

		List<Atom> atoms = null;
		TarArchiveInputStream tarInput = null;
		BufferedReader reader = null;

		String line;
		List<String> lines = new ArrayList<String>();

		try
			{
				// Define the input streams based on archive type: tar.bz2, tar.gz, or just .tar
				if (archive_name.endsWith("bz2")) tarInput = new TarArchiveInputStream(new BZip2CompressorInputStream(new FileInputStream(directory + archive_name)));
				else if (archive_name.endsWith("gz")) tarInput = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(directory + archive_name)));
				else
					tarInput = new TarArchiveInputStream(new FileInputStream(directory + archive_name));

				TarArchiveEntry currentEntry;
				int index = 0;
				while ((currentEntry = tarInput.getNextTarEntry()) != null)
					{
						if (currentEntry.isDirectory()) continue;

						PDB_file_names.add(currentEntry.getName());
						reader = new BufferedReader(new InputStreamReader(tarInput));
						// System.out.println("Processing " + currentEntry.getName());
						while ((line = reader.readLine()) != null && line.length() >= 1)
							{
								lines.add(line);
								// System.out.println(line);
							}
						atoms = parse_PDB_as_Lines_List(lines);

						if (Input_Parameters.doParityCheck)
							{
								// Verify that each PDB file read is equivalent to the reference PDB file:
								// Here, minimum parity means same atom type, same residue type, and same residue number.
								// Atom numbers may be different... This is needed for pooling mutants with INDELs.
								for (int k = 0; k < atoms.size(); k++)
									{
										Atom refAtom = ref_atoms_AA.get(k);
										Atom test = atoms.get(k);
										boolean ok = refAtom.minimumParity(refAtom, test);
										if (!ok)
											{
												System.err.println("The PDB file " + currentEntry.getName() + " differes from reference at position " + (k + 1));
												System.err.println(refAtom.toVeryShortString());
												System.err.println(test.toVeryShortString());
											}
									}
							}

						Matrix cv_aa = get_Cartesian_PDB_Coords(atoms);
						original_PDB_coordinates_AA.setMatrix(0, ROWS_AA - 1, index, index, cv_aa);
						lines.clear();
						index++;
					}
				tarInput.close();
				reader.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return original_PDB_coordinates_AA;
	}

	/**
	 * Returns the matrix of the all atom coordinates from all PDB files in specified ZIP archive
	 * 
	 * NOTE: The archive MUST be a SINGLE DIRECTORY containing all the PDB files to process (*.pdb)
	 * 
	 * TO PROCESS AN ARCHIVE, the 'doARCHIVE=TRUE' flag must be set in the Input Parameters File.
	 * 
	 * @return original_PDB_coordinates_AA
	 */
	public Matrix get_Matrix_of_Original_PDB_Coords_from_ZIP_Archive()
	{
		int ROWS_AA = (number_of_atoms_AA * 3);
		List<Atom> atoms = null;
		ZipFile zipFile = null;
		InputStream zip = null;
		BufferedReader reader = null;
		List<String> lines = new ArrayList<String>();
		String line;

		try
			{
				zipFile = new ZipFile(directory + archive_name);
				Enumeration<? extends ZipEntry> entries = zipFile.entries();
				int frames = 0;
				List<ZipEntry> zfe = zio.read_Entries_In_ZIP_Archive(directory + archive_name);
				for (ZipEntry z : zfe)
					{
						if (!(z.isDirectory())) frames++;
					}
				PDB_file_names = new ArrayList<>(frames);
				original_PDB_coordinates_AA = new Matrix(ROWS_AA, frames);

				ZipEntry entry;
				int index = 0;
				while (entries.hasMoreElements())
					{
						entry = entries.nextElement();

						if (entry.isDirectory()) continue;

						PDB_file_names.add(entry.getName());
						zip = zipFile.getInputStream(entry);
						reader = new BufferedReader(new InputStreamReader(zip, "UTF-8"));
						// System.out.println("Processing " + entry.getName());
						while ((line = reader.readLine()) != null && line.length() >= 1)
							{
								lines.add(line);
								// System.out.println(line);
							}
						atoms = parse_PDB_as_Lines_List(lines);

						if (Input_Parameters.doParityCheck)
							{
								// Verify that each PDB file read is equivalent to the reference PDB file:
								// Here, minimum parity means same atom type, same residue type, and same residue number.
								// Atom numbers may be different... This is needed for pooling mutants with INDELs.
								for (int k = 0; k < atoms.size(); k++)
									{
										Atom refAtom = ref_atoms_AA.get(k);
										Atom test = atoms.get(k);
										boolean ok = refAtom.minimumParity(refAtom, test);
										if (!ok)
											{
												System.err.println("The PDB file " + entry.getName() + " differes from reference at position " + (k + 1));
												System.err.println(refAtom.toVeryShortString());
												System.err.println(test.toVeryShortString());
											}
									}
							}

						Matrix cv_aa = get_Cartesian_PDB_Coords(atoms);
						original_PDB_coordinates_AA.setMatrix(0, ROWS_AA - 1, index, index, cv_aa);
						lines.clear();
						index++;
					}
				zipFile.close();
				zip.close();
				reader.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return original_PDB_coordinates_AA;
	}

	/* ************************************************* GET REFERENCE SUBSET METHODS **************************************************************************** */

	public List<Atom> get_Reference_Atom_Subset_List(List<Atom> ref_atoms, List<Residue_ID_Pair> id_pairs)
	// Returns a List of atoms that is a subset, based on specified residue list
	// Creates a List of the subset atom numbers from the specified residue lists and an adjusted atom numbers list to access the Matrix of PDB coordinates
	{
		int size = ref_atoms.size();
		ref_subset_atoms = new ArrayList<Atom>(size);
		SS_atom_numbers = new ArrayList<Integer>(size);
		SS_atom_numbers_original = new ArrayList<Integer>(size);

		for (Atom a : ref_atoms)
			{
				String id = a.getChainID(); // Get the chain ID of the atom:
				int val = a.getRes_number(); // Get the residue number of the atom:
				Residue_ID_Pair ref_id_pair = new Residue_ID_Pair(id, val);

				for (Residue_ID_Pair rip : id_pairs)
					{
						if (ref_id_pair.equals(rip))
							{
								Atom ra = new Atom(a); // Create new atom to add to the subset.
								ref_subset_atoms.add(ra);
								SS_atom_numbers_original.add(ra.getAtom_number());
								System.out.println("\t\t" + ra.toString());
							}
					}

			}
		for (int number : SS_atom_numbers_original) // Need to correct atom numbering to be consistent with Java Array Numbering
			{
				int i = atom_numbers_read_original.indexOf(number); // Get the index of the atom number in the original list of atom numbers from PDB file.
				SS_atom_numbers.add(i);
				// System.out.println("\t\tIndex of atom number in original list of atom numbers: " + i);
			}
		return ref_subset_atoms;
	}

	/* ************************************************* READ RESIDUE and ATOM LISTS METHODS **************************************************************************** */

	public List<Residue_ID_Pair> read_Residue_List(String res_list)
	// Reads a file containing a residue list to access records in a PDB file;
	// If no Chain ID is specified, the assumption is that the PDB file contains a single chain of residues...
	// In this case, the default chain ID of "A" supplied by the PDB Parser is used
	// Returns a residue list with residue numbering adjusted to access the X,Y,Z block-packing in the coordinates matrix.
	{
		try
			{
				// System.out.println("\t\tReading the Residue List: ");

				residue_list_orig = new ArrayList<Residue_ID_Pair>(1000); // preserves residue numbering for accessing the PDB file.
				residue_list = new ArrayList<Residue_ID_Pair>(1000); // adjusts residue numbering to access X,Y,Z packing in the coordinates matrix.

				File residues = new File(directory + res_list);
				BufferedReader residue_reader = new BufferedReader(new FileReader(residues));
				String line, Chain_ID;
				while ((line = residue_reader.readLine()) != null)
					{
						sToken = new StringTokenizer(line);
						if (sToken.countTokens() < 2)
							{
								Chain_ID = "A";
							}
						else
							Chain_ID = sToken.nextToken();
						// chain_ids.add(Chain_ID);
						Integer res = Integer.parseInt(sToken.nextToken());
						// residue_list_original.add(res);
						if (sToken.hasMoreTokens())
							{
								System.err.println("ERROR: Residue list must only have two columns: ChainID and Res#");
								System.err.println("Please check file format. Terminating program execution.");
								System.exit(0);
							}
						Residue_ID_Pair ID_Pair = new Residue_ID_Pair(Chain_ID, res);
						residue_list_orig.add(ID_Pair);

						// System.out.println("\t\t\t" + ID_Pair.toString());

						boolean test = false;
						for (Residue_ID_Pair rip : residue_ID_pairs_read)
							{
								if (ID_Pair.equals(rip))
									{
										// System.out.println(rip.toString());
										test = true;
									}
							}
						if (!test)
							{
								System.err.println("ERROR: Requested Chain_ID + Residue_Number Pair DOES NOT EXIST in the Reference PDB File: " + ID_Pair.toString());
								System.err.println("Terminating program execution.");
								System.exit(0);
							}
					}
				residue_reader.close();
				number_of_residues_SS = residue_list_orig.size();
			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not read the residue list file: " + directory + res_list);
				System.err.println("Terminating program execution.");
				io.printStackTrace();
				System.exit(0);
			}
		return residue_list_orig;
	}

	public void read_Atoms_List(String atomList)
	// Reads a file containing a list of atom numbers.
	// If no Chain ID are specified, the assumption is that the PDB file contains a single chain of residues...
	// In this case, the default chain ID of "A" supplied by the PDB Parser is used
	{
		try
			{
				File atm_list = new File(directory + atomList);
				BufferedReader atom_list_reader = new BufferedReader(new FileReader(atm_list));

				atom_list = new ArrayList<Integer>(100);
				atom_list_original = new ArrayList<Integer>(100);
				chain_ids0 = new ArrayList<String>(100);

				String line = null, element;
				while ((line = atom_list_reader.readLine()) != null)
					{
						sToken = new StringTokenizer(line);
						// This is to handle the case of a single chain PDB file using the default chain ID of "A"
						if (sToken.countTokens() == 1)
							{
								String chain_ID = "A";
								chain_ids0.add(chain_ID);
								element = sToken.nextToken();
								int atm = Integer.parseInt(element);
								atom_list_original.add(atm);
								Atom_ID_Pair pair = new Atom_ID_Pair(chain_ID, atm);
								int index = atom_ID_pairs_read.indexOf(pair);

								if (index == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID + atm);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list.add(index);
								// System.out.println("\t\t\t" + pair.toString());
							}
						else if (sToken.countTokens() == 2)
							{
								String chain_ID = sToken.nextToken();
								chain_ids0.add(chain_ID);
								element = sToken.nextToken();
								int atm = Integer.parseInt(element);
								atom_list_original.add(atm);
								Atom_ID_Pair pair = new Atom_ID_Pair(chain_ID, atm);
								int index = atom_ID_pairs_read.indexOf(pair);

								if (index == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID + atm);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list.add(index);
								// System.out.println("\t\t" + pair1.toString() + "\t" + pair2.toString());
							}
						else
							{
								System.err.println("ERROR: Atom list should have either 1 or 2 columns: (Atom#) OR (Atom#) and (ChainID)");
							}
					}
				atom_list_reader.close();
				number_of_atoms = atom_list.size();
			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not read the residue list file: " + directory + atomList);
				io.printStackTrace();
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
	}

	public void read_Atom_Pairs(String pairs)
	// Reads a file containing an atom pairs list for distance pair processing.
	// If no Chain ID are specified, the assumption is that the PDB file contains a single chain of residues...
	// In this case, the default chain ID of "A" supplied by the PDB Parser is used
	{
		try
			{
				// System.out.println("\t\tReading the Atom Pairs List: ");

				File atom_pairs = new File(directory + pairs);
				BufferedReader atom_pair_reader = new BufferedReader(new FileReader(atom_pairs));

				atom_list1 = new ArrayList<Integer>(100);
				atom_list1_original = new ArrayList<Integer>(100);
				atom_list2 = new ArrayList<Integer>(100);
				atom_list2_original = new ArrayList<Integer>(100);
				chain_ids1 = new ArrayList<String>(100);
				chain_ids2 = new ArrayList<String>(100);

				String line = null, element;
				while ((line = atom_pair_reader.readLine()) != null)
					{
						sToken = new StringTokenizer(line);
						// This is to handle the case of a single chain PDB file using the default chain ID of "A"
						if (sToken.countTokens() == 2)
							{
								String chain_ID1 = "A";
								String chain_ID2 = "A";
								chain_ids1.add(chain_ID1);
								chain_ids2.add(chain_ID2);

								element = sToken.nextToken();
								int atm1 = Integer.parseInt(element);
								atom_list1_original.add(atm1);
								Atom_ID_Pair pair1 = new Atom_ID_Pair(chain_ID1, atm1);
								int res_index1 = atom_ID_pairs_read.indexOf(pair1);

								if (res_index1 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID1 + atm1);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list1.add(res_index1);

								element = sToken.nextToken();
								if (sToken.hasMoreTokens()) System.err.println("ERROR, Too many colums in file.");
								int atm2 = Integer.parseInt(element);
								atom_list2_original.add(atm2);
								Atom_ID_Pair pair2 = new Atom_ID_Pair(chain_ID2, atm2);
								int res_index2 = atom_ID_pairs_read.indexOf(pair2);

								if (res_index2 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID2 + atm2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list2.add(res_index2);
								// System.out.println("\t\t\t" + pair1.toString() + "\t" + pair2.toString());
							}
						else if (sToken.countTokens() == 4)
							{
								// Reading the first list of pairs: Column 1: chain IDs and Column 2: atom numbers
								String chain_ID1 = sToken.nextToken();
								chain_ids1.add(chain_ID1);
								element = sToken.nextToken();
								int atm1 = Integer.parseInt(element);
								atom_list1_original.add(atm1);
								Atom_ID_Pair pair1 = new Atom_ID_Pair(chain_ID1, atm1);
								int res_index1 = atom_ID_pairs_read.indexOf(pair1);

								if (res_index1 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID1 + atm1);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list1.add(res_index1);

								// Reading the second pair of atoms: Column 3: chain IDs and Column 4: atom numbers
								String chain_ID2 = sToken.nextToken();
								chain_ids2.add(chain_ID2);
								element = sToken.nextToken();
								int atm2 = Integer.parseInt(element);
								atom_list2_original.add(atm2);
								Atom_ID_Pair pair2 = new Atom_ID_Pair(chain_ID2, atm2);
								int res_index2 = atom_ID_pairs_read.indexOf(pair2);

								if (res_index2 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID2 + atm2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list2.add(res_index2);
								// System.out.println("\t\t" + pair1.toString() + "\t" + pair2.toString());

							}
						else
							{
								System.err.println("ERROR: Atom pair list should only have 4 columns: 2 Res# and 2 ChainID");
							}
					}
				atom_pair_reader.close();
				number_of_atom_pairs = atom_list1.size();
			}
		catch (IOException io)
			{
				System.err.println("IOException thrown. Could not read the residue list file: " + directory + pairs);
				io.printStackTrace();
				System.err.println("Terminating program execution.");
				System.exit(0);
			}
	}

	/* ************************************************* COORDS FROM LIST METHOD **************************************************************************** */

	public Matrix get_Cartesian_PDB_Coords(List<Atom> data)
	// Returns a matrix COLUMN VECTOR of Cartesian coordinates with X,Y,Z block packing, based on the given list of atoms (a Single Conformation)
	{
		int i = 0;
		int number_of_atoms = data.size();

		Matrix col_vector = new Matrix(3 * number_of_atoms, 1);

		for (Atom a : data)
			{
				double x = a.getX();
				double y = a.getY();
				double z = a.getZ();

				col_vector.set(i, 0, x);
				col_vector.set((i + number_of_atoms), 0, y);
				col_vector.set((i + (2 * number_of_atoms)), 0, z);
				i++;
			}
		return col_vector;
	}

	/* ******************************************** GET RESIDUE ATOMs as LISTs of LISTs METHODs ***************************************************************** */

	public List<List<Atom>> get_Residue_Atoms_List(List<Atom> atoms, List<Residue_ID_Pair> id_pairs)
	{
		Residue_Atoms_List = new ArrayList<List<Atom>>(1000);
		for (int i = 0; i < id_pairs.size(); i++)
			{
				String chainID = id_pairs.get(i).getChain_ID();
				int residue_Number = id_pairs.get(i).getResidue_Number();
				List<Atom> residue_atoms = new ArrayList<Atom>(100);

				// System.out.println("\n\t Getting Residue Atoms:");
				// System.out.println("\t Chain ID = " + chainID);
				// System.out.println("\t Residue Number = " + residue_Number + "\n");

				for (Atom a : atoms)
					{
						if (a.chainID.equals(chainID))
							{
								if (a.res_number == residue_Number)
									{
										Atom ra = new Atom(a);
										residue_atoms.add(ra);
										// System.out.println("\t" + ra.toString());
									}
							}
					}
				Residue_Atoms_List.add(residue_atoms);
			}
		return Residue_Atoms_List;
	}

	/* ************************************************* GET RESIDUE COORDINATES METHODS *********************************************************************** */

	public List<Matrix> get_Residue_Coords(Matrix data, List<Integer> atoms_per_residue)
	{
		int offset = 0;
		int number_of_conformations = data.getColumnDimension();
		int num_atoms = data.getRowDimension() / 3;
		// System.out.println("Total Number of atoms: " + num_atoms);
		// System.out.println("residues size: " + atoms_per_residue.size());

		Residue_Coords_List = new ArrayList<Matrix>(atoms_per_residue.size());
		for (int i = 0; i < atoms_per_residue.size(); i++)
			{
				int number_Of_Atoms_In_Residue = atoms_per_residue.get(i);
				// System.out.println("Residue: " + residue_list_orig.get(i));
				// System.out.println("Number of atoms: " + number_Of_Atoms_In_Residue);

				Matrix Residue_coordinates = new Matrix(3 * number_Of_Atoms_In_Residue, number_of_conformations);

				for (int j = 0; j < number_Of_Atoms_In_Residue; j++)
					{
						for (int k = 0; k < number_of_conformations; k++)
							{
								double element_X = data.get((j + offset), k);
								double element_Y = data.get(((j + offset) + (num_atoms)), k);
								double element_Z = data.get(((j + offset) + (2 * num_atoms)), k);

								Residue_coordinates.set(j, k, element_X);
								Residue_coordinates.set((j + number_Of_Atoms_In_Residue), k, element_Y);
								Residue_coordinates.set((j + (2 * number_Of_Atoms_In_Residue)), k, element_Z);
							}
					}
				Residue_Coords_List.add(Residue_coordinates);
				offset += number_Of_Atoms_In_Residue;
				// System.out.println("offset " + offset);
			}
		return Residue_Coords_List;
	}

	public List<Matrix> get_Residue_Coords_HA(Matrix Coords, List<Integer> heavy_atoms_per_residue)
	{
		int offset = 0;
		int number_of_conformations = Coords.getColumnDimension();

		Matrix data = get_Subset_Coords(Coords, heavy_atom_numbers_read);
		int num_atoms = data.getRowDimension() / 3;
		// System.out.println("Total Number of atoms: " + num_atoms);
		// System.out.println("residues_read size: " + residues_read.size());

		Residue_Coords_HA_List = new ArrayList<Matrix>(residue_numbers_read_original.size());
		for (int i = 0; i < residue_numbers_read_original.size(); i++)
			{
				int number_Of_Atoms_In_Residue = heavy_atoms_per_residue.get(i);
				Matrix Residue_coordinates = new Matrix(3 * number_Of_Atoms_In_Residue, number_of_conformations);
				// System.out.println("Residue: " + residues_read.get(i));
				// System.out.println("Number of atoms: " + number_Of_Atoms_In_Residue);

				for (int j = 0; j < number_Of_Atoms_In_Residue; j++)
					{
						for (int k = 0; k < number_of_conformations; k++)
							{
								double element_X = data.get((j + offset), k);
								double element_Y = data.get(((j + offset) + (num_atoms)), k);
								double element_Z = data.get(((j + offset) + (2 * num_atoms)), k);

								Residue_coordinates.set(j, k, element_X);
								Residue_coordinates.set((j + number_Of_Atoms_In_Residue), k, element_Y);
								Residue_coordinates.set((j + (2 * number_Of_Atoms_In_Residue)), k, element_Z);
							}
					}
				Residue_Coords_HA_List.add(Residue_coordinates);
				offset += number_Of_Atoms_In_Residue;
			}
		return Residue_Coords_HA_List;
	}

	/* ************************************************* GET COORDINATE SUBSET METHODS **************************************************************************** */

	public Matrix get_Subset_Coords(Matrix coords, List<Integer> atom_numbers)
	// Returns a matrix with a subset of coordinates from the all atom matrix of PDB coordinates using the specified list of atom numbers, with X,Y,Z block packing;
	// For this method to work, the list of atom numbers must be adjusted to access the coordinates matrix, using Java Array Numbering (starting with 0).
	// This means that the atom numbers must be obtained from the `numbers_original` list, and then adjusted for both **continuity** and **JAN**.
	{
		int number_of_conformations = coords.getColumnDimension();
		int num_atoms = coords.getRowDimension() / 3;
		int num_atoms_SS = atom_numbers.size();

		// System.out.println("Number of atoms in Subset :" + num_atoms_SS);

		Matrix subset_PDB_coordinates = new Matrix(3 * num_atoms_SS, number_of_conformations);
		for (int i = 0; i < num_atoms_SS; i++) // Array Indices must reference ZERO as first element (Java Array Numbering)
			{
				for (int j = 0; j < number_of_conformations; j++)
					{
						double element_X = coords.get((atom_numbers.get(i)), j);
						double element_Y = coords.get((atom_numbers.get(i) + (num_atoms)), j);
						double element_Z = coords.get((atom_numbers.get(i) + (2 * num_atoms)), j);

						subset_PDB_coordinates.set(i, j, element_X);
						subset_PDB_coordinates.set((i + num_atoms_SS), j, element_Y);
						subset_PDB_coordinates.set((i + (2 * num_atoms_SS)), j, element_Z);
					}
			}
		return subset_PDB_coordinates;
	}

	/* ************************************************* EDIT PDB FILE METHODS **************************************************************************** */

	public Vector<Atom> Add_Missing_Chain_ID(Vector<Atom> atoms, String ID)
	{
		// System.out.println("\t Adding Missing Chain ID: ");

		int size = atoms.size();
		Vector<Atom> atoms_edited = new Vector<Atom>(size, 10);

		for (Atom a : atoms)
			{
				Atom ea = new Atom(a);
				if (ea.chainID.isEmpty())
					{
						ea.setChainID(ID);
						// System.out.println("Chain ID = " + a.chainID);
					}
				atoms_edited.add(ea);
				// System.out.println(ea.toString());
			}
		return atoms_edited;
	}

	public List<Atom> edit_B_Factors_with_RMSFs(List<Atom> atoms, List<Double> rmsfs)
	{
		// System.out.println("\t Editing Atom B-Factors: ");

		int size = atoms.size();
		List<Atom> editedAtoms = new ArrayList<Atom>(size);

		List<Double> sorted_res_rmsds = new ArrayList<>(rmsfs.size());
		sorted_res_rmsds.addAll(rmsfs);
		Collections.sort(sorted_res_rmsds, Collections.reverseOrder());

		double max_rmsd = sorted_res_rmsds.get(0);
		double min_rmsd = sorted_res_rmsds.get(rmsfs.size() - 1);
		if (min_rmsd < FLOOR) min_rmsd = FLOOR;
		double log_RR_max = Math.log10(max_rmsd);
		double log_RR_min = Math.log10(min_rmsd);
		double delta_x = ((log_RR_max) - (log_RR_min));
		double slope = (delta_y / delta_x);
		double y_min = (slope * log_RR_min);

		int i = 0;
		for (Atom a : atoms)
			{
				Atom ea = new Atom(a);

				double bff = (rmsfs.get(i));
				double log_bff = Math.log10(bff);
				double bf = ((slope * log_bff) - y_min);

				ea.setB_factor(bf);
				editedAtoms.add(ea);

				// System.out.println("Original atom: ");
				// System.out.println(atoms.get(i).toString());
				// System.out.println("Edited atom: ");
				// System.out.println(ea.toString());

				i++;
			}
		return editedAtoms;
	}


	/* ************************************************* SETTERS **************************************************************************** */

	public void setOut_dir(String out_dir)
	{
		this.out_dir = out_dir;
	}

	/* ************************************************* GETTERS **************************************************************************** */

	public int getNumber_of_atoms()
	{
		return number_of_atoms;
	}

	public int getNumber_of_atoms_AA()
	{
		return number_of_atoms_AA;
	}

	public int getNumber_of_residues()
	{
		return number_of_residues;
	}

	public int getNumber_of_residues_SS()
	{
		return number_of_residues_SS;
	}

	// ------------------------------------------------------------------------- //

	public Matrix getOriginal_PDB_coordinates_AA()
	{
		return original_PDB_coordinates_AA;
	}

	// ------------------------------------------------------------------------- //

	public Matrix getReference_PDB_coordinates_AA()
	{
		return reference_PDB_coordinates_AA;
	}

	public Matrix getReference_PDB_coordinates_BB()
	{
		return reference_PDB_coordinates_BB;
	}

	public Matrix getReference_PDB_coordinates_CA()
	{
		return reference_PDB_coordinates_CA;
	}

	public Matrix getReference_PDB_coordinates_HA()
	{
		return reference_PDB_coordinates_HA;
	}

	// ------------------------------------------------------------------------- //

	public List<Matrix> getReference_Residue_Coords_List()
	{
		return Reference_Residue_Coords_List;
	}

	public List<Matrix> getResidue_Coords_List()
	{
		return Residue_Coords_List;
	}

	// ------------------------------------------------------------------------- //

	public Matrix getReference_subset_PDB_coordinates_AA()
	{
		return reference_subset_PDB_coordinates_AA;
	}

	public Matrix getReference_subset_PDB_coordinates_BB()
	{
		return reference_subset_PDB_coordinates_BB;
	}

	public Matrix getReference_subset_PDB_coordinates_CA()
	{
		return reference_subset_PDB_coordinates_CA;
	}

	public Matrix getReference_subset_PDB_coordinates_HA()
	{
		return reference_subset_PDB_coordinates_HA;
	}

	// ------------------------------------------------------------------------- //

	public Matrix getSubset_PDB_coordinates_AA()
	{
		return subset_PDB_coordinates_AA;
	}

	public Matrix getSubset_PDB_coordinates_BB()
	{
		return subset_PDB_coordinates_BB;
	}

	public Matrix getSubset_PDB_coordinates_CA()
	{
		return subset_PDB_coordinates_CA;
	}

	public Matrix getSubset_PDB_coordinates_HA()
	{
		return subset_PDB_coordinates_HA;
	}

	// ------------------------------------------------------------------------- //

	public List<Atom> get_ref_Atoms_AA()
	{
		return ref_atoms_AA;
	}

	public List<Atom> get_ref_Atoms_BB()
	{
		return ref_atoms_BB;
	}

	public List<Atom> get_ref_Atoms_CA()
	{
		return ref_atoms_CA;
	}

	public List<Atom> get_ref_Atoms_HA()
	{
		return ref_atoms_HA;
	}

	// ------------------------------------------------------------------------- //

	public List<Atom> getAtoms_AA()
	{
		return atoms_AA;
	}

	public List<Atom> getAtoms_Backbone()
	{
		return atoms_BB;
	}

	public List<Atom> getAtoms_CA()
	{
		return atoms_CA;
	}

	public List<Atom> getAtoms_HA()
	{
		return atoms_HA;
	}

	// ------------------------------------------------------------------------- //

	public List<Integer> getAtom_list()
	{
		return atom_list;
	}

	public List<Integer> getAtom_list_original()
	{
		return atom_list_original;
	}

	public List<Integer> getAtom_list1()
	{
		return atom_list1;
	}

	public List<Integer> getAtom_list1_original()
	{
		return atom_list1_original;
	}

	public List<Integer> getAtom_list2()
	{
		return atom_list2;
	}

	public List<Integer> getAtom_list2_original()
	{
		return atom_list2_original;
	}

	public List<Integer> get_Atoms_Read()
	{
		return atom_numbers_read;
	}

	public List<Integer> get_Heavy_Atoms_Read()
	{
		return heavy_atom_numbers_read;
	}

	public List<Integer> get_Alpha_Carbons_Read()
	{
		return alpha_carbon_numbers_read;
	}

	public List<Integer> get_Residue_Numbers_Read()
	{
		return residue_numbers_read_original;
	}

	public List<Integer> get_SS_Atom_Numbers()
	{
		return SS_atom_numbers;
	}

	public List<Integer> get_Number_of_Atoms_in_Residues()
	{
		return number_of_Atoms_in_Residues;
	}

	public List<Integer> getNumber_of_Heavy_Atoms_in_Residues()
	{
		return number_of_Heavy_Atoms_in_Residues;
	}

	// ------------------------------------------------------------------------- //

	public List<String> get_PDB_file_names()
	{
		return PDB_file_names;
	}

	public List<String> get_Chain_IDs()
	{
		return chain_ids;
	}

	public List<String> get_Chain_IDs_Read()
	{
		return chain_IDs_read;
	}

	public List<String> getChain_ids0()
	{
		return chain_ids0;
	}

	public List<String> getChain_ids1()
	{
		return chain_ids1;
	}

	public List<String> getChain_ids2()
	{
		return chain_ids2;
	}

	// ------------------------------------------------------------------------- //

	public List<Residue_ID_Pair> get_Residue_List_Original()
	{
		return residue_list_orig;
	}

	public List<Residue_ID_Pair> getResidue_ID_pairs_read()
	{
		return residue_ID_pairs_read;
	}

	public List<Atom_ID_Pair> getAtom_ID_pairs_read()
	{
		return atom_ID_pairs_read;
	}

	// ------------------------------------------------------------------------- //
}
