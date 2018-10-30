package jedi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;

public class JEDi_PDB_Processing
	{
		int number_of_atoms, number_of_atoms_HA, number_of_atom_pairs, number_of_residues, number_of_residues_SS, number_of_residue_pairs, number_of_modes;
		final double FLOOR = 1.00E-3, delta_y = 99;
		String out_dir, directory, name, description, pdb_ref_file, rl_SS, rl_Cartesian_SS, rl_Cartesian_Residue_SS, rl_RGC_SS, atom_pairs_SS;
		List<Integer> atom_list1, atom_list1_original, atom_list2, atom_list2_original;
		List<Integer> atoms_read, residues_read, residue_list, residue_list_original, number_of_Atoms_in_Residues, number_of_Heavy_Atoms_in_Residues;
		List<String> chain_ids, chain_ids1, chain_ids2, chain_ids_read, residue_ID_pairs_read, atom_ID_pairs_read, PDB_file_names;
		Matrix reference_PDB_coordinates, reference_subset_PDB_coordinates_AA, original_PDB_coordinates_AA, subset_PDB_coordinates_AA;
		Matrix reference_PDB_coordinates_CA, reference_subset_PDB_coordinates_CA, original_PDB_coordinates_CA, subset_PDB_coordinates_CA;
		Matrix reference_PDB_coordinates_HA, reference_subset_PDB_coordinates_HA, original_PDB_coordinates_HA, subset_PDB_coordinates_HA;
		Matrix reference_RGC_coordinates, reference_subset_RGC_coordinates, original_RGC_coordinates, subset_Hierarchical_Coordinates;
		StringTokenizer sToken;
		Vector<Atom> atoms_AA, ref_atoms, atoms_CA, ref_atoms_CA, atoms_HA, ref_atoms_HA;
		Vector<Atom> ref_subset_atoms;
		List<Matrix> Residue_Coords_List, Residue_Coords_HA_List, Reference_Residue_Coords_List;
		PDB_IO pio;

		/* ************************************************* CONSTRUCTORS **************************************************************************** */

		public JEDi_PDB_Processing(String directory, String pdb_ref_file)
			{
				super();
				this.directory = directory;
				this.pdb_ref_file = pdb_ref_file;

				residue_list_original = new ArrayList<Integer>(); // preserves residue numbering for accessing the PDB file.
				residue_list = new ArrayList<Integer>(); // adjusts residue numbering to access X,Y,Z packing in the coordinates matrix.
				chain_ids = new ArrayList<String>(); // list of all chain identifiers specified: for Multi PDBs
				residue_ID_pairs_read = new ArrayList<String>(); // list of all chain---residue pairs specified: for Multi PDBs
				Reference_Residue_Coords_List = new ArrayList<Matrix>();
				ref_atoms = new Vector<Atom>();
				ref_atoms_CA = new Vector<Atom>();
				ref_atoms_HA = new Vector<Atom>();
			}

		public JEDi_PDB_Processing(String directory, String pdb_ref_file, String res_list)
			{
				super();
				this.directory = directory;
				this.pdb_ref_file = pdb_ref_file;
				this.rl_SS = res_list;

				residue_list_original = new ArrayList<Integer>(); // preserves residue numbering for accessing the PDB file.
				residue_list = new ArrayList<Integer>(); // adjusts residue numbering to access X,Y,Z packing in the coordinates matrix.
				chain_ids = new ArrayList<String>(); // list of all chain identifiers specified: for Multi PDBs
				residue_ID_pairs_read = new ArrayList<String>(); // list of all chain---residue pairs specified: for Multi PDBs
				Reference_Residue_Coords_List = new ArrayList<Matrix>();
			}

		/* ************************************************* METHODS **************************************************************************** */

		public void read_Reference_PDB() // for use with the MULTI residue lists for sub-setting purposes
			{
				pio = new PDB_IO(directory, pdb_ref_file);
				ref_atoms = pio.Read_PDB();
				ref_atoms_CA = pio.get_Atoms_CA();
				ref_atoms_HA = pio.get_Heavy_Atoms();
				atoms_read = pio.get_Atoms_read();
				residues_read = pio.get_Residues_read();
				chain_ids_read = pio.get_Chain_ids_read();
				residue_ID_pairs_read = pio.getResidue_ID_pairs();
				atom_ID_pairs_read = pio.getAtom_ID_pairs();
				number_of_Atoms_in_Residues = pio.get_Number_of_atoms_in_Residues();
				number_of_Heavy_Atoms_in_Residues = pio.get_Number_of_Heavy_Atoms_in_Residues();
				reference_PDB_coordinates = get_Cartesian_PDB_Coords(ref_atoms);
				reference_PDB_coordinates_CA = get_Cartesian_PDB_Coords(ref_atoms_CA);
				reference_PDB_coordinates_HA = get_Cartesian_PDB_Coords(ref_atoms_HA);
				Reference_Residue_Coords_List = get_Residue_Coords(reference_PDB_coordinates, number_of_Atoms_in_Residues);
				number_of_atoms = ref_atoms.size();
				number_of_residues = ref_atoms_CA.size();
				number_of_atoms_HA = ref_atoms_HA.size();
			}
		// ------------------------------------------------------------------------------------------------------------------------------------

		public void read_PDB(String directory, String pdb_file) // returns all atoms in specified PDB file, for building original_PDB_coordinates
			{
				pio = new PDB_IO(directory, pdb_file);
				atoms_AA = pio.Read_PDB();
				atoms_CA = pio.get_Atoms_CA();
				atoms_HA = pio.get_Heavy_Atoms();
			}

		public Matrix get_Original_PDB_Coords()
			{
				File pdb_file = new File(directory);
				FilenameFilter filter = new PDB_Filter(".pdb");
				String[] ls = pdb_file.list(filter);
				Arrays.sort(ls);
				int number_of_conformations = ls.length;
				int ROWS_AA = (number_of_atoms * 3);
				int ROWS_CA = (number_of_residues * 3);
				int ROWS_NCO = (number_of_atoms_HA * 3);

				PDB_file_names = new ArrayList<>();
				original_PDB_coordinates_AA = new Matrix(ROWS_AA, number_of_conformations);
				original_PDB_coordinates_CA = new Matrix(ROWS_CA, number_of_conformations);
				original_PDB_coordinates_HA = new Matrix(ROWS_NCO, number_of_conformations);

				for (int i = 0; ls != null && i < ls.length;)
					for (i = 0; i < ls.length; i++)
						{
							System.out.println("PDB file: " + ls[i]);
							PDB_file_names.add(ls[i]);
							read_PDB(directory, ls[i]);
							Matrix cv_aa = get_Cartesian_PDB_Coords(atoms_AA);
							Matrix cv_ca = get_Cartesian_PDB_Coords(atoms_CA);
							Matrix cv_ha = get_Cartesian_PDB_Coords(atoms_HA);
							original_PDB_coordinates_AA.setMatrix(0, ROWS_AA - 1, i, i, cv_aa);
							original_PDB_coordinates_CA.setMatrix(0, ROWS_CA - 1, i, i, cv_ca);
							original_PDB_coordinates_HA.setMatrix(0, ROWS_NCO - 1, i, i, cv_ha);
							if ((i % 10 == 0))
								System.gc();
						}
				return original_PDB_coordinates_AA;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public Vector<Atom> get_Reference_Subset_Single(Vector<Atom> atms, List<Integer> residue_list_original) // returns a subset of atoms from specified
																												 // vector of atoms, based on residue list
			{
				ref_subset_atoms = new Vector<Atom>();
				for (Atom a : atms)
					{
						int val = a.res_number;
						if (residue_list_original.indexOf(val) > -1)
							ref_subset_atoms.add(a);
					}
				return ref_subset_atoms;
			}

		public Vector<Atom> get_Reference_Subset_Multi(Vector<Atom> atms, List<Integer> residue_list_original, List<String> chain_ids) // returns a subset of
																																		 // atoms from specified
																																		 // vector of atoms
			{
				ref_subset_atoms = new Vector<Atom>();
				for (Atom a : atms)
					{
						int val = a.res_number;
						String chain_id = a.chainID;
						if (residue_list_original.indexOf(val) > -1 && chain_ids.indexOf(chain_id) > -1)
							ref_subset_atoms.add(a);
					}
				return ref_subset_atoms;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public List<Integer> read_residue_list_Single(String res_list) // reads a file containing a residue list to access Single PDB numbering;
		// returns a residue list with residue numbering to access X,Y,Z packing in the coordinates matrix.
			{
				try
					{
						residue_list = new ArrayList<Integer>();
						File residues = new File(directory + res_list);
						BufferedReader residue_reader = new BufferedReader(new FileReader(residues));
						String line;
						while ((line = residue_reader.readLine()) != null)
							{
								sToken = new StringTokenizer(line);
								String r = sToken.nextToken();
								if (sToken.hasMoreTokens())
									{
										System.err.println("ERROR: Residue list for SINGLE chain PDB file must have a SINGLE column of numbers.");
										System.err.println("Terminating program execution.");
										System.exit(0);

									}
								int res = Integer.parseInt(r);
								residue_list_original.add(res);
								int res_index = residues_read.indexOf(res);
								if (res_index == -1)
									{
										System.err.println("ERROR: Requested Residue DOES NOT EXIST in the Reference PDB File: " + res);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								residue_list.add(res_index);
							}
						number_of_residues = residue_list.size();
						residue_reader.close();
					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the residue list file: " + directory + res_list);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}
				return residue_list;

			}

		public List<Integer> getResidue_list_original()
			{
				return residue_list_original;
			}

		public List<Integer> read_residue_list_Multi(String res_list_multi) // reads a file containing a residue list to access Multi PDB numbering;
		// returns a residue list with residue numbering to access X,Y,Z packing in the coordinates matrix.
			{
				try
					{
						residue_list = new ArrayList<Integer>();
						rl_SS = res_list_multi;
						File residues = new File(directory + rl_SS);
						BufferedReader residue_reader = new BufferedReader(new FileReader(residues));
						String line, Chain_ID;
						while ((line = residue_reader.readLine()) != null)
							{
								sToken = new StringTokenizer(line);
								if (sToken.countTokens() < 2)
									{
										System.err.println("ERROR: Multi Chain PDB residue file list must have 2 columns: Res# and ChainID");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								Chain_ID = sToken.nextToken();
								chain_ids.add(Chain_ID);
								int res = Integer.parseInt(sToken.nextToken());
								if (sToken.hasMoreTokens())
									{
										System.err.println("ERROR: MULTI chain PDB residue list must only have two columns: Res# and ChainID");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								residue_list_original.add(res);
								String ID_Pair = Chain_ID + res;
								int res_index = residue_ID_pairs_read.indexOf(ID_Pair);
								if (res_index == -1)
									{
										System.err.println("ERROR: Requested Chain_ID + Residue Pair DOES NOT EXIST in the Reference PDB File: " + Chain_ID + res);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								residue_list.add(res_index);
							}
						residue_reader.close();
						number_of_residues = residue_list.size();
					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the residue list file: " + directory + rl_SS);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}
				return residue_list;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public void read_atom_pairs_Single(String pairs) // reads a file containing atom pairs for pair distance processing, formatted for Single PDBs
			{
				try
					{
						atom_pairs_SS = pairs;
						File atm_pairs = new File(directory + atom_pairs_SS);
						BufferedReader atom_pair_reader = new BufferedReader(new FileReader(atm_pairs));

						atom_list1 = new ArrayList<Integer>();
						atom_list1_original = new ArrayList<Integer>();
						atom_list2 = new ArrayList<Integer>();
						atom_list2_original = new ArrayList<Integer>();

						String line, element;

						while ((line = atom_pair_reader.readLine()) != null)
							{

								sToken = new StringTokenizer(line);
								if (sToken.countTokens() < 2)
									{
										System.err.println("ERROR: Single Chain PDB atom pair list must have 2 columns.");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								element = sToken.nextToken();
								int res1 = Integer.parseInt(element);
								atom_list1_original.add(res1);
								int res1_index = atoms_read.indexOf(res1);
								atom_list1.add(res1_index);
								element = sToken.nextToken();
								if (sToken.hasMoreTokens())
									System.err.println("ERROR: Atom list for SINGLE chain PDB file must have ONLY 2 columns of numbers.");
								int res2 = Integer.parseInt(element);
								atom_list2_original.add(res2);
								int res2_index = atoms_read.indexOf(res2);
								atom_list2.add(res2_index);
								if (res1_index == -1 || res2_index == -1)
									{
										System.err.println("ERROR: Requested Atom Pair DOES NOT EXIST in the Reference PDB File: " + res1 + "\t" + res2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
							}
						atom_pair_reader.close();
						number_of_atom_pairs = atom_list1.size();
					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the atom pair file: " + directory + rl_SS);
						System.err.println("Terminating program execution.");
						io.printStackTrace();
						System.exit(0);
					}
			}

		public void read_atom_pairs_Multi(String pairs) // reads a file containing atom pairs for pair distance processing, formatted for Multi PDBs
			{
				try
					{
						rl_SS = pairs;
						File residues = new File(directory + rl_SS);
						BufferedReader residue_reader = new BufferedReader(new FileReader(residues));

						atom_list1 = new ArrayList<Integer>();
						atom_list1_original = new ArrayList<Integer>();
						atom_list2 = new ArrayList<Integer>();
						atom_list2_original = new ArrayList<Integer>();
						chain_ids1 = new ArrayList<String>();
						chain_ids2 = new ArrayList<String>();

						String line = null, element;
						while ((line = residue_reader.readLine()) != null)
							{
								// Reading the first list of residues: Column 1: chain IDs and Column 2: residue numbers

								sToken = new StringTokenizer(line);
								if (sToken.countTokens() < 4)
									{
										System.err.println("ERROR: Multi Chain PDB atom pair list must have 4 columns: 2 Res# and 2 ChainID!");
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								String chain_ID1 = sToken.nextToken();
								chain_ids1.add(chain_ID1);
								element = sToken.nextToken();
								int res1 = Integer.parseInt(element);
								atom_list1_original.add(res1);
								String ID_Pair1 = chain_ID1 + res1;
								int res_index1 = atom_ID_pairs_read.indexOf(ID_Pair1);
								if (res_index1 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID1 + res1);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list1.add(res_index1);

								// Reading the second list of residues: Column 3: chain IDs and Column 4: residue numbers

								String chain_ID2 = sToken.nextToken();
								chain_ids2.add(chain_ID2);
								element = sToken.nextToken();
								if (sToken.hasMoreTokens())
									System.err.println("ERROR: MULTI chain PDB atom pair list should only have 4 columns: 2 Res# and 2 ChainID");
								int res2 = Integer.parseInt(element);
								atom_list2_original.add(res2);
								String ID_Pair2 = chain_ID2 + res2;
								int res_index2 = atom_ID_pairs_read.indexOf(ID_Pair2);
								if (res_index2 == -1)
									{
										System.err.println("ERROR: Requested Atom DOES NOT EXIST in the Reference PDB File: " + chain_ID2 + res2);
										System.err.println("Terminating program execution.");
										System.exit(0);
									}
								atom_list2.add(res_index2);
							}
						residue_reader.close();
						number_of_atom_pairs = atom_list1.size();
					} catch (IOException io)
					{
						System.err.println("IOException thrown. Could not read the residue list file: " + directory + rl_SS);
						io.printStackTrace();
						System.err.println("Terminating program execution.");
						System.exit(0);
					}
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public List<String> getChain_ids1()
			{
				return chain_ids1;
			}

		public List<String> getChain_ids2()
			{
				return chain_ids2;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public Matrix get_Cartesian_PDB_Coords(Vector<Atom> data) // returns a matrix of Cartesian coordinates with X,Y,Z packing based on the given vector of
																	 // atoms
			{
				int i = 0, number_of_atoms = data.size();
				Matrix col_vector = new Matrix(3 * number_of_atoms, 1);

				for (Atom a : data)
					{
						col_vector.set(i, 0, a.x);
						col_vector.set((i + number_of_atoms), 0, a.y);
						col_vector.set((i + (2 * number_of_atoms)), 0, a.z);
						i++;
					}
				return col_vector;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public List<Matrix> get_Residue_Coords(Matrix data, List<Integer> atoms_per_residue) // data should be the all atom original PDB coordinates; uses
																							 // residues_read to do all residues
			{
				int offset = 0;
				int number_of_conformations = data.getColumnDimension();
				int num_atoms = data.getRowDimension() / 3;
				Residue_Coords_List = new ArrayList<Matrix>();
				for (int i = 0; i < residues_read.size(); i++)
					{
						int number_Of_Atoms_In_Residue = atoms_per_residue.get(i);
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
						Residue_Coords_List.add(Residue_coordinates);
						offset += number_Of_Atoms_In_Residue;
					}
				return Residue_Coords_List;
			}

		public List<Matrix> get_Residue_Coords_HA(Matrix data, List<Integer> atoms_per_residue) // data should be the NCO PDB coordinates; uses residues_read to
																								 // do all residues
			{
				int offset = 0;
				int number_of_conformations = data.getColumnDimension();
				int num_atoms = data.getRowDimension() / 3;
				// System.out.println("Total Number of atoms: " + num_atoms);
				// System.out.println("residues_read size: " + residues_read.size());

				Residue_Coords_HA_List = new ArrayList<Matrix>();
				for (int i = 0; i < residues_read.size(); i++)
					{
						int number_Of_Atoms_In_Residue = atoms_per_residue.get(i);
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

		// ------------------------------------------------------------------------------------------------------------------------------------

		public Matrix get_subset_Cartesian_Coords_CA(Matrix CA_coords, List<Integer> residues) // returns a matrix of Cartesian sub-set coordinates with X,Y,Z
																								 // packing; expects the adjusted residue lists.
			{
				int number_of_conformations = CA_coords.getColumnDimension();
				number_of_residues = CA_coords.getRowDimension() / 3;
				number_of_residues_SS = residues.size();
				subset_PDB_coordinates_CA = new Matrix(3 * number_of_residues_SS, number_of_conformations);
				for (int i = 0; i < number_of_residues_SS; i++) // array indices must reference ZERO as first element (Java Array Numbering)
					{
						for (int j = 0; j < number_of_conformations; j++)
							{
								double element_X = CA_coords.get((residues.get(i)), j);
								double element_Y = CA_coords.get((residues.get(i) + (number_of_residues)), j);
								double element_Z = CA_coords.get((residues.get(i) + (2 * number_of_residues)), j);

								subset_PDB_coordinates_CA.set(i, j, element_X);
								subset_PDB_coordinates_CA.set((i + number_of_residues_SS), j, element_Y);
								subset_PDB_coordinates_CA.set((i + (2 * number_of_residues_SS)), j, element_Z);
							}
					}
				return subset_PDB_coordinates_CA;
			}

		public Matrix get_subset_Cartesian_Coords_HA(List<Matrix> res_coords, List<Integer> residues, List<Integer> number_of_Atoms_in_Residues)
			{
				int number_of_conformations = res_coords.get(0).getColumnDimension();
				int number_of_residues_SS = residues.size();
				int number_atoms_SS = 0;
				for (int res : residues)
					{
						System.out.println("Residue: " + res);
						int atom_count = number_of_Atoms_in_Residues.get(res);
						number_atoms_SS += atom_count;
						System.out.println("Residue atoms: " + atom_count);
						System.out.println("Total atoms: " + number_atoms_SS);
					}

				int ROWS_SS = (number_atoms_SS * 3);
				int COLS_SS = number_of_conformations;
				System.out.println("Rows in SS matrix: " + ROWS_SS);

				subset_PDB_coordinates_HA = new Matrix(ROWS_SS, COLS_SS);

				int offset = 0;
				for (int i = 0; i < number_of_residues_SS; i++) // array indices must reference ZERO as first element (Java Array Numbering)
					{
						int residue = residue_list.get(i);
						int count = number_of_Atoms_in_Residues.get(residue);
						// System.out.println("Count = " + count);
						// System.out.println("Offset = " + offset);
						Matrix coords = res_coords.get(residue);
						// System.out.println("Coordinates for Residue :" + residue);
						// System.out.println("Rows/3 :" + coords.getRowDimension() / 3);
						Matrix X = coords.getMatrix(0, count - 1, 0, number_of_conformations - 1);
						Matrix Y = coords.getMatrix(count, (2 * count - 1), 0, number_of_conformations - 1);
						Matrix Z = coords.getMatrix(2 * count, (3 * count - 1), 0, number_of_conformations - 1);
						int X1 = offset;
						int X2 = count + offset - 1;
						int Y1 = number_atoms_SS + offset;
						int Y2 = (count + number_atoms_SS + offset - 1);
						int Z1 = 2 * number_atoms_SS + offset;
						int Z2 = (count + 2 * number_atoms_SS + offset - 1);
						// System.out.println("X1,X2,Y1,Y2,Z1,Z2 " + X1 + "\t" + X2 + "\t" + Y1 + "\t" + Y2 + "\t" + Z1 + "\t" + Z2);
						subset_PDB_coordinates_HA.setMatrix(X1, X2, 0, number_of_conformations - 1, X);
						subset_PDB_coordinates_HA.setMatrix(Y1, Y2, 0, number_of_conformations - 1, Y);
						subset_PDB_coordinates_HA.setMatrix(Z1, Z2, 0, number_of_conformations - 1, Z);

						offset += count;
					}
				return subset_PDB_coordinates_HA;
			}

		public Matrix get_subset_Cartesian_Coords_All_Atom(List<Matrix> res_coords, List<Integer> residues, List<Integer> number_of_Atoms_in_Residues)
			{
				int number_of_conformations = res_coords.get(0).getColumnDimension();
				int number_of_residues_SS = residues.size();
				int number_atoms_SS = 0;
				for (int res : residues)
					{
						System.out.println("Residue: " + res);
						int atom_count = number_of_Atoms_in_Residues.get(res);
						number_atoms_SS += atom_count;
						System.out.println("Residue atoms: " + atom_count);
						System.out.println("Total atoms: " + number_atoms_SS);
					}

				int ROWS_SS = (number_atoms_SS * 3);
				int COLS_SS = number_of_conformations;
				System.out.println("Rows in SS matrix: " + ROWS_SS);

				subset_PDB_coordinates_AA = new Matrix(ROWS_SS, COLS_SS);

				int offset = 0;
				for (int i = 0; i < number_of_residues_SS; i++) // array indices must reference ZERO as first element (Java Array Numbering)
					{
						int residue = residue_list.get(i);
						int count = number_of_Atoms_in_Residues.get(residue);
						// System.out.println("Count = " + count);
						// System.out.println("Offset = " + offset);
						Matrix coords = res_coords.get(residue);
						// System.out.println("Coordinates for Residue :" + residue);
						// System.out.println("Rows/3 :" + coords.getRowDimension() / 3);
						Matrix X = coords.getMatrix(0, count - 1, 0, number_of_conformations - 1);
						Matrix Y = coords.getMatrix(count, (2 * count - 1), 0, number_of_conformations - 1);
						Matrix Z = coords.getMatrix(2 * count, (3 * count - 1), 0, number_of_conformations - 1);
						int X1 = offset;
						int X2 = count + offset - 1;
						int Y1 = number_atoms_SS + offset;
						int Y2 = (count + number_atoms_SS + offset - 1);
						int Z1 = 2 * number_atoms_SS + offset;
						int Z2 = (count + 2 * number_atoms_SS + offset - 1);
						// System.out.println("X1,X2,Y1,Y2,Z1,Z2 " + X1 + "\t" + X2 + "\t" + Y1 + "\t" + Y2 + "\t" + Z1 + "\t" + Z2);
						subset_PDB_coordinates_AA.setMatrix(X1, X2, 0, number_of_conformations - 1, X);
						subset_PDB_coordinates_AA.setMatrix(Y1, Y2, 0, number_of_conformations - 1, Y);
						subset_PDB_coordinates_AA.setMatrix(Z1, Z2, 0, number_of_conformations - 1, Z);

						offset += count;
					}
				return subset_PDB_coordinates_AA;
			}

		public Matrix get_Hierarchical_Subset_Coords(Matrix coords, List<Integer> residues, int number_of_modes_Residues)
			{
				int R1 = 0, R2 = 0, C1 = 0, C2 = 0, offset = 0, row_count = 0;
				int number_of_residues_SS = residues.size();
				int ROWS_SS = (number_of_residues_SS * number_of_modes_Residues);
				int COLS_SS = coords.getColumnDimension();
				int number_of_conformations = COLS_SS;
				subset_Hierarchical_Coordinates = new Matrix(ROWS_SS, COLS_SS);

				for (int i = 0; i < number_of_residues_SS; i++)
					{
						int residue_number = residues.get(i);
						offset = residue_number * number_of_modes_Residues;

						R1 = offset;
						R2 = offset + number_of_modes_Residues - 1;
						C1 = 0;
						C2 = number_of_conformations - 1;

						Matrix res_coords = coords.getMatrix(R1, R2, C1, C2);

						R1 = row_count;
						R2 = row_count + number_of_modes_Residues - 1;
						C1 = 0;
						C2 = number_of_conformations - 1;

						subset_Hierarchical_Coordinates.setMatrix(R1, R2, C1, C2, res_coords);

						row_count += number_of_modes_Residues;
					}
				return subset_Hierarchical_Coordinates;
			}

		// ------------------------------------------------------------------------------------------------------------------------------------

		public Vector<Atom> edit_B_Factors_with_RMSDs(Vector<Atom> atoms, List<Double> rmsds)
			{
				List<Double> sorted_res_rmsds = new ArrayList<>();
				sorted_res_rmsds.addAll(rmsds);
				Collections.sort(sorted_res_rmsds, Collections.reverseOrder());
				double max_rmsd = sorted_res_rmsds.get(0);
				double min_rmsd = sorted_res_rmsds.get(rmsds.size() - 1);
				if (min_rmsd < FLOOR)
					min_rmsd = FLOOR;
				double log_RR_max = Math.log10(max_rmsd);
				double log_RR_min = Math.log10(min_rmsd);
				double delta_x = ((log_RR_max) - (log_RR_min));
				double slope = (delta_y / delta_x);
				double y_min = (slope * log_RR_min);
				int i = 0;
				for (Atom a : atoms)
					{
						a.toString();
						double bff = (rmsds.get(i));
						double log_bff = Math.log10(bff);
						double bf = ((slope * log_bff) - y_min);
						a.b_factor = bf;
						i++;
					}
				return atoms;
			}

		/* ************************************************* GETTERS **************************************************************************** */

		public int getNumber_of_atoms()
			{
				return number_of_atoms;
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

		public Matrix getReference_PDB_coordinates()
			{
				return reference_PDB_coordinates;
			}

		public List<Matrix> getReference_Residue_Coords_List()
			{
				return Reference_Residue_Coords_List;
			}

		public Matrix getReference_subset_PDB_coordinates()
			{
				return reference_subset_PDB_coordinates_AA;
			}

		public Matrix getOriginal_PDB_coordinates()
			{
				return original_PDB_coordinates_AA;
			}

		public Matrix getSubset_PDB_coordinates()
			{
				return subset_PDB_coordinates_AA;
			}

		public Matrix getReference_PDB_coordinates_CA()
			{
				return reference_PDB_coordinates_CA;
			}

		public Matrix getReference_subset_PDB_coordinates_CA()
			{
				return reference_subset_PDB_coordinates_CA;
			}

		public Matrix getOriginal_PDB_coordinates_CA()
			{
				return original_PDB_coordinates_CA;
			}

		public Matrix getSubset_PDB_coordinates_CA()
			{
				return subset_PDB_coordinates_CA;
			}

		public Matrix getReference_PDB_coordinates_HA()
			{
				return reference_PDB_coordinates_HA;
			}

		public Matrix getReference_subset_PDB_coordinates_HA()
			{
				return reference_subset_PDB_coordinates_HA;
			}

		public Matrix getOriginal_PDB_coordinates_HA()
			{
				return original_PDB_coordinates_HA;
			}

		public Matrix getSubset_PDB_coordinates_HA()
			{
				return subset_PDB_coordinates_HA;
			}

		// ------------------------------------------------------------------------- //

		public Vector<Atom> get_ref_Atoms()
			{
				return ref_atoms;
			}

		public Vector<Atom> get_ref_Atoms_CA()
			{
				return ref_atoms_CA;
			}

		public Vector<Atom> get_ref_Atoms_HA()
			{
				return ref_atoms_HA;
			}

		public Vector<Atom> getAtoms()
			{
				return atoms_AA;
			}

		public Vector<Atom> getAtoms_CA()
			{
				return atoms_CA;
			}

		public Vector<Atom> getAtoms_HA()
			{
				return atoms_HA;
			}

		// ------------------------------------------------------------------------- //

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

		public List<Integer> getAtoms_read()
			{
				return atoms_read;
			}

		public List<Integer> getResidues_read()
			{
				return residues_read;
			}

		public List<Integer> getNumber_of_Atoms_in_Residues()
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

		public List<String> getChain_ids()
			{
				return chain_ids;
			}

		public List<String> getChain_ids_read()
			{
				return chain_ids_read;
			}

		// ------------------------------------------------------------------------- //

		public List<Matrix> getResidue_Coords_List()
			{
				return Residue_Coords_List;
			}

		/* ************************************************* SETTERS **************************************************************************** */

		public void setOut_dir(String out_dir)
			{
				this.out_dir = out_dir;
			}
	}
