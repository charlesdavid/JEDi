package jedi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

/**
 * JED class PDB_IO: Handles the reading and writing of PDB files using a Fortran Formatter and a PDB File Parser. Catches Exceptions thrown by the All Atom PDB File Parser class.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class PDB_IO
	{

		String read_directory, write_directory, read_name, write_name, path;
		Vector<Atom> atoms, atoms_CA, atoms_HA, residue_atoms, edited_atoms;
		List<String> chain_ids_read, residue_ID_pairs, atom_ID_pairs;
		List<Integer> atoms_read, residues_read, number_of_atoms_in_Residues, number_of_heavy_atoms_in_Residues;
		File pdb, edited_pdb;
		BufferedReader pdb_reader;
		BufferedWriter pdb_writer;
		FortranFormat formatter;
		All_Atom_PDB_File_Parser parser;

		/* ******************************************* CONSTRUCTORS ******************************************************************* */

		public PDB_IO(String dir, String name)
			{

				read_name = name;
				read_directory = dir;
				pdb = new File(read_directory + read_name);
				try
					{
						pdb_reader = new BufferedReader(new FileReader(pdb));
					} catch (FileNotFoundException e)
					{
						System.err.println("Could not find the file: " + read_directory + read_name);
						e.printStackTrace();
					}
				formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
				formatter.setAddReturn(true);
				parser = new All_Atom_PDB_File_Parser();

			}

		public PDB_IO(String path_to_file)
			{

				this.path = path_to_file;
				pdb = new File(path);
				try
					{
						pdb_reader = new BufferedReader(new FileReader(pdb));
					} catch (FileNotFoundException e)
					{
						System.err.println("Could not find the file: " + path_to_file);
						e.printStackTrace();
					}
				formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
				formatter.setAddReturn(true);
				parser = new All_Atom_PDB_File_Parser();

			}

		/* *********************************************** METHODS ******************************************************************* */

		public Vector<Atom> Read_PDB()
			{
				try
					{
						atoms = parser.parse_PDB(pdb_reader, formatter);
						atoms_CA = parser.get_Alpha_Carbons();
						atoms_HA = parser.get_NCO_Atoms();
						chain_ids_read = parser.get_chain_IDs_Read();
						residues_read = parser.get_Residues_Read();
						atoms_read = parser.get_Atoms_Read();
						residue_ID_pairs = parser.getResidue_ID_pairs();
						atom_ID_pairs = parser.getAtom_ID_pairs();
						number_of_atoms_in_Residues = parser.get_Number_of_atoms_in_Residues();
						number_of_heavy_atoms_in_Residues = parser.getNumber_of_atoms_NCO_in_Residues();
					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB_Add_Chain_ID(String chainID)
			{
				try
					{
						atoms = parser.parse_PDB_Add_Missing_Chain(pdb_reader, formatter, chainID);
						atoms_CA = parser.get_Alpha_Carbons();
						atoms_HA = parser.get_NCO_Atoms();
						chain_ids_read = parser.get_chain_IDs_Read();
						residues_read = parser.get_Residues_Read();

					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB(String path_to_PDB)
			{
				path = path_to_PDB;
				pdb = new File(path);
				try
					{
						pdb_reader = new BufferedReader(new FileReader(pdb));
					} catch (FileNotFoundException e)
					{
						System.err.println("Could not find the file: " + path);
						e.printStackTrace();
					}
				formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
				formatter.setAddReturn(true);
				parser = new All_Atom_PDB_File_Parser();
				try
					{
						atoms = parser.parse_PDB(pdb_reader, formatter);
						atoms_CA = parser.get_Alpha_Carbons();
						atoms_HA = parser.get_NCO_Atoms();

					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + path);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + path);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + path);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + path);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB_Residue(int residue_number)
			{

				try
					{
						residue_atoms = parser.parse_PDB_Residue(pdb_reader, formatter, residue_number);

					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return residue_atoms;
			}

		public Vector<Atom> Read_PDB_Atom_Subset(List<Integer> atm_list)
			{

				try
					{
						atoms = parser.parse_PDB_using_Atom_List(pdb_reader, formatter, atm_list);

					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB_Atom_Subset_Multi(List<String> chainIDs, List<Integer> atms)
			{
				try
					{
						atoms = parser.parse_PDB_using_ChainID_and_Atom_Lists(pdb_reader, formatter, chainIDs, atms);
					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB_Residue_Subset(List<Integer> residue_list)
			{

				try
					{
						atoms = parser.parse_PDB_using_Residue_List(pdb_reader, formatter, residue_list);

					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}

		public Vector<Atom> Read_PDB_Residue_Subset_Multi(List<String> chainIDs, List<Integer> res)
			{
				try
					{
						atoms = parser.parse_PDB_using_ChainID_and_Atom_Lists(pdb_reader, formatter, chainIDs, res);
					} catch (FileNotFoundException e)
					{
						System.err.println("'File Not Found' Exception Thrown: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (NumberFormatException e)
					{
						System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
						e.printStackTrace();
						System.exit(0);
					} catch (IOException e)
					{
						System.err.println("'IO' Exception Thrown. Could not read the file: " + read_directory + read_name);
						e.printStackTrace();
						System.exit(0);
					} catch (ArrayIndexOutOfBoundsException e)
					{
						System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Check that all PDB files in the working directory have the EXACT same number of residues (alpha carbons).");
						e.printStackTrace();
						System.exit(0);
					} catch (Exception e)
					{
						System.err.println("Exception thrown for PDB file: " + read_directory + read_name);
						System.err.println("Carefully Check the FORMAT of the PDB file.");
						e.printStackTrace();
						System.exit(0);
					}
				return atoms;
			}


		// ************************************************************************************************************************************************** //
		public void Write_PDB(String w_dir, String w_name, Vector<Atom> e_atoms)
			{

				write_directory = w_dir;
				write_name = w_name;
				edited_atoms = e_atoms;
				edited_pdb = new File(write_directory + write_name);
				try
					{
						pdb_writer = new BufferedWriter(new FileWriter(edited_pdb));
						parser.write_PDB(pdb_writer, edited_atoms, formatter);
						pdb_writer.close();
					} catch (IOException io)
					{
						System.err.println("IOException Thrown. Could not write the file: " + write_directory + write_name);
						io.printStackTrace();
					}
			}

		public static void Write_PDB(String path, Vector<Atom> atoms)
			{

				File pdb = new File(path);
				FortranFormat fm = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
				fm.setAddReturn(true);
				All_Atom_PDB_File_Parser ps = new All_Atom_PDB_File_Parser();

				try
					{
						BufferedWriter pdb_writer = new BufferedWriter(new FileWriter(pdb));
						ps.write_PDB(pdb_writer, atoms, fm);
						pdb_writer.close();
					} catch (IOException io)
					{
						System.err.println("IOException Thrown. Could not write the file: " + path);
						io.printStackTrace();
					}
			}

		/* ************************************************* GETTERS ******************************************************************* */

		public Vector<Atom> get_Atoms()
			{
				return atoms;
			}

		public List<Integer> get_Atoms_read()
			{
				return atoms_read;
			}

		public Vector<Atom> get_Atoms_CA()
			{
				return atoms_CA;
			}

		public Vector<Atom> get_Heavy_Atoms()
			{
				return atoms_HA;
			}

		public List<String> get_Chain_ids_read()
			{
				return chain_ids_read;
			}

		public List<Integer> get_Residues_read()
			{
				return residues_read;
			}

		public List<Integer> get_Number_of_atoms_in_Residues()
			{
				return number_of_atoms_in_Residues;
			}

		public List<Integer> get_Number_of_Heavy_Atoms_in_Residues()
			{
				return number_of_heavy_atoms_in_Residues;
			}

		public List<String> getResidue_ID_pairs()
			{
				return residue_ID_pairs;
			}

		public List<String> getAtom_ID_pairs()
			{
				return atom_ID_pairs;
			}
	}
