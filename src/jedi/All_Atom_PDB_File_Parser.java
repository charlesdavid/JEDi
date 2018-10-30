package jedi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

/**
 * JED class All_Atom_PDB_File_Parser: Parser class to read and write PDB files using Fortran Format. The exceptions thrown by this class are caught by the PDB_IO Class.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

class All_Atom_PDB_File_Parser
	{

		String line;
		final String CA = "CA", N = "N", O = "O", H = "H";
		StringBuilder sb1, sb2;
		int first, last, number_of_atoms;
		List<Integer> atoms_read, heavy_atoms_read, residues_read, atom_list, residue_list, number_of_atoms_in_Residues, number_of_heavy_atoms_in_Residues;
		List<String> lines, chainID_list, chain_IDs_read, residue_ID_pairs, atom_ID_pairs;
		Vector<Atom> atoms, atoms_CA, atoms_HA, residue_atoms;

		/* **************************** CONSTRUCTORS **************************************** */

		/**
		 * Constructs a PDB File Parser object, initializes all array lists.
		 */
		public All_Atom_PDB_File_Parser()
			{
				chain_IDs_read = new ArrayList<String>();
				chainID_list = new ArrayList<String>();
				atom_list = new ArrayList<Integer>();
				residue_list = new ArrayList<Integer>();
				number_of_atoms_in_Residues = new ArrayList<Integer>();
				number_of_heavy_atoms_in_Residues = new ArrayList<Integer>();
				heavy_atoms_read = new ArrayList<Integer>();
				atoms_read = new ArrayList<Integer>();
				residues_read = new ArrayList<Integer>();
				residue_ID_pairs = new ArrayList<String>();
				atom_ID_pairs = new ArrayList<String>();
				lines = new ArrayList<String>();
			}

		/* ******************************* METHODS ******************************************* */

		public Vector<Atom> parse_PDB(BufferedReader br, FortranFormat formatter) throws IOException
			{
				atoms = new Vector<Atom>();
				atoms_CA = new Vector<Atom>();
				atoms_HA = new Vector<Atom>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();
				int count_H = 0;

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{

								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);
								atoms.add(a);
								atoms_read.add(a.atom_number);
								atom_ID_pairs.add(a.chainID + "\t" + a.atom_number);

								if (!(a.symbol.contains((H))))
									{
										atoms_HA.add(a);
										heavy_atoms_read.add(a.atom_number);
									}

								if (a.symbol.contains(H))
									count_H++;

								if (a.symbol.equals(CA))
									{
										chain_IDs_read.add(a.chainID);
										residues_read.add(a.res_number);
										residue_ID_pairs.add(a.chainID + "\t" + a.res_number);
										atoms_CA.add(a);
									}

								if (a.symbol.equals(N))
									{
										first = a.atom_number;
									}

								if (a.symbol.equals(O))
									{
										last = a.atom_number;
										number_of_atoms = (last - first) + 1;
										// System.out.println("The number of atoms in residue " + a.res_number + " is " + number_of_atoms);
										number_of_atoms_in_Residues.add(number_of_atoms);
										number_of_heavy_atoms_in_Residues.add(number_of_atoms - count_H);
										// System.out.println("The number of NCO atoms in residue " + a.res_number + " is " + (number_of_atoms - count_H));
										count_H = 0;
									}

								if (a.symbol.equals("OXT"))
									{
										int last = number_of_atoms_in_Residues.size() - 1;
										int val = number_of_atoms_in_Residues.get(last);
										int val2 = number_of_heavy_atoms_in_Residues.get(last);
										number_of_atoms_in_Residues.set(last, val + 1);
										number_of_heavy_atoms_in_Residues.set(last, val2);
										count_H = 0;
									}
							}

						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		public Vector<Atom> parse_PDB_Add_Missing_Chain(BufferedReader br, FortranFormat formatter, String ID) throws IOException
			{
				atoms = new Vector<>();
				atoms_CA = new Vector<>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								if (a.chainID.isEmpty())
									{
										a.chainID = ID;
										// System.out.println("Chain ID = " + a.chainID);
									}
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);

								if (a.symbol.equals(CA))
									{
										chain_IDs_read.add(a.chainID);
										residues_read.add(a.res_number);
										atoms_CA.add(a);
									}
								atoms.add(a);
								atoms_read.add(a.atom_number);
							}

						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		public Vector<Atom> parse_PDB_Residue(BufferedReader br, FortranFormat formatter, int residue_number) throws IOException
			{
				residue_atoms = new Vector<>();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								int test = (Integer) objects.get(6);
								if (test == residue_number)
									{
										Atom a = new Atom();
										a.header = (String) objects.get(0);
										a.atom_number = (Integer) objects.get(1);
										a.symbol = (String) objects.get(2);
										a.res_type = (String) objects.get(4);
										a.chainID = (String) objects.get(5);
										a.res_number = (Integer) objects.get(6);
										a.x = (Double) objects.get(8);
										a.y = (Double) objects.get(9);
										a.z = (Double) objects.get(10);
										if (objects.get(11) != null)
											a.occupancy = (Double) objects.get(11);
										a.b_factor = (Double) objects.get(12);
										a.code = (String) objects.get(13);
										residue_atoms.add(a);
									}
								if (test > residue_number)
									{
										break;
									}
							}
					}
				return residue_atoms;
			}

		public Vector<Atom> parse_PDB_using_Atom_List(BufferedReader br, FortranFormat formatter, List<Integer> atm_list) throws IOException
			{
				atoms = new Vector<Atom>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();

				atom_list = atm_list;
				int num_of_atoms = atom_list.size();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);
								int key = a.atom_number;
								for (int i = 0; i < num_of_atoms; i++)
									{
										if (atom_list.get(i).equals(key))
											{
												atoms.add(a);
												atoms_read.add(a.atom_number);
											}
									}
							}
						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		public Vector<Atom> parse_PDB_using_ChainID_and_Atom_Lists(BufferedReader br, FortranFormat formatter, List<String> chain_ids, List<Integer> res) throws IOException
			{
				atoms = new Vector<Atom>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();

				chainID_list = chain_ids;
				atom_list = res;
				int num_of_atoms = atom_list.size();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);

								String ID = a.chainID;
								int key = a.res_number;

								for (int i = 0; i < num_of_atoms; i++)
									{
										if (chainID_list.get(i).equals(ID) && residue_list.get(i).equals(key))
											{
												atoms.add(a);
												atoms_read.add(a.atom_number);
											}
									}
							}

						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		public Vector<Atom> parse_PDB_using_Residue_List(BufferedReader br, FortranFormat formatter, List<Integer> residues) throws IOException
			{
				atoms = new Vector<>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();

				residue_list = residues;
				int num_of_residues = residue_list.size();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);
								int key = a.res_number;
								for (int res_indx = 0; res_indx < num_of_residues; res_indx++)
									{
										if (residue_list.get(res_indx).equals(key))
											{
												atoms.add(a);
												atoms_read.add(a.atom_number);
											}
									}
							}
						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		public Vector<Atom> parse_PDB_using_ChainID_and_Residue_Lists(BufferedReader br, FortranFormat formatter, List<String> chain_ids, List<Integer> res) throws IOException
			{
				atoms = new Vector<>();
				sb1 = new StringBuilder();
				sb2 = new StringBuilder();

				chainID_list = chain_ids;
				residue_list = res;
				int num_of_residues = residue_list.size();

				while ((line = br.readLine()) != null)
					lines.add(line);
				br.close();

				for (String file_line : lines)
					{
						if (file_line.startsWith("SEQRES"))
							sb1.append(file_line + "\n");

						if (file_line.startsWith("ATOM"))
							{
								Vector<Object> objects = formatter.parse(file_line);
								Atom a = new Atom();
								a.header = (String) objects.get(0);
								a.atom_number = (Integer) objects.get(1);
								a.symbol = (String) objects.get(2);
								a.res_type = (String) objects.get(4);
								a.chainID = (String) objects.get(5);
								a.res_number = (Integer) objects.get(6);
								a.x = (Double) objects.get(8);
								a.y = (Double) objects.get(9);
								a.z = (Double) objects.get(10);
								if (objects.get(11) != null)
									a.occupancy = (Double) objects.get(11);
								a.b_factor = (Double) objects.get(12);
								a.code = (String) objects.get(13);

								String ID = a.chainID;
								int key = a.res_number;

								for (int res_indx = 0; res_indx < num_of_residues; res_indx++)
									{
										if (chainID_list.get(res_indx).equals(ID) && residue_list.get(res_indx).equals(key))
											{
												atoms.add(a);
												atoms_read.add(a.atom_number);
											}
									}
							}

						if (file_line.startsWith("TER"))
							sb2.append(file_line + "\n");

						if (file_line.startsWith("CONECT"))
							sb2.append(file_line + "\n");
					}
				return atoms;
			}

		/* ----------------------------------------------------------------------------------------------------------------------------------------------------- */

		public void write_PDB(BufferedWriter writer, Vector<Atom> atoms, FortranFormat formatter) throws IOException
			{
				Vector<Object> objects = new Vector<>(15);
				if (sb1 != null)
					writer.write(sb1.toString());
				for (Atom a : atoms)
					{
						objects.clear();
						objects.add(a.header + (a.header.length() == 1 ? "   " : "  "));
						objects.add(a.atom_number);
						objects.add(a.symbol + (a.symbol.length() == 1 ? "   " : "  "));
						objects.add("");
						objects.add(a.res_type);
						objects.add(a.chainID);
						objects.add(a.res_number);
						objects.add(null);
						objects.add(a.x);
						objects.add(a.y);
						objects.add(a.z);
						objects.add(a.occupancy);
						objects.add(a.b_factor);
						objects.add(a.code);
						writer.append(formatter.format(objects));
					}
				if (sb2 != null)
					writer.write(sb2.toString());
				writer.write("END");
				writer.close();
			}

		/* ****************************** GETTERS ************************************************************** */

		public List<String> get_chain_IDs_Read()
			{

				return chain_IDs_read;
			}

		public List<Integer> get_Atoms_Read()
			{

				return atoms_read;
			}

		public List<Integer> get_Atoms_NCO_Read()
			{

				return heavy_atoms_read;
			}

		public List<Integer> get_Residues_Read()
			{

				return residues_read;
			}

		public Vector<Atom> get_Alpha_Carbons()
			{
				return atoms_CA;
			}

		public Vector<Atom> get_NCO_Atoms()
			{
				return atoms_HA;
			}

		public List<Integer> get_Number_of_atoms_in_Residues()
			{
				return number_of_atoms_in_Residues;
			}

		public List<Integer> getNumber_of_atoms_NCO_in_Residues()
			{
				return number_of_heavy_atoms_in_Residues;
			}

		public Vector<Atom> getAtoms()
			{
				return atoms;
			}

		public Vector<Atom> get_Residue_atoms()
			{
				return residue_atoms;
			}

		public int get_Number_of_atoms()
			{
				return number_of_atoms;
			}

		public List<Integer> get_Atom_list()
			{
				return atom_list;
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
