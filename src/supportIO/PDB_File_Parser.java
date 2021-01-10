package supportIO;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import support.Atom;
import support.Atom_ID_Pair;
import support.Residue_ID_Pair;

/**
 * Class All_Atom_PDB_File_Parser: Parser class to read and write PDB files using Fortran Format.
 * 
 * The exceptions thrown by this class are caught by the PDB_IO Class.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class PDB_File_Parser
{
	String line;
	final String C = "C", CA = "CA", N = "N", O = "O", H = "H", P = "P", S = "S", M = "M";
	final List<Atom> atoms_AA, atoms_BB, atoms_CA, atoms_HA;
	final List<Integer> atoms_read, backbone_atoms_read, heavy_atoms_read, alpha_carbons_read, residues_read;
	final List<String> lines, chain_IDs_read;
	final List<Residue_ID_Pair> residue_ID_pairs;
	final List<Atom_ID_Pair> atom_ID_pairs;

	/* ******************************************* CONSTRUCTOR *************************************************************************** */
	/**
	 * Constructs a PDB File Parser object, initializes all array lists and vectors.
	 */
	public PDB_File_Parser()
	{
		lines = new ArrayList<String>(100000);

		residue_ID_pairs = new ArrayList<Residue_ID_Pair>(1000);
		atom_ID_pairs = new ArrayList<Atom_ID_Pair>(100000);
		chain_IDs_read = new ArrayList<String>(1000);

		atoms_read = new ArrayList<Integer>(100000);
		backbone_atoms_read = new ArrayList<Integer>(10000);
		heavy_atoms_read = new ArrayList<Integer>(50000);
		alpha_carbons_read = new ArrayList<Integer>(1000);

		residues_read = new ArrayList<Integer>(1000);

		atoms_AA = new ArrayList<Atom>(100000);
		atoms_BB = new ArrayList<Atom>(5000);
		atoms_CA = new ArrayList<Atom>(1000);
		atoms_HA = new ArrayList<Atom>(50000);
	}

	/* ********************************************* METHODS ******************************************************************************* */

	public List<Atom> parse_PDB(BufferedReader br, FortranFormat formatter)
	{
		try
			{
				while ((line = br.readLine()) != null) lines.add(line);
				br.close();
			}
		catch (FileNotFoundException e)
			{
				System.err.println("'File Not Found' Exception Thrown. ");
				e.printStackTrace();
				System.exit(0);
			}
		catch (NumberFormatException e)
			{
				System.err.println("'Number Format' Exception Thrown. The format of the PDB file is non-standard.");
				e.printStackTrace();
				System.exit(0);
			}
		catch (IOException e)
			{
				System.err.println("'IO' Exception Thrown. Could not read the file. ");
				e.printStackTrace();
				System.exit(0);
			}
		catch (ArrayIndexOutOfBoundsException e)
			{
				System.err.println("'Array Index Out Of Bounds' Exception thrown for PDB file. ");
				System.err.println("Check that all PDB files have the EXACT same NUMBER of atoms, and in the same ORDER.");
				e.printStackTrace();
				System.exit(0);
			}
		catch (Exception e)
			{
				System.err.println("Exception thrown for PDB file. ");
				System.err.println("Carefully Check the FORMAT of the PDB file.");
				e.printStackTrace();
				System.exit(0);
			}
		for (String line : lines)
			{
				if (line.startsWith("ATOM") || line.startsWith("HETATM"))
					{
						Vector<Object> objects = formatter.parse(line);

						Atom a = new Atom();
						a.header = (String) objects.get(0);
						a.atom_number = (Integer) objects.get(1);
						a.symbol = (String) objects.get(2);
						final String SYMBOL = a.symbol.trim();
						a.res_type = (String) objects.get(4);
						a.chainID = (String) objects.get(5);
						if (a.chainID.equals("") || a.chainID.equals(null)) a.setChainID("A");
						a.res_number = (Integer) objects.get(6);
						a.x = (Double) objects.get(8);
						a.y = (Double) objects.get(9);
						a.z = (Double) objects.get(10);
						if (objects.get(11) != null) a.occupancy = (Double) objects.get(11);
						a.b_factor = (Double) objects.get(12);
						a.code = (String) objects.get(13);

						atoms_AA.add(a);
						atoms_read.add(a.atom_number);
						Atom_ID_Pair atmIDpair = new Atom_ID_Pair(a.chainID, a.atom_number);
						atom_ID_pairs.add(atmIDpair);

						/* ATOM Block:---------------------------------------------------------- */
						if (a.header.equals("ATOM"))
							{
								if (SYMBOL.startsWith(N))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(N))
											{
												backbone_atoms_read.add(a.atom_number);
												Atom bbN = new Atom(a);
												atoms_BB.add(bbN);
											}
									}

								if (SYMBOL.startsWith(C))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(C))
											{
												if (!backbone_atoms_read.contains(a.atom_number)) backbone_atoms_read.add(a.atom_number);
													{
														Atom bbC = new Atom(a);
														atoms_BB.add(bbC);
													}
											}

										if (SYMBOL.equals(CA))
											{
												chain_IDs_read.add(a.chainID);
												residues_read.add(a.res_number);
												Residue_ID_Pair resIDpair = new Residue_ID_Pair(a.chainID, a.res_number);

												if (!residue_ID_pairs.contains(resIDpair)) residue_ID_pairs.add(resIDpair);

												if (!backbone_atoms_read.contains(a.atom_number))
													{
														backbone_atoms_read.add(a.atom_number);
														Atom bbCA = new Atom(a);
														atoms_BB.add(bbCA);
													}
												if (!alpha_carbons_read.contains(a.atom_number) & SYMBOL.contains(CA))
													{
														alpha_carbons_read.add(a.atom_number);
														Atom ac = new Atom(a);
														atoms_CA.add(ac);
													}
											}
									}

								if (SYMBOL.startsWith(O))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(O) || SYMBOL.equals("OXT") || SYMBOL.equals("OC1") || SYMBOL.equals("OC2"))
											{
												backbone_atoms_read.add(a.atom_number);
												Atom bbO = new Atom(a);
												atoms_BB.add(bbO);
											}
									}

								if (SYMBOL.startsWith(S) || SYMBOL.startsWith("P"))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);
									}
							}
						/* HETATM Block: ------------------------------------------------------------------------------------------------------------------------------------ */
						if (a.header.equals("HETATM"))
							{
								if (SYMBOL.startsWith(C) || SYMBOL.startsWith(N) || SYMBOL.startsWith(O) || SYMBOL.startsWith(P) || SYMBOL.startsWith(S) || SYMBOL.startsWith(M))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);
										Residue_ID_Pair resIDpair = new Residue_ID_Pair(a.chainID, a.res_number);
										if (!residue_ID_pairs.contains(resIDpair)) residue_ID_pairs.add(resIDpair);
									}
							}
					}
			}
		return atoms_AA;
	}

	public List<Atom> parse_PDB(List<String> lines, FortranFormat formatter)
	{
		for (String line : lines)
			{
				if (line.startsWith("ATOM") || line.startsWith("HETATM"))
					{
						Vector<Object> objects = formatter.parse(line);

						Atom a = new Atom();
						a.header = (String) objects.get(0);
						a.atom_number = (Integer) objects.get(1);
						a.symbol = (String) objects.get(2);
						final String SYMBOL = a.symbol.trim();
						a.res_type = (String) objects.get(4);
						a.chainID = (String) objects.get(5);
						if (a.chainID.equals("") || a.chainID.equals(null)) a.setChainID("A");
						a.res_number = (Integer) objects.get(6);
						a.x = (Double) objects.get(8);
						a.y = (Double) objects.get(9);
						a.z = (Double) objects.get(10);
						if (objects.get(11) != null) a.occupancy = (Double) objects.get(11);
						a.b_factor = (Double) objects.get(12);
						a.code = (String) objects.get(13);

						atoms_AA.add(a);
						atoms_read.add(a.atom_number);
						Atom_ID_Pair atmIDpair = new Atom_ID_Pair(a.chainID, a.atom_number);
						atom_ID_pairs.add(atmIDpair);

						/* ATOM Block:---------------------------------------------------------- */
						if (a.header.equals("ATOM"))
							{
								if (SYMBOL.startsWith(N))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(N))
											{
												backbone_atoms_read.add(a.atom_number);
												Atom bbN = new Atom(a);
												atoms_BB.add(bbN);
											}
									}

								if (SYMBOL.startsWith(C))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(C))
											{
												if (!backbone_atoms_read.contains(a.atom_number)) backbone_atoms_read.add(a.atom_number);
													{
														Atom bbC = new Atom(a);
														atoms_BB.add(bbC);
													}
											}

										if (SYMBOL.equals(CA))
											{
												chain_IDs_read.add(a.chainID);
												residues_read.add(a.res_number);
												Residue_ID_Pair resIDpair = new Residue_ID_Pair(a.chainID, a.res_number);

												if (!residue_ID_pairs.contains(resIDpair)) residue_ID_pairs.add(resIDpair);

												if (!backbone_atoms_read.contains(a.atom_number))
													{
														backbone_atoms_read.add(a.atom_number);
														Atom bbCA = new Atom(a);
														atoms_BB.add(bbCA);
													}
												if (!alpha_carbons_read.contains(a.atom_number) & SYMBOL.contains(CA))
													{
														alpha_carbons_read.add(a.atom_number);
														Atom ac = new Atom(a);
														atoms_CA.add(ac);
													}
											}
									}

								if (SYMBOL.startsWith(O))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);

										if (SYMBOL.equals(O) || SYMBOL.equals("OXT") || SYMBOL.equals("OC1") || SYMBOL.equals("OC2"))
											{
												backbone_atoms_read.add(a.atom_number);
												Atom bbO = new Atom(a);
												atoms_BB.add(bbO);
											}
									}

								if (SYMBOL.startsWith(S) || SYMBOL.startsWith("P"))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);
									}
							}
						/* HETATM Block: ------------------------------------------------------------------------------------------------------------------------------------ */
						if (a.header.equals("HETATM"))
							{
								if (SYMBOL.startsWith(C) || SYMBOL.startsWith(N) || SYMBOL.startsWith(O) || SYMBOL.startsWith(P) || SYMBOL.startsWith(S) || SYMBOL.startsWith(M))
									{
										heavy_atoms_read.add(a.atom_number);
										Atom ha = new Atom(a);
										atoms_HA.add(ha);
										Residue_ID_Pair resIDpair = new Residue_ID_Pair(a.chainID, a.res_number);
										if (!residue_ID_pairs.contains(resIDpair)) residue_ID_pairs.add(resIDpair);
									}
							}
					}
			}
		return atoms_AA;
	}

	/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------- */

	public void write_PDB(BufferedWriter writer, Vector<Atom> atoms, FortranFormat formatter) throws IOException
	{
		Vector<Object> objects = new Vector<>(15);
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
		writer.write("END");
		writer.close();
	}

	public void write_PDB(BufferedWriter writer, List<Atom> atoms, FortranFormat formatter) throws IOException
	{
		Vector<Object> objects = new Vector<>(15);
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
		writer.write("END");
		writer.close();
	}


	/* *********************************************************** SETTERS ******************************************************************************************** */

	/* *********************************************************** GETTERS ******************************************************************************************** */

	public List<Atom> get_Alpha_Carbons()
	{
		return atoms_CA;
	}

	public List<Atom> get_Heavy_Atoms()
	{
		return atoms_HA;
	}

	public List<Atom> get_Backbone_Atoms()
	{
		return atoms_BB;
	}

	// ----------------------------------------------------

	public List<Integer> get_Atoms_Read()
	{
		return atoms_read;
	}

	public List<Integer> get_Backbone_Atoms_Read()
	{
		return backbone_atoms_read;
	}

	public List<Integer> get_Heavy_Atoms_Read()
	{
		return heavy_atoms_read;
	}

	public List<Integer> get_alpha_carbons_read()
	{
		return alpha_carbons_read;
	}

	public List<Integer> get_Residues_Read()
	{
		return residues_read;
	}

	// ----------------------------------------------------

	public List<String> get_chain_IDs_Read()
	{
		return chain_IDs_read;
	}

	public List<Residue_ID_Pair> get_Residue_ID_pairs()
	{
		return residue_ID_pairs;
	}

	public List<Atom_ID_Pair> get_Atom_ID_pairs()
	{
		return atom_ID_pairs;
	}
}
