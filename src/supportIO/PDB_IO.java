package supportIO;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.List;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import support.Atom;
import support.Atom_ID_Pair;
import support.Residue_ID_Pair;

/**
 * JED class PDB_IO: Handles the reading and writing of PDB files using a Fortran Formatter and a PDB File Parser. Catches Exceptions thrown by the All Atom PDB File Parser class.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class PDB_IO
{
	List<Atom> atoms_AA, atoms_BB, atoms_CA, atoms_HA, residue_atoms;
	List<String> chain_ids_read, Lines;
	List<Integer> atoms_read, atoms_backbone_read, atoms_heavy_read, atoms_alpha_carbon_read, residues_read, number_of_atoms_in_Residues, number_of_heavy_atoms_in_Residues;
	List<Residue_ID_Pair> residue_ID_pairs;
	List<Atom_ID_Pair> atom_ID_pairs;
	BufferedReader pdb_reader;
	BufferedWriter pdb_writer;
	StringBuilder SEQRES, CONECT;
	final String path;
	final File pdb;
	final FortranFormat formatter;
	final PDB_File_Parser parser;

	/* ******************************************* CONSTRUCTORS ******************************************************************* */

	public PDB_IO(String path_to_file)
	{
		this.path = path_to_file;

		pdb = new File(path);
		try
			{
				pdb_reader = new BufferedReader(new FileReader(pdb));
			}
		catch (FileNotFoundException e)
			{
				System.err.println("Could not find the file: " + path_to_file);
				e.printStackTrace();
			}
		formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		formatter.setAddReturn(true);
		parser = new PDB_File_Parser();
	}

	public PDB_IO(List<String> lines)
	{
		this.Lines = lines;
		this.pdb = null;
		this.path = null;

		formatter = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		formatter.setAddReturn(true);
		parser = new PDB_File_Parser();
	}

	/* *********************************************** METHODS ******************************************************************* */

	public List<Atom> parse_PDB_File()
	{
		atoms_AA = parser.parse_PDB(pdb_reader, formatter);
		atoms_BB = parser.get_Backbone_Atoms();
		atoms_HA = parser.get_Heavy_Atoms();
		atoms_CA = parser.get_Alpha_Carbons();

		chain_ids_read = parser.get_chain_IDs_Read();
		residue_ID_pairs = parser.get_Residue_ID_pairs();
		atom_ID_pairs = parser.get_Atom_ID_pairs();

		residues_read = parser.get_Residues_Read(); // list of RESIDUE NUMBERS read from the PDB file
		atoms_read = parser.get_Atoms_Read(); // list of ATOM NUMBERS for all atoms read from the PDB file
		atoms_backbone_read = parser.get_Backbone_Atoms_Read(); // list of ATOM NUMBERS for all BACKBONE atoms read from the PDB file
		atoms_heavy_read = parser.get_Heavy_Atoms_Read(); // list of HEAVY ATOM NUMBERS read from the PDB file
		atoms_alpha_carbon_read = parser.get_alpha_carbons_read(); // list of ALPHA CARBON ATOM NUMBERS read from the PDB file

		return atoms_AA;
	}

	public List<Atom> parse_PDB_Lines_List()
	{
		atoms_AA = parser.parse_PDB(Lines, formatter);
		atoms_BB = parser.get_Backbone_Atoms();
		atoms_HA = parser.get_Heavy_Atoms();
		atoms_CA = parser.get_Alpha_Carbons();

		chain_ids_read = parser.get_chain_IDs_Read();
		residue_ID_pairs = parser.get_Residue_ID_pairs();
		atom_ID_pairs = parser.get_Atom_ID_pairs();

		residues_read = parser.get_Residues_Read(); // list of RESIDUE NUMBERS read from the PDB file
		atoms_read = parser.get_Atoms_Read(); // list of ATOM NUMBERS for all atoms read from the PDB file
		atoms_backbone_read = parser.get_Backbone_Atoms_Read(); // list of ATOM NUMBERS for all BACKBONE atoms read from the PDB file
		atoms_heavy_read = parser.get_Heavy_Atoms_Read(); // list of HEAVY ATOM NUMBERS read from the PDB file
		atoms_alpha_carbon_read = parser.get_alpha_carbons_read(); // list of ALPHA CARBON ATOM NUMBERS read from the PDB file

		return atoms_AA;
	}

	// ************************************************************************************************************************************************** //

	public static void Write_PDB(String path, List<Atom> atoms)
	{

		File pdb = new File(path);
		FortranFormat fm = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		fm.setAddReturn(true);
		PDB_File_Parser ps = new PDB_File_Parser();

		try
			{
				BufferedWriter pdb_writer = new BufferedWriter(new FileWriter(pdb));
				ps.write_PDB(pdb_writer, atoms, fm);
				pdb_writer.close();
			}
		catch (IOException io)
			{
				System.err.println("IOException Thrown. Could not write the file: " + path);
				io.printStackTrace();
			}
	}

	public static void Write_BZ2_PDB(String path, List<Atom> atoms)
	{
		FortranFormat fm = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		fm.setAddReturn(true);
		PDB_File_Parser ps = new PDB_File_Parser();

		try
			{
				FileOutputStream output = new FileOutputStream(path);
				BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new BZip2CompressorOutputStream(output), "UTF-8"));
				ps.write_PDB(writer, atoms, fm);
				writer.close();
				output.close();
			}
		catch (IOException io)
			{
				System.err.println("IOException Thrown. Could not write the file: " + path);
				io.printStackTrace();
			}
	}

	public static void Write_PDB(String path, List<Atom> atoms, StringBuilder sb1, StringBuilder sb2)
	{
		File pdb = new File(path);
		FortranFormat FF = new FortranFormat("(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)");
		FF.setAddReturn(true);
		PDB_File_Parser fp = new PDB_File_Parser();

		try
			{
				BufferedWriter pdb_writer = new BufferedWriter(new FileWriter(pdb));
				fp.write_PDB(pdb_writer, atoms, FF);
				pdb_writer.close();
			}
		catch (IOException io)
			{
				System.err.println("IOException Thrown. Could not write the file: " + path);
				io.printStackTrace();
			}
	}

	/* ************************************************* SETTERS ******************************************************************* */

	/* ************************************************* GETTERS ******************************************************************* */

	public List<Atom> get_Atoms_AA()
	{
		return atoms_AA;
	}

	public List<Atom> get_Backbone_Atoms()
	{
		return atoms_BB;
	}

	public List<Atom> get_Heavy_Atoms()
	{
		return atoms_HA;
	}

	public List<Atom> get_Atoms_CA()
	{
		return atoms_CA;
	}

	// --------------------------------------------------------------------------------------------------------------

	public List<Integer> get_Atom_Numbers_Read()
	{
		return atoms_read;
	}

	public List<Integer> get_Heavy_Atom_Numbers_Read()
	{
		return atoms_heavy_read;
	}

	public List<Integer> get_Backbone_Atom_Numbers_Read()
	{
		return atoms_backbone_read;
	}

	public List<Integer> get_Alpha_Carbon_Numbers_Read()
	{
		return atoms_alpha_carbon_read;
	}

	public List<Integer> get_Residue_Numbers_Read()
	{
		return residues_read;
	}

	// --------------------------------------------------------------------------------------------------------------

	public List<Integer> get_Number_of_Atoms_in_Residues()
	{
		return number_of_atoms_in_Residues;
	}

	public List<Integer> get_Number_of_Heavy_Atoms_in_Residues()
	{
		return number_of_heavy_atoms_in_Residues;
	}

	// --------------------------------------------------------------------------------------------------------------

	public List<String> get_Chain_IDs_Read()
	{
		return chain_ids_read;
	}

	public List<Residue_ID_Pair> getResidue_ID_pairs()
	{
		return residue_ID_pairs;
	}

	public List<Atom_ID_Pair> getAtom_ID_pairs()
	{
		return atom_ID_pairs;
	}

	public StringBuilder getSb1()
	{
		return SEQRES;
	}

	public StringBuilder getSb2()
	{
		return CONECT;
	}
}
