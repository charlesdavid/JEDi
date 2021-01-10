package support;

import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import supportIO.Input_Parameters;
import supportIO.List_IO;

/**
 * Class for selecting atoms based on statistical thresholds: variance, skew, or kurtosis.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Select_Variables_by_Statistical_Threshold
{
	final int COLS, ROWS, numberAtoms;
	final double variance_threshold, skew_threshold, kurtosis_threshold;
	final Matrix stats;
	final List<Atom> atoms;
	final String type;
	List<Integer> atomNumbersVariance, atomNumbersSkew, atomNumbersKurtosis;
	List<Integer> residueNumbersVariance, residueNumbersSkew, residueNumbersKurtosis;
	List<Atom> atoms_VAR, atoms_SKEW, atoms_KURTOSIS;
	NumberFormat nf3;
	RoundingMode rm;

	/* *************************** CONSTRUCTOR *************************************************************************************** */

	/**
	 * Constructor takes a List of Atoms and a 1x4 matrix of variable (coordinate) stats: Mean,Var,Skew,Kurtosis
	 * 
	 * @param original_atoms The List of atoms
	 * @param stat_data      The matrix of stats
	 */
	public Select_Variables_by_Statistical_Threshold(String typeOfAnalysis, List<Atom> original_atoms, Matrix stat_data)
	{
		this.stats = stat_data;
		this.atoms = original_atoms;
		this.type = typeOfAnalysis;

		this.COLS = stats.getColumnDimension();
		this.ROWS = stats.getRowDimension();
		this.numberAtoms = atoms.size();

		this.variance_threshold = Input_Parameters.VARIANCE_THRESHOLD;
		this.skew_threshold = Input_Parameters.SKEW_THRESHOLD;
		this.kurtosis_threshold = Input_Parameters.KURTOSIS_THRESHOLD;

		atoms_VAR = new ArrayList<Atom>();
		atoms_SKEW = new ArrayList<Atom>();
		atoms_KURTOSIS = new ArrayList<Atom>();

		nf3 = NumberFormat.getInstance();
		rm = RoundingMode.HALF_UP;
		nf3.setMaximumFractionDigits(3);
		nf3.setMinimumFractionDigits(3);
		nf3.setRoundingMode(rm);
	}

	/* *************************** METHODS *************************************************************************************** */

	public void doThresholding()
	{
		if (variance_threshold > 0) getAtomsVariance();
		if (skew_threshold > 0) getAtomsSkew();
		if (kurtosis_threshold > 0) getAtomsKurtosis();
	}

	private List<Atom> getAtomsVariance()
	{
		int i = 0;
		for (Atom atm : atoms)
			{
				double var1 = stats.get(i, 1);
				double var2 = stats.get(i + numberAtoms, 1);
				double var3 = stats.get(i + 2 * numberAtoms, 1);

				if (var1 > variance_threshold | var2 > variance_threshold | var3 > variance_threshold)
					{
						Atom a = new Atom(atm);
						atoms_VAR.add(a);
					}
				i++;
				// System.out.println("Mean: " + mean + " Variance: " + var + " Skew: " + skew + " Kurtosis: " + kurtosis);
			}
		atomNumbersVariance = getAtomNumbers(atoms_VAR);
		residueNumbersVariance = getResidueNumbers(atoms_VAR);

		String name = type + "_High_Variance_Atoms_" + nf3.format(variance_threshold) + ".txt";
		String path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(atomNumbersVariance, path);
		name = type + "_High_Variance_Residues_" + nf3.format(variance_threshold) + ".txt";
		path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(residueNumbersVariance, path);

		return atoms_VAR;
	}

	private List<Atom> getAtomsSkew()
	{
		int i = 0;
		for (Atom atm : atoms)
			{
				double skew1 = Math.abs(stats.get(i, 2));
				double skew2 = Math.abs(stats.get(i + numberAtoms, 2));
				double skew3 = Math.abs(stats.get(i + 2 * numberAtoms, 2));

				if (skew1 > skew_threshold | skew2 > skew_threshold | skew3 > skew_threshold)
					{
						Atom a = new Atom(atm);
						atoms_SKEW.add(a);
					}
				i++;
				// System.out.println("Mean: " + mean + " Variance: " + var + " Skew: " + skew + " Kurtosis: " + kurtosis);
			}
		atomNumbersSkew = getAtomNumbers(atoms_SKEW);
		residueNumbersSkew = getResidueNumbers(atoms_SKEW);

		String name = type + "_High_Skew_Atoms_" + nf3.format(skew_threshold) + ".txt";
		String path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(atomNumbersSkew, path);
		name = type + "_High_Skew_Residues_" + nf3.format(skew_threshold) + ".txt";
		path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(residueNumbersSkew, path);

		return atoms_SKEW;
	}

	private List<Atom> getAtomsKurtosis()
	{
		int i = 0;
		for (Atom atm : atoms)
			{
				double kurtosis1 = stats.get(i, 3);
				double kurtosis2 = stats.get(i + numberAtoms, 3);
				double kurtosis3 = stats.get(i + 2 * numberAtoms, 3);

				if (kurtosis1 > kurtosis_threshold | kurtosis2 > kurtosis_threshold | kurtosis3 > kurtosis_threshold)
					{
						Atom a = new Atom(atm);
						atoms_KURTOSIS.add(a);
					}
				i++;
				// System.out.println("Mean: " + mean + " Variance: " + var + " Skew: " + skew + " Kurtosis: " + kurtosis);
			}
		atomNumbersKurtosis = getAtomNumbers(atoms_KURTOSIS);
		residueNumbersKurtosis = getResidueNumbers(atoms_KURTOSIS);

		String name = type + "_High_Kurtosis_Atoms_" + nf3.format(kurtosis_threshold) + ".txt";
		String path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(atomNumbersKurtosis, path);
		name = type + "_High_Kurtosis_Residues_" + nf3.format(kurtosis_threshold) + ".txt";
		path = Input_Parameters.OUT_DIR + name;
		List_IO.write_Integer_List(residueNumbersKurtosis, path);

		return atoms_KURTOSIS;
	}

	private List<Integer> getAtomNumbers(List<Atom> atms)
	{
		List<Integer> atmNums = new ArrayList<Integer>();
		for (Atom atm : atms)
			{
				atmNums.add(atm.getAtom_number());
			}
		return atmNums;
	}

	private List<Integer> getResidueNumbers(List<Atom> atms)
	{
		List<Integer> resNums = new ArrayList<Integer>();
		for (Atom atm : atms)
			{
				int residue_number = atm.getRes_number();
				if (!resNums.contains(residue_number)) resNums.add(residue_number);
			}
		return resNums;
	}
}
