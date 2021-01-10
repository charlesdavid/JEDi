package jedi;

import java.util.List;

import Jama.Matrix;

/**
 * JED class JED_Get_Distances_for_Residue_Pairs: Constructs the matrix of distances for the Residue Pairs subset. Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Distances_for_Atom_Pairs
{
	final int ROWS, ROWS_dp, COLS, number_of_atoms, number_of_atom_pairs;
	Matrix X_vectors, ref_distances, distances;
	List<Integer> atom_List1, atom_List2;
	boolean exist;

	/* ***************************************** CONSTRUCTOR ********************************************************************* */

	public JEDi_Get_Distances_for_Atom_Pairs(Matrix data, List<Integer> atm_list1, List<Integer> atm_list2)
	{
		this.X_vectors = data;
		this.ROWS = X_vectors.getRowDimension();
		this.COLS = X_vectors.getColumnDimension();
		this.number_of_atoms = (ROWS / 3);
		this.atom_List1 = atm_list1;
		this.atom_List2 = atm_list2;

		this.number_of_atom_pairs = atom_List1.size();
		this.ROWS_dp = atom_List1.size();
	}


	/* ******************************************* METHODS ********************************************************************* */

	public Matrix Get_Ref_Distances()
	{

		ref_distances = new Matrix(ROWS_dp, 1);

		for (int i = 0; i < number_of_atom_pairs; i++)
			{
				int atom1 = atom_List1.get(i);
				int atom2 = atom_List2.get(i);

				double x = X_vectors.get(atom1, 0);
				double y = X_vectors.get(atom1 + number_of_atoms, 0);
				double z = X_vectors.get(atom1 + 2 * number_of_atoms, 0);

				double xr = X_vectors.get(atom2, 0);
				double yr = X_vectors.get(atom2 + number_of_atoms, 0);
				double zr = X_vectors.get(atom2 + 2 * number_of_atoms, 0);

				double xxr = Math.pow((x - xr), 2);
				double yyr = Math.pow((y - yr), 2);
				double zzr = Math.pow((z - zr), 2);

				double sum_xyz = xxr + yyr + zzr;
				double d = Math.sqrt(sum_xyz);
				ref_distances.set(i, 0, d);
			}
		return ref_distances;
	}

	public Matrix get_Distances()
	{

		distances = new Matrix(ROWS_dp, COLS);
		for (int frame = 0; frame < COLS; frame++)
			{
				for (int i = 0; i < number_of_atom_pairs; i++)
					{
						int atom1 = atom_List1.get(i);
						int atom2 = atom_List2.get(i);

						double x = X_vectors.get(atom1, frame);
						double y = X_vectors.get(atom1 + number_of_atoms, frame);
						double z = X_vectors.get(atom1 + 2 * number_of_atoms, frame);

						double xr = X_vectors.get(atom2, frame);
						double yr = X_vectors.get(atom2 + number_of_atoms, frame);
						double zr = X_vectors.get(atom2 + 2 * number_of_atoms, frame);

						double xxr = Math.pow((x - xr), 2);
						double yyr = Math.pow((y - yr), 2);
						double zzr = Math.pow((z - zr), 2);

						double sum_xyz = xxr + yyr + zzr;
						double d = Math.sqrt(sum_xyz);
						distances.set(i, frame, d);
					}
			}
		return distances;
	}
}
