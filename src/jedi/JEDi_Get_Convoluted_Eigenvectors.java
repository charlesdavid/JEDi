package jedi;

import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

/**
 * JED class JEDi_Get_Convoluted_Eigenvectors:
 * 
 * Class for convoluting the Residue eigenvectors (V) with the Hierarchical Eigenvectors (U).
 * 
 * Copyright (C) 2012 Dr. Charles David
 *
 * @author Dr. Charles David
 */

public class JEDi_Get_Convoluted_Eigenvectors
{
	int number_of_atoms;
	final int number_of_residues, number_of_modes_Residue, number_of_modes_Hierarchical;
	Matrix convoluted_Eigenvectors, convoluted_DVs, convoluted_DVPs;
	final Matrix hierarchical_Eigenvectors;
	final List<Matrix> convoluted_residue_Eigenvectors, residue_Eigenvectors;

	NumberFormat df;
	RoundingMode rm;

	// ********************************************* CONSTRUCTOR ******************************************************************* */

	public JEDi_Get_Convoluted_Eigenvectors(Matrix hierarchical_evects, List<Matrix> residue_evects_list, int modes_eigen_residues)
	{
		this.hierarchical_Eigenvectors = hierarchical_evects;
		this.residue_Eigenvectors = residue_evects_list;
		this.number_of_residues = residue_Eigenvectors.size();
		this.number_of_modes_Hierarchical = hierarchical_Eigenvectors.getColumnDimension();
		this.number_of_modes_Residue = modes_eigen_residues; // added due to separation of variables

		this.convoluted_residue_Eigenvectors = new ArrayList<Matrix>();
	}

	// ************************************************ METHODS ******************************************************************* */

	public Matrix get_Convoluted_Eigenvectors() // Convolutes the Eigen-Residues (V) with the corresponding Residue Generalized Coordinate-Eigenvectors (U).
	{
		int U_offset = 0, G_offset = 0;
		for (Matrix V : residue_Eigenvectors) // Iterate through the residue eigenvectors
		{
			int rows_V = V.getRowDimension();
			Matrix G = new Matrix(rows_V, number_of_modes_Hierarchical);

			for (int j = 0; j < number_of_modes_Hierarchical; j++) // Iterating across the U matrix
			{
				Matrix Uk_col = hierarchical_Eigenvectors.getMatrix(U_offset, U_offset + number_of_modes_Residue - 1, j, j); // the part of U_col that corresponds to the residue.
				Matrix g = Uk_col.transpose().times(V.transpose()); // Apply the weights using matrix multiplication
				G.setMatrix(0, rows_V - 1, j, j, g.transpose());
			}
			convoluted_residue_Eigenvectors.add(G);
			U_offset += number_of_modes_Residue;
			G_offset += rows_V;
		}
		convoluted_Eigenvectors = new Matrix(G_offset, number_of_modes_Hierarchical);
		number_of_atoms = (G_offset / 3);

		int big_G_offset = 0;
		for (int i = 0; i < number_of_residues; i++)
		{
			Matrix gK = convoluted_residue_Eigenvectors.get(i);
			int num_of_atoms_Res = gK.getRowDimension() / 3;

			for (int j = 0; j < number_of_modes_Hierarchical; j++)
			{
				for (int k = 0; k < num_of_atoms_Res; k++)
				{
					double X = gK.get(k, j);
					double Y = gK.get(k + num_of_atoms_Res, j);
					double Z = gK.get(k + 2 * num_of_atoms_Res, j);

					convoluted_Eigenvectors.set(k + big_G_offset, j, X);
					convoluted_Eigenvectors.set(k + big_G_offset + number_of_atoms, j, Y);
					convoluted_Eigenvectors.set(k + big_G_offset + 2 * number_of_atoms, j, Z);
				}
			}
			big_G_offset += num_of_atoms_Res;
		}
		return convoluted_Eigenvectors;
	}

	public int getNumber_of_atoms()
	{
		return number_of_atoms;
	}
}