package jedi;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

public class JEDi_Get_Big_G
	{
		int number_of_residues, number_modes_hierarchical, number_of_atoms, number_of_modes_SS, number_of_modes_Residues, number_of_conformations;
		String directory, out_dir, description;
		Matrix eigenvectors, big_G, big_DVs, DVPs;
		List<Matrix> residue_Gs, residue_Eigenvectors;
		boolean exist;
		boolean success;

		// ********************************************* CONSTRUCTOR ******************************************************************* */

		public JEDi_Get_Big_G(String dir, String desc, int number_of_modes_SS, int modes_res, Matrix U_evects, List<Matrix> residue_evects)
			{
				this.directory = dir;
				this.description = desc;
				this.number_of_modes_SS = number_of_modes_SS;
				this.number_of_modes_Residues = modes_res;
				this.eigenvectors = U_evects;
				this.residue_Eigenvectors = residue_evects;
				this.number_of_residues = residue_Eigenvectors.size();
				residue_Gs = new ArrayList<Matrix>();

			}

		// ************************************************ METHODS ******************************************************************* */

		public Matrix get_big_G_Row() // weights the eigenResidues by the corresponding components of the rgcEigenvectors.
			{
				int U_offset = 0;
				int G_offset = 0;
				number_modes_hierarchical = eigenvectors.getColumnDimension();
				for (Matrix V : residue_Eigenvectors) // Get the residue eigenvectors
					{
						int rows_V = V.getRowDimension();

						Matrix G = new Matrix(rows_V, number_of_modes_Residues);
						for (int i = 0; i < number_of_modes_Residues; i++) // moving down the U matrix
							{
								Matrix Vk_col = V.getMatrix(0, rows_V - 1, i, i);
								// System.out.println("Vk: " + (i + 1));
								// Vk_col.print(9, 3);
								Matrix Uk_row = eigenvectors.getMatrix(i + U_offset, i + U_offset, 0, number_modes_hierarchical - 1);
								// System.out.println("Uk_row: " + (i + 1));
								// Uk_row.print(9, 6);
								Matrix g = new Matrix(rows_V, 1);
								for (int j = 0; j < number_modes_hierarchical; j++) // moving across the U matrix
									{
										double weight = Uk_row.get(0, j);
										Matrix gk = Vk_col.times(weight);
										g.plusEquals(gk);
										// System.out.println("Weight " + j + " = " + weight);
									}
								// System.out.println("gk: " + (i + 1));
								// g.print(9, 3);
								G.setMatrix(0, rows_V - 1, i, i, g);
							}

						residue_Gs.add(G);
						// G.print(9, 3);
						U_offset += number_of_modes_Residues;
						G_offset += rows_V;
						// System.out.println("U_offset: " + U_offset);
						// System.out.println("G_offset: " + G_offset);
					}
				big_G = new Matrix(G_offset, number_of_modes_Residues);
				number_of_atoms = (G_offset / 3);
				// System.out.println("number_of_atoms: " + number_of_atoms);

				int big_G_offset = 0;
				for (int i = 0; i < number_of_residues; i++)
					{
						Matrix gK = residue_Gs.get(i);
						int num_of_atoms_Res = gK.getRowDimension() / 3;
						// System.out.println("num_of_atoms_Res: " + num_of_atoms_Res);
						for (int j = 0; j < number_of_modes_Residues; j++)
							{
								for (int k = 0; k < num_of_atoms_Res; k++)
									{
										double X = gK.get(k, j);
										double Y = gK.get(k + num_of_atoms_Res, j);
										double Z = gK.get(k + 2 * num_of_atoms_Res, j);
										big_G.set(k + big_G_offset, j, X);
										big_G.set(k + big_G_offset + number_of_atoms - 1, j, Y);
										big_G.set(k + big_G_offset + 2 * number_of_atoms - 1, j, Z);
									}
							}
						big_G_offset += num_of_atoms_Res;
						// System.out.println("Big-G_offset: " + big_G_offset);
					}
				// System.out.println("Big G");
				// big_G.print(12, 6);
				return big_G;
			}

		public Matrix get_big_G_Col() // weights the eigenResidues by the corresponding components of the rgcEigenvectors.
			{
				int U_offset = 0;
				int G_offset = 0;
				int i = 0;
				int j = 0;
				number_modes_hierarchical = eigenvectors.getColumnDimension();
				for (Matrix V : residue_Eigenvectors) // Get the residue eigenvectors
					{
						int rows_V = V.getRowDimension();
						Matrix G = new Matrix(rows_V, number_of_modes_Residues);

						for (j = 0; j < number_modes_hierarchical; j++) // moving across the U matrix
							{
								Matrix Uk_col = eigenvectors.getMatrix(U_offset, U_offset + number_of_modes_Residues - 1, j, j);
								// System.out.println("Uk: " + (j + 1));
								// Uk_col.print(12, 9); //////////////////////////////
								Matrix g = new Matrix(rows_V, 1);
								for (i = 0; i < number_of_modes_Residues; i++) // moving across the V matrix
									{
										Matrix Vk_col = V.getMatrix(0, rows_V - 1, i, i);
										// System.out.println("Vk: " + (i + 1));
										// Vk_col.print(9, 3); ///////////////////////////////////
										double weight = Uk_col.get(i, 0);
										Matrix gk = Vk_col.times(weight); // each residue eigenvector times the weight from the hierarchical eigenvector
										g.plusEquals(gk);
										// System.out.println("Weight " + j + " = " + weight);
										// System.out.println("gk: " + (i + 1));
									}
								// V.print(12, 9);///////////////////////////
								// g.print(12, 9); ///////////////////////////
								G.setMatrix(0, rows_V - 1, j, j, g);
							}

						residue_Gs.add(G);
						// G.print(12, 9);
						U_offset += number_of_modes_Residues;
						G_offset += rows_V;
						// System.out.println("U_offset: " + U_offset);
						// System.out.println("G_offset: " + G_offset);
					}
				big_G = new Matrix(G_offset, number_of_modes_Residues);
				number_of_atoms = (G_offset / 3);
				// System.out.println("number_of_atoms: " + number_of_atoms);

				int big_G_offset = 0;
				for (i = 0; i < number_of_residues; i++)
					{
						Matrix gK = residue_Gs.get(i);
						int num_of_atoms_Res = gK.getRowDimension() / 3;
						// System.out.println("num_of_atoms_Res: " + num_of_atoms_Res);
						for (j = 0; j < number_of_modes_Residues; j++)
							{
								for (int k = 0; k < num_of_atoms_Res; k++)
									{
										double X = gK.get(k, j);
										double Y = gK.get(k + num_of_atoms_Res, j);
										double Z = gK.get(k + 2 * num_of_atoms_Res, j);
										big_G.set(k + big_G_offset, j, X);
										big_G.set(k + big_G_offset + number_of_atoms - 1, j, Y);
										big_G.set(k + big_G_offset + 2 * number_of_atoms - 1, j, Z);
									}
							}
						big_G_offset += num_of_atoms_Res;
						// System.out.println("Big-G_offset: " + big_G_offset);
					}
				// System.out.println("Big G");
				// big_G.print(12, 9);
				return big_G;
			}

		public int getNumber_of_atoms()
			{
				return number_of_atoms;
			}
	}