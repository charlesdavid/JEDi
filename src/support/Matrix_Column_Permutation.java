package support;

import java.util.ArrayList;
import java.util.Collections;

import Jama.Matrix;


public class Matrix_Column_Permutation
{

	final int ROWS, COLS;
	final Matrix data_in, data_permuted;
	final ArrayList<Integer> indices, indices_Permuted;
	final ArrayList<Double> Components, Components_Permuted;

	public Matrix_Column_Permutation(Matrix m)
	{
		data_in = m;
		ROWS = data_in.getRowDimension();
		COLS = data_in.getColumnDimension();
		data_permuted = new Matrix(ROWS, COLS);

		indices = new ArrayList<Integer>();
		indices_Permuted = new ArrayList<Integer>();

		Components = new ArrayList<Double>();
		Components_Permuted = new ArrayList<Double>();
	}

	public Matrix Get_Permuted_Matrix_by_Col() // Performs a systematic column permutation
	{
		for (int j = 0; j < COLS; j++)
		{
			int row_index_1 = 0;
			int row_index_2 = ROWS - 1;
			int col_index_1 = j;
			int col_index_2 = j;

			Matrix col = data_in.getMatrix(row_index_1, row_index_2, col_index_1, col_index_2);
			Matrix col_permuted = new Matrix(ROWS, 1);
			for (int k = 0; k < ROWS; k++)
			{
				double value = col.get(k, 0);
				if ((k + 1) < ROWS)
				{
					col_permuted.set(k + 1, 0, value);
				} else
				{
					col_permuted.set(0, 0, value);
				}
			}
			data_permuted.setMatrix(row_index_1, row_index_2, col_index_1, col_index_2, col_permuted);
		}
		return data_permuted;
	}

	public Matrix Get_Random_Permuted_Matrix_by_Col() // shuffles EACH column randomly
	{

		for (int j = 0; j < COLS; j++)
		{
			int row_index_1 = 0;
			int row_index_2 = ROWS - 1;
			int col_index_1 = j;
			int col_index_2 = j;

			Matrix col = data_in.getMatrix(row_index_1, row_index_2, col_index_1, col_index_2);
			Matrix col_permuted = new Matrix(ROWS, 1);
			for (int i = 0; i < ROWS; i++)
			{
				double element = col.get(i, 0);
				Components.add(element);
				Components_Permuted.add(element);
			}
			Collections.shuffle(Components_Permuted);

			for (int k = 0; k < ROWS; k++)
			{
				double value = Components_Permuted.get(k);
				col_permuted.set(k, 0, value);
			}
			data_permuted.setMatrix(row_index_1, row_index_2, col_index_1, col_index_2, col_permuted);
		}
		return data_permuted;
	}

	public Matrix Get_Random_Permuted_Matrix() // shuffles EVERY column according to the SAME random map
	{
		for (int i = 0; i < ROWS; i++)
		{
			indices.add(i);
			indices_Permuted.add(i);
		}
		Collections.shuffle(indices_Permuted);

		for (int j = 0; j < COLS; j++)
		{
			int row_index_1 = 0;
			int row_index_2 = ROWS - 1;
			int col_index_1 = j;
			int col_index_2 = j;

			Matrix col = data_in.getMatrix(row_index_1, row_index_2, col_index_1, col_index_2);
			Matrix col_permuted = new Matrix(ROWS, 1);

			for (int i = 0; i < ROWS; i++)
			{
				double element = col.get(i, 0);
				int index = indices_Permuted.get(i);
				col_permuted.set(index, 0, element);
			}
			data_permuted.setMatrix(row_index_1, row_index_2, col_index_1, col_index_2, col_permuted);
		}
		return data_permuted;
	}
}
