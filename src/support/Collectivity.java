package support;

import Jama.Matrix;

public class Collectivity
{
	public static final double MACH_EPS = 1.00E-16;

	public static Matrix get_Collectivity(Matrix sq_modes)
	{
		int number_of_modes = sq_modes.getColumnDimension();
		int rows = sq_modes.getRowDimension();

		Matrix ev_collectivity = new Matrix(number_of_modes, 1);

		for (int i = 0; i < number_of_modes; i++)
		{
			int row_index_1 = 0;
			int row_index_2 = rows - 1;
			int col_index_1 = i;
			int col_index_2 = i;

			double collectivity = 0;
			double sum_info = 0;

			Matrix col = sq_modes.getMatrix(row_index_1, row_index_2, col_index_1, col_index_2);
			for (int j = 0; j < rows; j++)
			{
				double element = col.get(j, 0);
				if (Math.abs(element) < MACH_EPS)
				{
					element = MACH_EPS;
				}
				double info = element * (Math.log(element));
				sum_info += info;
			}
			collectivity = ((Math.exp(-sum_info)) / rows);
			ev_collectivity.set(i, 0, collectivity);
		}
		return ev_collectivity;
	}
}
