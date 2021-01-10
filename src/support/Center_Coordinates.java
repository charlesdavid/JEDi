package support;

import Jama.Matrix;

/**
 * Class Center_Coordinates:
 * 
 * Centers the Cartesian coordinates in a PDB file using the X, Y, and Z centroids.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */
public class Center_Coordinates
	{

		final int COLS = 1;
		int ROWS, number_of_alpha_carbons;
		public double x_centroid;
		public double y_centroid;
		public double z_centroid;
		Matrix coords, centered_coords;

		/**
		 * Centers the coordinates that are in the column vector.
		 * 
		 * @param col_vector A single frame (of alpha carbons) obtained from a PDB file.
		 */
		public Center_Coordinates(Matrix col_vector)
			{

				coords = col_vector;
				ROWS = coords.getRowDimension();
				number_of_alpha_carbons = (ROWS / 3);

			}

		/**
		 * Returns the centered coordinates.
		 * 
		 * @return centered_coords
		 */
		public Matrix get_centered_coordinates()
			{

				double sum_x = 0;
				double sum_y = 0;
				double sum_z = 0;
				for (int i = 0; i < number_of_alpha_carbons; i++)
					{
						double x = coords.get(i, 0);
						sum_x += x;
						double y = coords.get((i + number_of_alpha_carbons), 0);
						sum_y += y;
						double z = coords.get((i + (2 * number_of_alpha_carbons)), 0);
						sum_z += z;
					}
				x_centroid = (sum_x / number_of_alpha_carbons);
				y_centroid = (sum_y / number_of_alpha_carbons);
				z_centroid = (sum_z / number_of_alpha_carbons);
				centered_coords = new Matrix(ROWS, COLS);
				for (int i = 0; i < number_of_alpha_carbons; i++)
					{
						double xrc = (coords.get(i, 0) - x_centroid);
						double yrc = (coords.get((i + number_of_alpha_carbons), 0) - y_centroid);
						double zrc = (coords.get((i + (2 * number_of_alpha_carbons)), 0) - z_centroid);
						centered_coords.set(i, 0, xrc);
						centered_coords.set((i + number_of_alpha_carbons), 0, yrc);
						centered_coords.set((i + (2 * number_of_alpha_carbons)), 0, zrc);
					}
				return centered_coords;
			}

		/**
		 * Returns the X centroid of the frame
		 * 
		 * @return
		 */
		public double getX_centroid()
			{

				return x_centroid;
			}

		/**
		 * Returns the Y centroid of the frame
		 * 
		 * @return
		 */
		public double getY_centroid()
			{

				return y_centroid;
			}

		/**
		 * Returns the Z centroid of the frame
		 * 
		 * @return
		 */
		public double getZ_centroid()
			{

				return z_centroid;
			}
	}
