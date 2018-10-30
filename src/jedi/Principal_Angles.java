package jedi;

import java.util.ArrayList;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

/**
 * JED class Principal_Angles: This class computes the principle angles from a matrix of overlaps (projections) using SVD. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

public class Principal_Angles
	{

		double ss_dim, max_angle;
		Matrix input_matrix1, input_matrix2, projections;
		ArrayList<Double> principle_angles, principle_angles_rads, cosine_products;

		/* ***************************** CONSTRUCTORS *************************************** */

		/**
		 * Constructor to compute the principal angles between two equidimensional subspaces of a vector space. This class takes as input the matrix of inner products between the
		 * two orthonormal bases. If the two subspaces are 20 dimensional, the matrix of inner products will be a 20 x 20 square symmetric matrix
		 * 
		 * @param m The matrix of inner products (projections).
		 */
		public Principal_Angles(Matrix m)
			{
				projections = m;
				ss_dim = projections.getColumnDimension();
				principle_angles = new ArrayList<>();
				principle_angles_rads = new ArrayList<>();
				cosine_products = new ArrayList<>();
				get_Principle_Angles();
			}

		/* ***************************** METHODS *************************************** */

		/**
		 * This method computes the Principal Angles (PAs) using the SVD of the matrix of inner products. The assumption is that the matrix of inner products was produced from two
		 * orthonormal sets of eigenvectors. IF the sets of eigenvectors are NOT orthonormal, then the singular values will NOT be in the range [0,1] This method will set any
		 * values out of range of the cosine function to the closest value in range.
		 * 
		 * @param projections The matrix of inner products
		 * 
		 * @return The PAs in degrees
		 */
		private ArrayList<Double> get_Principle_Angles()
			{

				SingularValueDecomposition svd = new SingularValueDecomposition(projections);
				double[] singular_values = svd.getSingularValues();
				int num_Of_Angles = singular_values.length;

				for (int i = 0; i < num_Of_Angles; i++)
					{
						double d = singular_values[i];
						if (d > 1)
							d = 1.000;
						if (d < 0)
							d = 0.000;

						// double sqrt_d = Math.sqrt(d);
						double angle = Math.acos(d) * (180 / Math.PI);

						principle_angles_rads.add(d);
						principle_angles.add(angle);

						double cosine_product = 1;
						for (int index = 0; index <= i; index++)
							{
								cosine_product = cosine_product * principle_angles_rads.get(index);
							}
						cosine_products.add(cosine_product);
					}
				return principle_angles;
			}

		/* ******************************************* GETTERS ************************************************************ */

		/**
		 * @return The maximum hyper angle between subspaces of the specified dimension
		 */
		public double get_Max_Angle()
			{
				max_angle = (Math.sqrt(ss_dim) * 90);
				if (ss_dim == 2)
					max_angle = 90;
				return max_angle;
			}

		/**
		 * @return The list of cosine products; for subspaces of dimension n, there will be n products where each is the product of the first k PAs (in rads).
		 */
		public ArrayList<Double> get_Cosine_Products()
			{
				return cosine_products;
			}

		/**
		 * @return The Principle Angles in degrees
		 */
		public ArrayList<Double> get_Principle_Angles_Degrees()
			{
				return principle_angles;
			}
	}
