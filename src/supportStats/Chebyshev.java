/*
 * To change this license header, choose License Headers in Project Properties. To change this template file, choose Tools | Templates and open the template in the editor.
 */
package supportStats;

import java.util.ArrayList;

/**
 *
 * @author jenny
 */
public class Chebyshev
{

	int size;
	double[] dz;

	ArrayList<double[]> termsT = new ArrayList<double[]>();

	public void initialize(double dzLocal[], int sizeLocal)
	{
		size = sizeLocal;
		dz = dzLocal;
		double[] zeroT = new double[size];
		double[] oneT = new double[size];
		double[] twoT = new double[size];

		for (int z = 0; z < size; z++)
			{
				double x = -1 + dz[z] * 2;
				zeroT[z] = 1;
				oneT[z] = x;
				twoT[z] = 2 * x * oneT[z] - 1;
			}
		termsT.add(zeroT);
		termsT.add(oneT);
		termsT.add(twoT);
	}

	public double[] getTerms(int mode)
	{
		try
			{
				return termsT.get(mode);
			}
		catch (IndexOutOfBoundsException e)
			{
				return addMode(mode);
			}
	}

	public ArrayList<double[]> getAllTerms(int mode)
	{
		for (int i = 3; i < mode; i++)
			{
				if (termsT.size() <= i)
					{
						addMode(i);
					}
			}
		return termsT;
	}


	private double[] addMode(int mode)
	{

		double T[] = termsT.get(mode - 1);
		double Tprev[] = termsT.get(mode - 2);
		double Tnext[] = new double[size];
		double x = 0;

		for (int z = 0; z < size; z++)
			{
				x = -1 + 2 * dz[z];
				Tnext[z] = 2 * x * T[z] - Tprev[z];
			}
		termsT.add(Tnext);
		return Tnext;
	}
}
