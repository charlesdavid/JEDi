package support;

import Jama.Matrix;

/**
 * Class KMO_MSA: Uses the CORR(R), and PCORR(P) matrices to calculate the Kaiser-Meyer-Olkin measure and the Measure of Sampling Adequacy.
 * 
 * Copyright (C) 2019 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class KMO_MSA
{
	public static double KMO;

	public static Matrix get_MSA(Matrix CORR, Matrix PCORR)
	{
		int ROWS = CORR.getRowDimension();
		Matrix MSA = new Matrix(ROWS, 1);
		double sum_sum_rr = 0;
		double sum_sum_pp = 0;
		for (int j = 0; j < ROWS; j++)
			{
				double sum_rr = 0;
				double sum_pp = 0;
				for (int k = j + 1; k < ROWS; k++)
					{
						double rUT = CORR.get(j, k);
						double rLT = CORR.get(k, j);
						double rrUT = (rUT * rUT);
						double rrLT = (rLT * rLT);
						sum_rr += (rrUT + rrLT);
						sum_sum_rr += sum_rr;

						double pUT = PCORR.get(j, k);
						double pLT = PCORR.get(k, j);
						double ppUT = (pUT * pUT);
						double ppLT = (pLT * pLT);
						sum_pp += (ppUT + ppLT);
						sum_sum_pp += sum_pp;
					}
				if (j == ROWS - 1)
					{
						for (int k = 0; k < j; k++)
							{
								double r = CORR.get(j, k);
								double rr = (r * r);
								sum_rr += (rr);
								sum_sum_rr += sum_rr;

								double p = PCORR.get(j, k);
								double pp = (p * p);
								sum_pp += (pp);
								sum_sum_pp += sum_pp;
							}

					}
				double MSAj = (sum_rr) / (sum_rr + sum_pp);
				MSA.set(j, 0, MSAj);
			}
		KMO = (sum_sum_rr) / (sum_sum_rr + sum_sum_pp);
		return MSA;
	}

	public static double getKMO()
	{
		return KMO;
	}
}