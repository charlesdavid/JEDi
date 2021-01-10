/*
 * 
 */
package supportStats;

import java.util.ArrayList;

/**
 *
 * @author jenny
 */
public class Results
{

	OutputControl out = new OutputControl();

	ArrayList<Double> x = new ArrayList<Double>();
	ArrayList<Double> PDF = new ArrayList<Double>();
	ArrayList<Double> CDF = new ArrayList<Double>();
	ArrayList<Double> px = new ArrayList<Double>();
	ArrayList<Double> SQR = new ArrayList<Double>();
	ArrayList<Double> L = new ArrayList<Double>();


	public void createSolution(InputParameters input, InputData data, MinimizeScore solution, Score score)
	{
		boolean power = false;

		double max = data.maximumCalc;
		double min = data.minimumCalc;
		double normFactor = 1;

		double dzSize = data.nPoints;
		double dzUniform = (max - min) * 1.0 / dzSize;
		dzSize++;

		double[] dr;
		dr = data.dz;

		double[] dz;
		dz = new double[(int) dzSize];
		for (int i = 0; i < dzSize; i++)
		{
			dz[i] = dzUniform;
		}

		double[] termsT = new double[solution.mode + 1];
		double[] termsP = new double[(int) dzSize];
		;
		double p = 0;
		double termsSum = 0;
		double z = 0;
		double q = min;
		double[] lagrange = solution.getLagrange();
		for (int k = 0; k < dzSize; k++)
		{
			z = (2 * q - max - min) / (max - min);

			termsT[0] = 1.0;
			termsT[1] = z;

			p = lagrange[0];
			if (power)
			{
				p = 1.0 / Math.pow(q, lagrange[1]);
			} else
			{
				for (int t = 1; t < solution.mode; t++)
				{
					p += termsT[t] * lagrange[t];
					termsT[t + 1] = 2 * z * termsT[t] - termsT[t - 1];
				}
				p = Math.exp(p);// + solution->normalize);
				p /= (max - min) / 2;
			}
			termsP[k] = p;
			termsSum += p * dz[k];
			q += dz[k];
		}
		termsSum /= normFactor;

		double lambdaZero = -Math.log(termsSum);// solution->normalize;//

		L.add(lambdaZero);
		for (int j = 1; j < solution.mode; j++)
		{
			L.add(lagrange[j]);
			out.print("Lagrange   ", L.get(j));
		}

		for (int k = 0; k < dzSize; k++)
		{
			double pk = termsP[k] / termsSum;
			termsP[k] = pk;
		}
		double[] pdf = new double[(int) dzSize];
		double dzBig[] = new double[(int) dzSize];
		termsSum = 0;
		int count = 0;

		for (int k = 0; k < dzSize; k++)
		{
			termsSum += termsP[k] * dz[k];
			CDF.add(termsSum);
			dzBig[k] = dz[k];
			pdf[k] = termsP[k];
			count++;
		}

		q = min;

		for (int k = 0; k < count; k++)
		{
			x.add(q);
			PDF.add(pdf[k]);
			px.add(pdf[k] * dr[k]); // added this to get probs for JEDi
			// px.add(pdf[k] * dzUniform); // added this to get probs for JEDi

			if (k < count)
			{
				q += dzBig[k];
			}
		}
	}
}
