/*
 *
 */
package supportStats;

import java.util.ArrayList;

/**
 *
 * @author jenny
 */
@SuppressWarnings("unused")
public class MinimizeScore
{

	OutputControl out = new OutputControl();

	int mode;
	float duration;
	double bestScore;
	double normalize;
	double[] trialRandom;
	int N;

	private int nPoints;
	private int maxLagrange;
	private int seed;
	private boolean useLast;
	private double y2;
	private Chebyshev cheby;
	private double[] inverse;
	private double[] z;
	private double[] dz;
	private double[] doubleInverse;
	private double[] xUntransform;
	private double[] bestRandom;
	private double[] bestLagrange;
	private double[] rawDataPartition;
	public ArrayList<double[]> T;


	public MinimizeScore()
	{

		normalize = 0;
		useLast = false;
		y2 = 0;
		seed = 12345678;
		nPoints = 0;
		N = 0;
		maxLagrange = 0;
		mode = 0;
		duration = 0;
		bestScore = 0;
	}

	public double[] getLagrange()
	{
		return bestLagrange;
	}

	public boolean minimize(InputParameters input, InputData data, Score score)
	{

		int minLagrange = input.minLagrange;
		maxLagrange = input.maxLagrange;
		int nLagrangeAdd = input.nLagrangeAdd;
		double fractionLagrangeAdd = input.fractionLagrangeAdd;
		int loopMax = input.loopMax;
		double initSigma = input.initSigma;
		double finalSigma = input.finalSigma;
		double decayFactor = input.decayFactor;
		int partitionSize = input.initPartitionSize;

		boolean phaseOne = true;

		int targetPartition = data.N;
		boolean incPartition = true;

		if (phaseOne)
			{
				targetPartition = partitionSize;
				incPartition = false;
			}


		this.inverse = data.inverse;
		this.dz = data.dz;
		this.doubleInverse = data.doubleInverse;
		this.xUntransform = data.xUntransform;
		double[] transformedZeroOne = data.transformedZeroOne;
		this.cheby = data.cheby;
		this.nPoints = data.nPointsAdjust;
		this.N = data.N;
		if (partitionSize > N)
			{
				partitionSize = N;
				targetPartition = N;
			}
		if (partitionSize == 0)
			{
				partitionSize = N;
				targetPartition = N;
			}

		int loopCount = 0;
		mode = minLagrange;
		int inc;
		double sigmaFactor = 1;
		if (phaseOne)
			{
				sigmaFactor = Math.sqrt(partitionSize * 1.0 / targetPartition) * Math.sqrt(1.0 / mode);
			}
		double originalFinalSigma = finalSigma;
		double originalInitSigma = initSigma;
		double currentSigma = initSigma * sigmaFactor;
		finalSigma *= sigmaFactor;
		initSigma *= sigmaFactor;


		double[] trialLagrange = new double[maxLagrange];
		bestLagrange = new double[maxLagrange];

		for (int i = 0; i < maxLagrange; i++)
			{
				trialLagrange[i] = 0;
				bestLagrange[i] = 0;
			}

		T = cheby.getAllTerms(maxLagrange);

		double lastLagrangeScore = 0;
		int lagrangeAddCount = 0;

		trialRandom = new double[N];
		bestRandom = new double[N];
		rawDataPartition = new double[N];

		for (int c = 0; c < N; c++)
			{
				trialRandom[c] = 0;
				bestRandom[c] = 0;
				rawDataPartition[c] = 0;
			}

		double trialScore = 0;
		bestScore = -Double.MAX_VALUE;
		double targetScore = score.targetScore;
		double minimumScore = score.minimumScore;
		double maximumScore = score.maximumScore;

		boolean funnelFinished = false;
		boolean solutionNotFound = false;


		ArrayList<Integer> indices = score.setIndices(N, N);
		for (int i = 0; i < N; i++)
			{
				rawDataPartition[i] = transformedZeroOne[i];
			}
		bestScore = score.calculateScore(rawDataPartition, N, N);
		out.print("initial score SURD", score.SURD);
		out.print("initial variance SURD", score.QZVariance);

		if ((score.getLikelihood() > targetScore) && (score.getLikelihood() < maximumScore))
			{
				out.print("*Uniform Solution found");
				return false;
			}

		indices = score.getIndices(N, partitionSize);
		targetPartition = partitionSize;
		for (int i = 0; i < partitionSize; i++)
			{
				rawDataPartition[i] = transformedZeroOne[indices.get(i)];
			}

		indices = score.setIndices(targetPartition, partitionSize);
		double[] cdf;
		cdf = new double[nPoints];
		calculatePDF(cdf, trialLagrange, mode);
		map(trialRandom, cdf, rawDataPartition, partitionSize);
		boolean continueLooking = true;


		bestScore = score.calculateScore(trialRandom, targetPartition, partitionSize);
		out.print("initial score", score.SURD);
		out.print("initial variance", score.QZVariance);
		out.print("partition size", partitionSize);
		out.print("actual partition size", partitionSize);
		out.print("target size", targetPartition);

		if (score.getLikelihood() > targetScore)
			{
				if (mode == 1)
					{
						if (maxLagrange > 1)
							{
								trialLagrange[0] = 0;
								bestLagrange[0] = 0;
								trialLagrange[1] = 0;
								bestLagrange[1] = 0;
								mode = 2;
							}
						else
							{
								out.error("*Maximum number of lagrange multipliers exceeded", maxLagrange);
							}
					}
				else
					{

					}
			}

		bestRandom[0] = trialRandom[0];
		while (partitionSize <= N)
			{
				while (continueLooking)
					{
						loopCount++;
						cdf = new double[nPoints];
						calculatePDF(cdf, trialLagrange, mode);
						map(trialRandom, cdf, rawDataPartition, partitionSize);

						trialScore = score.calculateScore(trialRandom, targetPartition, partitionSize);
						if (trialScore > bestScore)
							{

								String strOut = "SURD score: " + score.SURD + ";  qz var: " + score.QZVariance + ";  partition size: " + partitionSize + "; target:  "
										+ targetPartition;
								out.print(strOut);
								if (score.getLikelihood() < maximumScore)
									{
										bestScore = trialScore;
									}
								bestLagrange[0] = normalize;
								for (int k = 1; k < mode; k++)
									{
										bestLagrange[k] = trialLagrange[k];
									}
								bestRandom[0] = trialRandom[0];
								if ((score.getLikelihood() > targetScore) && (score.getLikelihood() < maximumScore))
									{
										out.print("*Solution found");
										break;
									}
							}

						if (loopCount > loopMax)
							{
								currentSigma /= decayFactor;
								if (currentSigma < finalSigma) funnelFinished = true;
								loopCount = 0;
							}


						if (funnelFinished)
							{
								if (mode < 5) inc = 1;
								else
									inc = 2;
								mode += inc;
								if (mode > maxLagrange)
									{
										if (score.getLikelihood() > minimumScore)
											{
												out.print("*Lower threshold accepted", score.getLikelihood());
											}
										else
											{
												out.error("PDF Estimator: Maximum number of lagrange multipliers exceeded", maxLagrange);
												solutionNotFound = true;
												continueLooking = false;
											}
										mode -= inc;
										break;
									}
								out.print("*Adding lagrange", mode);
								funnelFinished = false;

								if (Math.abs(bestScore - lastLagrangeScore) / Math.abs(lastLagrangeScore) < fractionLagrangeAdd)
									{
										if (lagrangeAddCount > nLagrangeAdd)
											{
												if ((score.getLikelihood() > minimumScore) && (score.getLikelihood() < maximumScore))
													{
														out.print("*Improvement not found in required number of attempts");
														out.print("*Lower threshold accepted", score.getLikelihood());
														targetScore = score.getLikelihood();
													}
												else
													{
														out.print("*Improvement not found in required number of attempts");
														solutionNotFound = true;
													}
												break;
											}
										else
											{
												lagrangeAddCount++;
											}
									}
								else
									{
										lastLagrangeScore = bestScore;
										lagrangeAddCount = 0;
									}
								currentSigma = initSigma;
								funnelDiffusion(bestLagrange, trialLagrange, mode, currentSigma, 1);
							}
						else
							{
								funnelDiffusion(bestLagrange, trialLagrange, mode, currentSigma, 1);
							}
					}

				if (!continueLooking)
					{
						break;
					}
				if (incPartition)
					{
						if (partitionSize == N) break;
						partitionSize = (partitionSize - 1) * 2 + 1;
						if (partitionSize > N)
							{
								partitionSize = N;
							}
						indices = score.setIndices(N, partitionSize);
						for (int i = 0; i < partitionSize; i++)
							{
								rawDataPartition[i] = transformedZeroOne[indices.get(i)];
							}
						bestScore = -Double.MAX_VALUE;
					}
				else
					{
						if (partitionSize == N) break;
						targetPartition = (targetPartition - 1) * 2 + 1;

						if (targetPartition >= N)
							{
								targetPartition = N;
								incPartition = true;
							}
						indices = score.setIndices(targetPartition, partitionSize);
						bestScore = -Double.MAX_VALUE;
					}

				if (phaseOne)
					{
						double tempSigmaFactor = Math.sqrt(partitionSize * 1.0 / targetPartition) * Math.sqrt(1.0 / mode);
						initSigma = originalInitSigma * tempSigmaFactor;
						finalSigma = originalFinalSigma * tempSigmaFactor;
						currentSigma = initSigma;
					}
			}

		if (solutionNotFound)
			{
				// out.error("Solution not found");
			}

		return solutionNotFound;

	}


	public void calculatePDF(double cdf[], double lagrange[], int modes)
	{
		int pdfPoints = nPoints * 2 - 1;
		double[] pdf;
		pdf = new double[pdfPoints];
		double[] x;
		x = new double[pdfPoints];
		for (int i = 0; i < pdfPoints; i++)
			{
				x[i] = 0;
			}

		for (int k = 0; k < pdfPoints; k++)
			{
				for (int n = 0; n < modes; n++)
					{
						x[k] += lagrange[n] * T.get(n)[k];
					}
				pdf[k] = Math.exp(x[k]);
			}

		int count = 1;
		cdf[0] = (pdf[0]) * dz[0] / 2;

		normalize = 0;
		for (int k = 1; k < nPoints; k++)
			{
				cdf[k] = cdf[k - 1] + (pdf[count - 1] + 4 * pdf[count] + pdf[count + 1]) * dz[k - 1] / 6.0;
				normalize += pdf[count] * dz[k];
				count += 2;
			}

		normalize = -Math.log(normalize);

		double constant = cdf[nPoints - 1];
		for (int k = 0; k < nPoints; k++)
			{
				cdf[k] /= constant;
			}
	}



	public void map(double r[], double cdf[], double rawDataPartition[], int partitionSize)
	{
		double z1;
		double z2;
		double zCalc1;
		double endpointCDF;
		int startPoint = 0;

		for (int k = 0; k < partitionSize; k++)
			{
				int j = startPoint;
				for (j = startPoint; j < nPoints; j++)
					{
						if (rawDataPartition[k] < inverse[j]) break;
					}
				startPoint = j;
				if (j == 0)
					{
						z1 = 0;
					}
				else
					{
						z1 = inverse[j - 1];
					}

				if (j >= nPoints)
					{
						z2 = 1.0;
						endpointCDF = 1.0;
					}
				else
					{
						z2 = inverse[j];
						endpointCDF = cdf[j];
					}
				zCalc1 = (rawDataPartition[k] - z1) / (z2 - z1);
				if (j == 0)
					{
						r[k] = zCalc1 * (endpointCDF);
					}
				else
					{
						r[k] = cdf[j - 1] + zCalc1 * (endpointCDF - cdf[j - 1]);
					}
				if (r[k] < 0)
					{
						out.error("ERROR: random number is negative\n");
					}
			}
	}


	public void funnelDiffusion(double original[], double updated[], int arraySize, double currentSigmaMu, int startIndex)
	{
		for (int j = startIndex; j < arraySize; j++)
			{
				updated[j] = random(original[j], (currentSigmaMu * (0.1 * Math.abs(original[j]) + 1.0) / 2));
			}
	}


	public void funnelDiffusion(double original[], double updated[], int arraySize, double currentSigmaMu)
	{
		funnelDiffusion(original, updated, arraySize, currentSigmaMu, 1);
	}

	public double random(double m, double s)
	{
		double x1, x2, w = 2, y1;

		if (useLast)
			{
				y1 = y2;
				useLast = false;
			}
		else
			{
				do
					{
						x1 = 2.0 * ranX() - 1;
						x2 = 2.0 * ranX() - 1;
						w = x1 * x1 + x2 * x2;
					}
				while (w >= 1.0);

				w = Math.sqrt((-2.0 * Math.log(w)) / w);
				y1 = x1 * w;
				y2 = x2 * w;
				useLast = true;
			}
		return (m + y1 * s);
	}

	public double ranX()
	{
		seed = seed * 1566083941 + 1;
		return seed * 2.328306e-10 + 0.5;
	}


}
