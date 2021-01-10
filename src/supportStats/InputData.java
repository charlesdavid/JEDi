/*
 * 
 */
package supportStats;

import java.util.ArrayList;
import java.util.Collections;


/**
 *
 * @author jenny
 */
@SuppressWarnings("unused")
public class InputData
{
	OutputControl out = new OutputControl();

	double[] inverse;
	double[] doubleInverse;
	double[] xUntransform;
	double[] transformedZeroOne;
	double[] dz;

	int N = 0;
	int nPoints = 0;

	double minimumRaw = 0;
	double maximumRaw = 0;
	double minimumCalc = 0;
	double maximumCalc = 0;

	int nRightOutliers = 0;
	int nLeftOutliers = 0;
	int nPointsAdjust = 0;

	Chebyshev cheby = new Chebyshev();

	private InputParameters input;

	private boolean leftOutliers;
	private boolean rightOutliers;

	ArrayList<Double> rawData = new ArrayList<Double>();
	ArrayList<Double> transformedData = new ArrayList<Double>();


	public InputData(InputParameters input)
	{
		this.input = input;
	}

	public void setData(ArrayList<Double> data)
	{
		rawData = data;
		processData();
	}

	public boolean processData()
	{
		nPoints = input.integrationPoints;
		if (nPoints == -1)
			{
				// nPoints = (int) (200 + rawData.size() / 200.0);
				nPoints = rawData.size();
				// if (nPoints > 1500) nPoints = 1500;
			}
		Collections.sort(rawData);

		minimumRaw = rawData.get(0);
		maximumRaw = rawData.get(rawData.size() - 1);
		if (minimumRaw == maximumRaw)
			{
				out.error("All input data has the same value ", minimumRaw);
				return false;
			}
		identifyOutliers();
		if (!transformData())
			{
				return false;
			}
		setAdaptiveDz();
		cheby.initialize(doubleInverse, 2 * nPointsAdjust - 1);
		return true;
	}

	public void identifyOutliers()
	{
		double q1 = 0;
		double q3 = 0;

		int nValues = rawData.size();

		int middle = nValues / 2;
		int quarter = middle / 2;

		if (nValues % 2 == 0)
			{
				if (middle % 2 == 0)
					{
						q1 = (rawData.get(quarter - 1) + rawData.get(quarter)) / 2;
						q3 = (rawData.get(quarter + middle - 1) + rawData.get(quarter + middle)) / 2;
					}
				else
					{
						q1 = rawData.get(quarter);
						q3 = rawData.get(quarter + middle);
					}
			}
		else
			{
				if (middle % 2 == 0)
					{
						q1 = (rawData.get(quarter - 1) + rawData.get(quarter)) / 2;
						q3 = (rawData.get(quarter + middle) + rawData.get(quarter + middle + 1)) / 2;
					}
				else
					{
						q1 = rawData.get(quarter);
						q3 = rawData.get(quarter + middle) + 1;
					}
			}
		double iqr = input.outlierCutoff * (q3 - q1);
		double leftOutlier = q1 - iqr;
		double rightOutlier = q3 + iqr;

		if (input.upperBoundSpecified)
			{
				maximumCalc = input.upperBound;
			}
		else
			{
				double max = rawData.get(nValues - 1);
				maximumCalc = max + (max - rawData.get(rawData.size() - 5));
				if (maximumCalc > rightOutlier)
					{
						maximumCalc = rightOutlier;
						rightOutliers = true;
					}
			}

		if (input.lowerBoundSpecified)
			{
				minimumCalc = input.lowerBound;
			}
		else
			{
				double min = rawData.get(0);
				minimumCalc = min + (min - rawData.get(4));
				if (minimumCalc < leftOutlier)
					{
						minimumCalc = leftOutlier;
						leftOutliers = true;
					}
			}
	}

	public boolean transformData()
	{
		int nValues = rawData.size();
		ArrayList<Double> tempData = new ArrayList<Double>();

		for (int i = 0; i < rawData.size(); i++)
			{
				if (rawData.get(i) >= minimumCalc)
					{
						if (rawData.get(i) <= maximumCalc)
							{
								tempData.add(rawData.get(i));
							}
						else
							{
								nRightOutliers++;
							}
					}
				else
					{
						nLeftOutliers++;
					}
			}

		nValues = tempData.size();
		if (nValues == 0)
			{
				out.error("No data within specified boundaries");
				return false;
			}

		transformedZeroOne = new double[nValues];
		for (int i = 0; i < nValues; i++)
			{
				transformedData.add((2 * (tempData.get(i)) - maximumCalc - minimumCalc) / (maximumCalc - minimumCalc));
				transformedZeroOne[i] = (transformedData.get(i) + 1) / 2.0;
			}
		return true;
	}

	public void setAdaptiveDz()
	{
		ArrayList<Double> dzVector = new ArrayList<Double>();
		N = transformedData.size();

		double dzMax = 2.0 / (nPoints - 1);

		int skip = N / (nPoints - 1);
		if (skip == 0) skip = 1;

		double last = -1.0;
		double next;

		for (int b = skip; b <= (N + skip); b += skip)
			{
				if (b >= (N))
					{
						next = transformedData.get(N - 1);
					}
				else
					{
						next = transformedData.get(b);
					}
				double test = next - last;
				double difference = Math.abs(test);
				if (difference > dzMax)
					{
						double steps = difference / dzMax;
						int iSteps = (int) steps;
						for (int k = 0; k < (iSteps + 1); k++)
							{
								dzVector.add(difference / (iSteps + 1));
							}
					}
				else
					{
						dzVector.add(difference);
					}
				last = next;
			}

		int dzSize = dzVector.size();
		inverse = new double[dzSize + 1];
		inverse[0] = 0;
		for (int j = 1; j <= dzSize; j++)
			{
				inverse[j] = inverse[j - 1] + dzVector.get(j - 1) / 2.0;
			}

		double difference = 1.0 - inverse[dzSize];
		if (difference > dzMax)
			{
				double steps = (difference) / dzMax;
				int iSteps = (int) steps;
				for (int k = 0; k < 2 * (iSteps + 1); k++)
					{
						dzVector.add(difference / (iSteps + 1));
					}
			}
		else
			{
				dzVector.add(difference);
			}
		dzSize = dzVector.size();
		dz = new double[dzSize];
		inverse = new double[dzSize];
		doubleInverse = new double[2 * dzSize - 1];
		inverse[0] = dzVector.get(0) / 2.0;
		dz[0] = dzVector.get(0) / 2.0;
		int count = 0;
		for (int j = 1; j < dzSize; j++)
			{
				dz[j] = dzVector.get(j) / 2.0;
				inverse[j] = inverse[j - 1] + dzVector.get(j) / 2.0;
				doubleInverse[count] = inverse[j - 1];
				doubleInverse[count + 1] = (inverse[j - 1] + inverse[j]) / 2.0;
				count += 2;
			}
		doubleInverse[count] = (inverse[dzSize - 1] + 1.0) / 2.0;

		xUntransform = new double[2 * dzSize - 1];

		for (int i = 0; i < 2 * dzSize - 1; i++)
			{
				xUntransform[i] = doubleInverse[i] * 2.0 - 1;
			}
		for (int i = 0; i < 2 * dzSize - 1; i++)
			{
				xUntransform[i] = (maximumCalc - minimumCalc) * xUntransform[i] + minimumCalc + maximumCalc;
				xUntransform[i] /= 2;
			}
		nPointsAdjust = dzSize;
	}
}
