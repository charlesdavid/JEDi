/*
 *
 */

package supportStats;

import java.util.ArrayList;
import java.util.Collections;
//import java.util.List;

/**
 *
 * @author jenny
 */
public abstract class Score
{
	double targetScore;
	double minimumScore;
	double maximumScore;

	double SURD;
	double QZVariance = 0;

	protected double likelihood;
	protected int spacing;
	protected ArrayList<Double> scores = new ArrayList<Double>();
	protected ArrayList<Double> SURDs = new ArrayList<Double>();
	protected ArrayList<Integer> indices = new ArrayList<Integer>();

	public double getLikelihood()
	{
		return likelihood;
	}

	public ArrayList<Integer> setIndices(int N, int p)
	{
		return null;
	}

	public ArrayList<Integer> getIndices(int N, int p)
	{
		return null;
	}

	public double calculateScore(double r[], int N, int p)
	{
		return 0;
	}

	public double getTargetScore(double SURD)
	{

		int index = -Collections.binarySearch(SURDs, SURD / 100);

		if (index == SURDs.size())
			{
				return scores.get(index - 1);
			}
		if (index == 0)
			{
				return scores.get(index);
			}

		if (index >= scores.size())
			{
				return scores.get(scores.size() - 1);
			}
		else
			{
				double E1 = scores.get(index - 1);
				double E2 = scores.get(index);
				double P1 = SURDs.get(index - 1);
				double P2 = SURDs.get(index);
				double E = E1 + (SURD / 100 - P1) * (E2 - E1) / (P2 - P1);
				return E;
			}
	}

	double getConfidence(double score)
	{

		int index = Collections.binarySearch(scores, score / 100);

		if (index == scores.size())
			{
				return SURDs.get(index - 1);
			}
		if (index == 0)
			{
				return SURDs.get(index);
			}

		double E1 = scores.get(index - 1);
		double E2 = scores.get(index);
		double P1 = SURDs.get(index - 1);
		double P2 = SURDs.get(index);
		double P = P1 + (score - E1) * (P2 - P1) / (E2 - E1);
		return P * 100;

	}
}
