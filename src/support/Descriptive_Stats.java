package support;

import java.util.Arrays;

/**
 * JED class Descriptive_Stats: Computes descriptive stats for arrays of numbers. Copyright (C) 2012 Dr. Charles David
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
public class Descriptive_Stats
{
	public final double NORM = Math.sqrt(1.00 / (2.00 * Math.PI)), b = 1.4826;

	/* ************************************** CONSTRUCTOR ************************************************************************** */

	public Descriptive_Stats()
	{
		super();
	}

	/* ************************************** METHODS ************************************************************************** */

	/**
	 * Returns the mean of the double array sample.
	 * 
	 * @param sample a double array
	 * @return the mean
	 */
	public double get_mean(double[] sample)
	{
		long n = 0;
		double s = 0.0;

		for (double x : sample)
			{
				n++;
				s += x;
			}
		double mean = (s / n);
		return mean;
	}

	/**
	 * Returns the median of the double array sample.
	 * 
	 * @param sample a double array
	 * @return the median
	 */
	public double get_median(double[] sample)
	{
		Arrays.sort(sample);
		double median = 0;
		if (sample.length == 1) median = sample[0];
		else if (sample.length % 2 == 0) median = (sample[sample.length / 2] + sample[sample.length / 2 - 1]) / 2;
		else
			median = sample[sample.length / 2];
		return median;
	}

	/**
	 * Returns the median absolute deviation of the double array sample.
	 * 
	 * @param sample a double array
	 * @return the median absolute deviation
	 */
	public double get_mad(double[] sample)
	{
		double median = get_median(sample);
		double[] median_subtracted_array = new double[sample.length];
		int i = 0;
		for (double d : sample)
			{
				double value = Math.abs(d - median);
				median_subtracted_array[i] = value;
				i++;
			}
		double mad = get_median(median_subtracted_array) * b;
		return mad;
	}

	/**
	 * Returns the median absolute deviation scores for the double array sample. The parameter 'b' is used to specify the 0.75Q of the underlying distribution, 1.4826 for Normal.
	 * Here, Normal a Distribution is assumed. This calculation is NOT VALID for NON-NORMAL Distributions!
	 * 
	 * @param sample a double array
	 * @return the median absolute deviation scores
	 */
	public double[] get_mad_scores(double[] sample, double median, double mad)
	{
		double[] median_absolute_deviations = new double[sample.length];
		int i = 0;
		for (double d : sample)
			{
				double score = (d - median) / (mad);
				median_absolute_deviations[i] = score;
				i++;
			}
		return median_absolute_deviations;
	}

	/**
	 * Returns the sum of the squared deviations of the double array sample.
	 * 
	 * @param sample the double array
	 * @param mean   the mean of the sample
	 * @return the sum of squared deviations
	 */
	public double get_sum_of_squared_deviations(double[] sample, double mean)
	{
		double sum_squared_deviations = 0.0;
		for (double x : sample)
			{
				double delta = (x - mean) * (x - mean);
				sum_squared_deviations += delta;
			}
		return sum_squared_deviations;
	}

	/**
	 * Returns the sum of the cubed deviations of the double array sample.
	 * 
	 * @param sample the double array
	 * @param mean   the mean of the sample
	 * @return the sum of cubed deviations
	 */
	public double get_sum_of_cubed_deviations(double[] sample, double mean)
	{
		double sum_cubed_deviations = 0.0;
		for (double x : sample)
			{
				double delta = (x - mean) * (x - mean) * (x - mean);
				sum_cubed_deviations += delta;
			}
		return sum_cubed_deviations;
	}

	/**
	 * Returns the sum of the quartic deviations of the double array sample.
	 * 
	 * @param sample the double array
	 * @param mean   the mean of the sample
	 * @return the sum of quartic deviations
	 */
	public double get_sum_of_quartic_deviations(double[] sample, double mean)
	{
		double sum_quartic_deviations = 0.0;
		for (double x : sample)
			{
				double delta = (x - mean) * (x - mean) * (x - mean) * (x - mean);
				sum_quartic_deviations += delta;
			}
		return sum_quartic_deviations;
	}

	/**
	 * 
	 * Returns the variance of the double array sample.**
	 * 
	 * @param sample the double array*
	 * @param mean   the mean of the sample*@return
	 */

	public double get_variance(double[] sample, double mean)
	{
		int n = sample.length;
		double sum_squared_deviations = get_sum_of_squared_deviations(sample, mean);
		double variance = (sum_squared_deviations / (n - 1));
		return variance;
	}

	/**
	 * Returns the standard deviation of the double array sample.
	 * 
	 * @param sample the double array
	 * @return the standard deviation
	 */
	public double get_standard_deviation(double[] sample)
	{
		int n = sample.length;
		double mean = get_mean(sample);
		double sum_squared_deviations = get_sum_of_squared_deviations(sample, mean);
		double variance = (sum_squared_deviations / (n - 1));
		double std_dev = Math.sqrt(variance);
		return std_dev;
	}

	/**
	 * Returns the standard deviation of the sample given its mean.
	 * 
	 * @param sample
	 * @param mean
	 * @return
	 */
	public double get_standard_deviation(double[] sample, double mean)
	{
		int n = sample.length;
		double sum_squared_deviations = get_sum_of_squared_deviations(sample, mean);
		double variance = (sum_squared_deviations / (n - 1));
		double std_dev = Math.sqrt(variance);
		return std_dev;
	}

	/**
	 * Returns the standard deviation of the sample given its mean and the sum of squared deviations.
	 * 
	 * @param sample
	 * @param mean
	 * @param sum_sq_devs
	 * @return
	 */
	public double get_standard_deviation(double[] sample, double mean, double sum_sq_devs)
	{
		int n = sample.length;
		double sum_squared_deviations = sum_sq_devs;
		double variance = (sum_squared_deviations / (n - 1));
		double std_dev = Math.sqrt(variance);
		return std_dev;
	}

	/**
	 * Returns the Z scores of the sample given its mean and standard deviation.
	 * 
	 * @param sample
	 * @param mean
	 * @param std_dev
	 * @return
	 */
	public double[] get_Z_scores(double[] sample, double mean, double std_dev)
	{

		double[] z_scores = new double[sample.length];
		int i = 0;
		for (double d : sample)
			{
				double z = ((d - mean) / std_dev);
				z_scores[i] = z;
				i++;
			}
		return z_scores;
	}

	/**
	 * Returns the skew of the double array 'sample'.
	 * 
	 * @param sample the double array
	 * @return the skew
	 */
	public double get_skew(double[] sample)
	{
		int n = sample.length;
		double mean = get_mean(sample);
		double sum_squared_deviations = get_sum_of_squared_deviations(sample, mean);
		double variance = (sum_squared_deviations / (n - 1));
		double std_dev = Math.sqrt(variance);
		double sum_cubed_deviations = get_sum_of_cubed_deviations(sample, mean);
		double adjustment = Math.sqrt(n * (n - 1)) / (n - 2);
		double skew = adjustment * (sum_cubed_deviations) / (n * std_dev * std_dev * std_dev);

		return skew;
	}

	/**
	 * Returns the kurtosis of the double array 'sample'.
	 * 
	 * @param sample the double array
	 * @return the kurtosis
	 */
	public double get_kurtosis(double[] sample)
	{
		int n = sample.length;
		double mean = get_mean(sample);
		double sum_squared_deviations = get_sum_of_squared_deviations(sample, mean);
		double variance = (sum_squared_deviations / (n - 1));
		double sum_quartic_deviations = get_sum_of_quartic_deviations(sample, mean);
		double kurtosis = (sum_quartic_deviations) / (n * variance * variance);

		return kurtosis;
	}

	/**
	 * Returns a double array of probabilities given the data, its mean, and its standard deviation. This method assumes a normal distribution parameterized by mean and sigma. For
	 * non-Gaussian distributions, please use KDE
	 * 
	 * @param data
	 * @param mean
	 * @param sigma
	 * @return
	 */
	public double[] get_probabilities(double[] data, double mean, double sigma)
	{

		int n = 0;
		double[] probabilities = new double[data.length];
		for (double x : data)
			{

				double prob = 0;
				double coeff = (NORM / (sigma));
				double arg = (-0.50000000) * ((x - mean) / sigma) * ((x - mean) / sigma);
				prob = coeff * Math.exp(arg);
				probabilities[n] = prob;
				n++;
			}
		return probabilities;
	}
}
