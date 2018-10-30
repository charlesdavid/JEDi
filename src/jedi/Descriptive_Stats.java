package jedi;

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

	public static double mean, variance, std_dev, sum_squared_deviations;
	public static final double NORM = Math.sqrt(1.00 / (2.00 * Math.PI));

	/**
	 * Returns the mean of the double array sample.
	 * 
	 * @param sample
	 *            a double array
	 * @return the mean
	 */
	public static double get_mean(double[] sample)
		{

			long n = 0;
			mean = 0;
			double s = 0.0;

			for (double x : sample)
				{
					n++;
					s += x;
				}
			return mean = (s / n);
		}

	/**
	 * Returns the sum of the squared deviation of the double array sample.
	 * 
	 * @param sample
	 *            the double array
	 * @param mean
	 *            the mean of the sample
	 * @return the sum of squared deviations
	 */
	public static double get_sum_of_squared_deviations(double[] sample, double mean)
		{
			sum_squared_deviations = 0.0;
			for (double x : sample)
				{
					double delta = (x - mean) * (x - mean);
					sum_squared_deviations += delta;
				}
			return sum_squared_deviations;
		}

	/**
	 * Returns the variance of the double array sample.
	 * 
	 * @param sample
	 *            the double array
	 * @param mean
	 *            the mean of the sample
	 * @return
	 */
	public static double get_variance(double[] sample, double mean)
		{

			int n = sample.length;
			get_sum_of_squared_deviations(sample, mean);
			variance = (sum_squared_deviations / (n - 1));
			return variance;
		}

	/**
	 * Returns the standard deviation of the double array sample.
	 * 
	 * @param sample
	 *            the double array
	 * @return the standard deviation
	 */
	public static double get_standard_deviation(double[] sample)
		{

			int n = sample.length;
			mean = get_mean(sample);
			get_sum_of_squared_deviations(sample, mean);
			variance = (sum_squared_deviations / (n - 1));
			std_dev = Math.sqrt(variance);
			return std_dev;
		}

	/**
	 * Returns the standard deviation of the sample given its mean.
	 * 
	 * @param sample
	 * @param mean
	 * @return
	 */
	public static double get_standard_deviation(double[] sample, double mean)
		{

			int n = sample.length;
			get_sum_of_squared_deviations(sample, mean);
			variance = (sum_squared_deviations / (n - 1));
			std_dev = Math.sqrt(variance);
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
	public static double get_standard_deviation(double[] sample, double mean, double sum_sq_devs)
		{

			int n = sample.length;
			sum_squared_deviations = sum_sq_devs;
			variance = (sum_squared_deviations / (n - 1));
			std_dev = Math.sqrt(variance);
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
	public static double[] get_Z_scores(double[] sample, double mean, double std_dev)
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
	 * Returns a double array of probabilities given the data, its mean, and its standard deviation.
	 * This method assumes a normal distribution parameterized by mean and sigma.
	 * For non-Gaussian distributions, please use KDE
	 * 
	 * @param data
	 * @param mean
	 * @param sigma
	 * @return
	 */
	public static double[] get_probabilities(double[] data, double mean, double sigma)
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
