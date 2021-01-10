package supportKDE;

import bits.fft.FastCosineTransform2d;
import supportIO.Input_Parameters;

/**
 * Class to transform 2D points into a probability distribution function.
 * 
 * This class performs that task with a bivariate kernel density estimation algorithm.
 * 
 * Essentially, each point is convolved with a kernel. The kernel is the normal function.
 * 
 * The main difficulty here is choosing the bandwidth of the kernel...
 * 
 * This implementation computes a table of the KDE PDF within some range.
 * 
 * Values provided by this function are not exact, but a based on bilinear interpolation from this table.
 *
 */
public class KDE_2D
{
	final double kernel_resolution = Input_Parameters.KDE_RESOLUTION;
	final double margin = Input_Parameters.KDE_MARGIN;
	final double kde_cell = Input_Parameters.KDE_CELL;
	final double PI2 = Math.PI * 2.0;
	final double PIPI = Math.PI * Math.PI;
	final double INV_SQRT_PI2 = 1.0 / Math.sqrt(PI2);
	final int numberPoints, offset;
	Double cellSize;
	double[] data, bounds, bandwidths;
	BiLinear_Sampler bs;
	Function21 kernel;

	/* *************************************************************************************************************************** */

	public KDE_2D(double[] points, int off, int numPoints, double[] bounds, Double cellsize, Function21 kernel)
	{
		this.data = points;
		this.offset = off;
		this.numberPoints = numPoints;
		this.bounds = bounds;
		this.cellSize = cellsize;
		this.kernel = kernel;
	}

	/* *************************************************************************************************************************** */

	public BiLinear_Sampler compute()
	{
		if (bounds == null)
			{
				bounds = computeBounds();
			}
		if (kernel == null)
			{
				bandwidths = getBandwidths(bounds);
				// bandwidths = getOptimalBandwidths(data, offset, 100, bounds, 1 << 8); // Method never completes...
				kernel = selectKernel(bandwidths);
			}
		if (cellSize == null)
			{
				double maxDim = Math.max(bounds[2] - bounds[0], bounds[3] - bounds[1]);
				cellSize = (maxDim / kde_cell);
			}
		// Construct table of KDE values along a grid.
		// Note that the grid may exceed bounds.
		// Samples of the PDF are taken at each grid corner.
		final int cols = (int) Math.ceil((bounds[2] - bounds[0]) / cellSize) + 1;
		final int rows = (int) Math.ceil((bounds[3] - bounds[1]) / cellSize) + 1;
		final double[] table = new double[cols * rows];

		// Iterate through each corner of cell.
		for (int y = 0; y < rows; y++)
			{
				for (int x = 0; x < cols; x++)
					{
						final double px = bounds[0] + cellSize * x;
						final double py = bounds[1] + cellSize * y;
						double sum = 0.0;
						// Compute contribution of each point .
						for (int p = 0; p < numberPoints; p++)
							{
								double dx = data[offset + p * 2] - px;
								double dy = data[offset + p * 2 + 1] - py;
								sum += kernel.apply(dx, dy);
							}
						// Insert value into table.
						table[x + y * cols] = sum;
					}
			}
		// Compute volume of grid.
		double tableVolume = 0.0;
		double cellArea = cellSize * cellSize;

		for (int y = 0; y < rows - 1; y++)
			{
				final int y0 = y * cols;
				final int y1 = y0 + cols;

				for (int x = 0; x < cols - 1; x++)
					{
						double v0 = table[x + y0];
						double v1 = table[x + y1];
						double v2 = table[x + 1 + y0];
						double v3 = table[x + 1 + y1];

						tableVolume += cellArea * (v0 + v1 + v2 + v3) * 0.25;
					}
			}

		// Normalize table.
		for (int i = 0; i < table.length; i++)
			{
				table[i] /= tableVolume;
			}
		mKernel = kernel;
		bs = new BiLinear_Sampler(table, cols, rows, bounds[0], bounds[1], bounds[0] + cellSize * cols, bounds[1] + cellSize * rows);
		return bs;
	}

	public double[] getBandwidths(double[] bounds)
	{
		double[] band_widths = new double[] { 0, 0, 0, 0 };
		band_widths[0] = ((bounds[2] - bounds[0]) / kernel_resolution);
		band_widths[3] = ((bounds[3] - bounds[1]) / kernel_resolution);
		return band_widths;
	}

	public Function21 selectKernel(double[] bandWidth2x2)
	{
		return Gaussian2.fromSigma(bandWidth2x2[0], bandWidth2x2[3], bandWidth2x2[1]);
	}

	public double[] computeBounds()
	{
		double x0 = -1d;
		double x1 = 1d;
		double y0 = -1d;
		double y1 = 1d;
		double mx = (x1 - x0) * margin;
		double my = (y1 - y0) * margin;

		bounds = new double[] { x0 - mx, y0 - my, x1 + mx, y1 + my };
		return bounds;
	}

	public double[] computeBoundsReal(double minX, double maxX, double minY, double maxY)
	{
		double x0 = minX;
		double x1 = maxX;
		double y0 = minY;
		double y1 = maxY;
		double mx = (x1 - x0) * margin;
		double my = (y1 - y0) * margin;

		bounds = new double[] { x0 - mx, y0 - my, x1 + mx, y1 + my };
		return bounds;
	}

	public double[] getOptimalBandwidths(double[] points, int off, int numPoints, double[] bounds, int quant)
	{
		quant = Pots.ceilPot(quant);

		double[] hist = hist(points, off, numPoints, bounds, quant);
		double[] freq = new double[quant * quant];

		FastCosineTransform2d trans = new FastCosineTransform2d(quant);
		trans.apply(hist, 0, false, freq, 0);

		double[] bigI = new double[quant];
		double[] bigA = new double[quant * quant];

		for (int x = 0; x < quant; x++)
			{
				bigI[x] = x * x;

				for (int y = 0; y < quant; y++)
					{
						double n = freq[y + x * quant];
						bigA[y + x * quant] = n * n;
					}
			}

		EvolveFunc func = new EvolveFunc(quant, numPoints, bigI, bigA);
		double tStar = Double.NaN;

		try
			{
				tStar = FZero.findZeroIn(func, 0.0, 0.1); // original code value is: 0.1
			}
		catch (MathException ignore)
			{
			}

		if (Double.isNaN(tStar))
			{
				try
					{
						tStar = FZero.findZeroNear(func, 0.15);
					}
				catch (MathException ex)
					{
						System.err.println("WARNING: Kernel Bandwidths could not be optimized: Possibly not enough points were used.");
						System.err.println("Using specified Kernel Resolution to compute bandwidths.");
						return bandwidths = getBandwidths(bounds);
					}
			}
		double p02 = func.sumFunc(0, 2, tStar);
		double p20 = func.sumFunc(2, 0, tStar);
		double p11 = func.sumFunc(1, 1, tStar);

		double ty = Math.pow(Math.pow(p02, 0.75) / (4.0 * Math.PI * numPoints * Math.pow(p20, 0.75) * (p11 + Math.sqrt(p20 * p02))), 1.0 / 3.0);
		double tx = Math.pow(Math.pow(p20, 0.75) / (4.0 * Math.PI * numPoints * Math.pow(p02, 0.75) * (p11 + Math.sqrt(p20 * p02))), 1.0 / 3.0);

		double[] band_widths = new double[] { 0, 0, 0, 0 };
		band_widths[0] = Math.sqrt(tx) * (bounds[2] - bounds[0]);
		band_widths[3] = Math.sqrt(ty) * (bounds[3] - bounds[1]);
		return band_widths;
	}

	private double[] hist(double[] points, int off, int len, double[] bounds, int quant)
	{
		double[] ret = new double[quant * quant];

		double addX = -bounds[0];
		double addY = -bounds[1];

		double scaleX = (1.0 - 1E-10) * quant / (bounds[2] - bounds[0]);
		double scaleY = (1.0 - 1E-10) * quant / (bounds[3] - bounds[1]);
		double weight = 1.0 / len;

		for (int i = 0; i < len; i++)
			{
				int p0 = (int) ((points[off + i * 2] + addX) * scaleX);
				int p1 = (int) ((points[off + i * 2 + 1] + addY) * scaleY);
				if (p0 < 0 || p1 < 0 || p0 >= quant || p1 >= quant)
					{
						continue;
					}

				ret[p0 + p1 * quant] += weight;
			}
		return ret;
	}

	private class EvolveFunc implements Function11
	{
		private final int mDim;
		private final int mLen;

		private final double[] mBigI;
		private final double[] mBigA;

		private final double[] mWorkA;
		private final double[] mWorkB;

		EvolveFunc(int dim, int len, double[] bigI, double[] bigA)
		{
			mDim = dim;
			mLen = len;

			mBigI = bigI;
			mBigA = bigA;

			mWorkA = new double[dim];
			mWorkB = new double[dim];
		}

		@Override
		public double apply(double t)
		{
			double sf = sumFunc(0, 2, t) + sumFunc(2, 0, t) + 2.0 * sumFunc(1, 1, t);

			double time = Math.pow(PI2 * mLen * sf, -1.0 / 3.0);
			return t - (t - time) / time;
		}

		double psi(int s0, int s1, double time)
		{
			final int dim = mDim;
			final double[] work0 = mWorkA;
			final double[] work1 = mWorkB;
			final double[] bigI = mBigI;
			final double[] bigA = mBigA;

			for (int i = 0; i < dim; i++)
				{
					double ww = Math.exp(-bigI[i] * PIPI * time) * (i == 0 ? 1.0 : 0.5);
					work0[i] = ww * Math.pow(bigI[i], s0);
					work1[i] = ww * Math.pow(bigI[i], s1);
				}

			double sum = 0.0;

			for (int i = 0; i < dim; i++)
				{
					double vy = work1[i];
					for (int j = 0; j < dim; j++)
						{
							sum += vy * bigA[i + j * dim] * work0[j];
						}
				}
			return (((s0 + s1) & 1) * -2 + 1) * sum * Math.pow(Math.PI, 2.0 * (s0 + s1));
		}

		double sumFunc(int s0, int s1, double t)
		{
			if (s0 + s1 <= 4)
				{
					double vs = sumFunc(s0 + 1, s1, t) + sumFunc(s0, s1 + 1, t);
					// double vc = (1 + 1.0 / Math.pow(2.0, s0 + s1 + 1)) / 3.0;
					double vc = (1.0 + 1.0 / (1 << (s0 + s1 + 1))) / 3.0;
					double time = Math.pow(-2.0 * vc * k(s0) * k(s1) / mLen / vs, 1.0 / (2.0 + s0 + s1));
					return psi(s0, s1, time);
				}
			else
				{
					return psi(s0, s1, t);
				}
		}

		double k(int s)
		{
			double prod = 1.0;

			for (int i = 3; i <= 2 * s - 1; i += 2)
				{
					prod *= i;
				}
			return ((s & 1) * -2 + 1) * prod * INV_SQRT_PI2;
		}
	}

	public void setBounds(double[] bounds1x4)
	{
		this.bounds = bounds1x4;
	}

	private Function21 mKernel;

	public double getProb(double x, double y)
	{
		return bs.apply(x, y);
	}

	public Function21 kernel()
	{
		return mKernel;
	}

	public double[] get_band_Widths()
	{
		return bandwidths;
	}

	public double get_cell_Size()
	{
		return cellSize;
	}

	public double[] getBounds()
	{
		return bounds;
	}
}
