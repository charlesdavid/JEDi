/*
 * Copyright (c) 2014, Massachusetts Institute of Technology Released under the BSD 2-Clause License http://opensource.org/licenses/BSD-2-Clause
 */
package supportKDE;

/**
 * So you've got some 2d points and you want to transform that into a probability distribution function. This class performs that task with a bivariate kernel density estimation
 * algorithm. Essentially, each point is convolved with a kernel. Usually, the kernel is the normal function. The main difficulty here is selecting the bandwidth of the kernel.
 * This class will do so automatically. For more information, see DiagonalBandwidthSelector2d. This implementation precomputes a table of the KDE PDF within some range. Values
 * provided by this function are not exact, but a based on bilinear interpolation from this table.
 *
 * @author Philip DeCamp
 */
public class KernelDensityEstimate2d extends BilinearSampler implements Function21
{
	/**
	 * Provides a 2D Kernel Density Estimate of a set of 2D points.
	 *
	 * @param points    Array containing points: [x0, y0, x1, y1...]
	 * @param off       Offset into point array.
	 * @param numPoints Number of points to use.
	 * @return New KernelDensityEstimate of the given data points.
	 * @throws Exception
	 */
	public static KernelDensityEstimate2d compute(double[] points, int off, int numPoints) throws Exception
	{
		return compute(points, off, numPoints, null, null, null);
	}

	/**
	 * Generates a 2D Kernel Density Estimate of a set of 2D points.
	 *
	 * @param points    Array containing points: [x0, y0, x1, y1...]
	 * @param off       Offset into point array.
	 * @param numPoints Number of points to use.
	 * @param bounds    Bounds over which the PDF table will be computed. (Optional)
	 * @param cellSize  Size of the side of each cell in table: IE, the sampling interval. (Optional)<br>
	 *                  Use a low value for more approximate results, but faster initialization and less memory use. <br>
	 *                  Use a high value for more accurate results, but slower initialization and more memory use. <br>
	 *                  Use <code>null</code> to select a default cell size.
	 * @param kernel    Kernel function (Optional) <br>
	 *                  Set to <code>null</code> to select a kernel automatically.
	 * @return New KernelDensityEstimate of the given data points.
	 * @throws Exception
	 */
	public static KernelDensityEstimate2d compute(double[] points, int off, int numPoints, double[] bounds, Double cellSize, Function21 kernel)
	{
		if (bounds == null)
			{
				bounds = selectBounds(points, off, numPoints, 0.25);
			}
		if (kernel == null)
			{
				double[] bandWidth = selectBandwidth(points, off, numPoints, bounds);
				kernel = selectKernel(bandWidth);
			}
		if (cellSize == null)
			{
				double maxDim = Math.max(bounds[2] - bounds[0], bounds[3] - bounds[1]);
				// cellSize = maxDim / 256; // original value
				cellSize = maxDim / 64;
			}
		// Construct table of KDE values along a grid.
		// Note that the grid may exceed bounds.
		// Samples of the PDF are taken at each grid corner.
		final double d = cellSize;
		final int cols = (int) Math.ceil((bounds[2] - bounds[0]) / d) + 1;
		final int rows = (int) Math.ceil((bounds[3] - bounds[1]) / d) + 1;
		final double[] table = new double[cols * rows];

		// Iterate through each corner of cell.
		for (int y = 0; y < rows; y++)
			{
				for (int x = 0; x < cols; x++)
					{
						final double px = bounds[0] + d * x;
						final double py = bounds[1] + d * y;
						double sum = 0.0;
						// Compute contribution of each point .
						for (int p = 0; p < numPoints; p++)
							{
								double dx = points[off + p * 2] - px;
								double dy = points[off + p * 2 + 1] - py;
								sum += kernel.apply(dx, dy);
								// System.out.println("sum = " + sum);
							}
						// Insert value into table.
						table[x + y * cols] = sum;
					}
			}
		// Compute volume of grid.
		double tableVolume = 0.0;
		double cellArea = d * d;

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

		return new KernelDensityEstimate2d(table, cols, rows, bounds[0], bounds[1], bounds[0] + d * cols, bounds[1] + d * rows, kernel);
	}

	/**
	 * Selects the bandwidth for a set of points. Pretty much just for convenience, so you won't have to look at DiagonalBandwidthSelector2d.
	 *
	 * @param points    Array containing points: [x0, y0, x1, y1...]
	 * @param off       Offset into point array.
	 * @param numPoints Number of points to use.
	 * @param bounds2x2 Bounds being used for KDE function. (Optional)
	 * @return 2x2 bandwidth matrix
	 * @throws Exception
	 */
	public static double[] selectBandwidth(double[] points, int off, int numPoints, double[] bounds2x2)
	{
		if (bounds2x2 == null) bounds2x2 = selectBounds(points, off, numPoints, 0.25);
		double[] bandWidth = DiagonalBandwidthSelector2d.computeBandwidth(points, off, numPoints, bounds2x2);
		return bandWidth;
	}

	/**
	 * Creates a kernel. It's a normal function.
	 *
	 * @param bandWidth2x2 2x2 Bandwidth matrix
	 * @return a normal kernel function.
	 */
	public static Function21 selectKernel(double[] bandWidth2x2)
	{
		return Gaussian2.fromSigma(bandWidth2x2[0], bandWidth2x2[3], bandWidth2x2[1]);
	}

	/**
	 * Selects bounds for a set of points. Recommended margin is 0.25.
	 *
	 * @param points    Array containing points: [x0, y0, x1, y1...]
	 * @param off       Offset into point array.
	 * @param numPoints Number of points to use.
	 * @param margin    Margin around points, in units of the span of the points.
	 * @return 4x1 matrix [minX, minY, maxX, maxY]
	 */
	public static double[] selectBounds(double[] points, int off, int numPoints, double margin)
	{
		double x0 = Double.POSITIVE_INFINITY;
		double x1 = Double.NEGATIVE_INFINITY;
		double y0 = Double.POSITIVE_INFINITY;
		double y1 = Double.NEGATIVE_INFINITY;

		for (int i = 0; i < numPoints; i++)
			{
				double v = points[i * 2 + off];

				if (v < x0)
					{
						x0 = v;
					}
				if (v > x1)
					{
						x1 = v;
					}

				v = points[i * 2 + 1 + off];

				if (v < y0)
					{
						y0 = v;
					}
				if (v > y1)
					{
						y1 = v;
					}
			}

		double mx = (x1 - x0) * margin;
		double my = (y1 - y0) * margin;

		return new double[] { x0 - mx, y0 - my, x1 + mx, y1 + my };
	}

	private final Function21 mKernel;

	private KernelDensityEstimate2d(double[] tableRef, int width, int height, double minX, double minY, double maxX, double maxY, Function21 kernel)
	{
		super(tableRef, width, height, minX, minY, maxX, maxY);
		mKernel = kernel;
	}

	/**
	 * @return Kernel used by this KernelDensityEstimate to generate PDF table.
	 */
	public Function21 kernel()
	{
		return mKernel;
	}
}
