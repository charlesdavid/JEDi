/*
 * Copyright (c) 2014, Massachusetts Institute of Technology Released under the BSD 2-Clause License http://opensource.org/licenses/BSD-2-Clause
 */

package jedi;

import bits.fft.FastCosineTransform2d;
import jed.Pots;

/**
 * So you've got some 2d points and you want to convolve them with a gaussian kernel to generate a pdf. What bandwidth of kernel do you use? Well, there's some criteria, like MSIE
 * or some such, that you should usually try to minimize. That is what this class does. Give it some points and it will give you the optimal kernel bandwidth.
 * <p>
 * Currently limited to diagonal bandwidth matrices (Gaussians that are symmetric about each axis).
 * <p>
 * This is based on matlab code by Zdravko Botev:
 * <p>
 * Z. I. Botev, J. F. Grotowski and D. P. Kroese "KERNEL DENSITY ESTIMATION VIA DIFFUSION", Submitted to the Annals of Statistics, 2009
 * 
 * @author Philip DeCamp
 */
public class DiagonalBandwidthSelector2d
	{

		/**
		 * @param points    Set of 2d points, closely packed: [x0,y0,x1,y1...]. <code>points.length > off + numPoints * 2</code>
		 * @param off       Offset into points array.
		 * @param numPoints Number of points.
		 * @param bounds    A length-4 array containing the bounds [x0,y0,x1,y1] of the region of interest. (Optional)
		 * @return A 2x2 matrix, [exx, exy; exy, eyy] containing the bandwidth matrix. Will always be a diagonal matrix.
		 * @throws MathException If method fails to estimate optimal bandwidth.
		 */
		static double[] bandwidths;

		public DiagonalBandwidthSelector2d()
			{

			}

		public static double[] computeBandwidth(double[] points, int off, int numPoints, double[] bounds) throws MathException
			{
				return computeBandwidth(points, off, numPoints, bounds, 1 << 8);
			}

		/**
		 * @param points    Set of 2d points, closely packed: [x0,y0,x1,y1...]. <code>points.length > off + numPoints * 2</code>
		 * @param off       Offset into points array.
		 * @param numPoints Number of points.
		 * @param bounds    A length-4 array containing the bounds [x0,y0,x1,y1] of the region of interest. [Optional]
		 * @param quant     Determines level of quantization of point positions. Default is 256.
		 * @return A 2x2 matrix, [exx, exy; exy, eyy] containing the bandwidth matrix. Will always be diagonal.
		 * @throws MathException If method fails to estimate optimal bandwidth.
		 */
		public static double[] computeBandwidth(double[] points, int off, int numPoints, double[] bounds, int quant) throws MathException
			{
				quant = Pots.ceilPot(quant);

				if (bounds == null)
					{
						bounds = computeBounds(points, off, numPoints, 0.25);
					}

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
						tStar = FZero.findZeroIn(func, 0.0, 0.1); // was 0.1
					} catch (MathException ignore)
					{
					}

				if (Double.isNaN(tStar))
					{
						try
							{
								tStar = FZero.findZeroNear(func, 0.15);
							} catch (MathException ex)
							{
								// throw new MathException("Failed to minimize objective function. Possibly not enough points.");
								bandwidths = new double[] { 0, 0, 0, 0 };
								bandwidths[0] = ((bounds[2] - bounds[0]) / 256);
								bandwidths[3] = ((bounds[3] - bounds[1]) / 256);
								System.err.println("WARNING: Kernel Bandwidths could not be optimized: Possibly not enough points were used.");
								System.err.println("Kernel Bandwidths have been set to: " + bandwidths[0] + ", " + bandwidths[3]);
								return bandwidths;
							}
					}

				double p02 = func.sumFunc(0, 2, tStar);
				double p20 = func.sumFunc(2, 0, tStar);
				double p11 = func.sumFunc(1, 1, tStar);

				double ty = Math.pow(Math.pow(p02, 0.75) / (4.0 * Math.PI * numPoints * Math.pow(p20, 0.75) * (p11 + Math.sqrt(p20 * p02))), 1.0 / 3.0);
				double tx = Math.pow(Math.pow(p20, 0.75) / (4.0 * Math.PI * numPoints * Math.pow(p02, 0.75) * (p11 + Math.sqrt(p20 * p02))), 1.0 / 3.0);

				bandwidths = new double[] { 0, 0, 0, 0 };
				bandwidths[0] = Math.sqrt(tx) * (bounds[2] - bounds[0]);
				bandwidths[3] = Math.sqrt(ty) * (bounds[3] - bounds[1]);
				return bandwidths;
			}

		public static double[] get_Bandwidths()
			{
				return bandwidths;
			}

		private static final double PI2 = Math.PI * 2.0;
		private static final double PIPI = Math.PI * Math.PI;
		private static final double INV_SQRT_PI2 = 1.0 / Math.sqrt(PI2);

		private static double[] computeBounds(double[] points, int off, int len, double margin)
			{
				double x0 = Double.POSITIVE_INFINITY;
				double x1 = Double.NEGATIVE_INFINITY;
				double y0 = Double.POSITIVE_INFINITY;
				double y1 = Double.NEGATIVE_INFINITY;

				for (int i = 0; i < len; i++)
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

		private static double[] hist(double[] points, int off, int len, double[] bounds, int quant)
			{
				double[] ret = new double[quant * quant];

				double addX = -bounds[0];
				double addY = -bounds[1];

				double scaleX = (1.0 - 1E-10) * quant / (bounds[2] - bounds[0]);
				double scaleY = (1.0 - 1E-10) * quant / (bounds[3] - bounds[1]);
				double weight = 1.0 / len;

				for (int i = 0; i < len; i++)
					{
						// double v = points[off + i];
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

		private static class EvolveFunc implements Function11
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

						// return Math.pow(-1.0, s0 + s1) * sum * Math.pow(Math.PI, 2.0 * (s0 + s1));
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

							} else
							{
								return psi(s0, s1, t);
							}
					}

				static double k(int s)
					{
						double prod = 1.0;

						for (int i = 3; i <= 2 * s - 1; i += 2)
							{
								prod *= i;
							}

						// return Math.pow(-1.0, s) * prod * INV_SQRT_PI2;
						return ((s & 1) * -2 + 1) * prod * INV_SQRT_PI2;
					}

			}

	}
