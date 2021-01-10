/*
 * Copyright (c) 2014, Massachusetts Institute of Technology Released under the BSD 2-Clause License http://opensource.org/licenses/BSD-2-Clause
 */
package supportKDE;

/**
 * Port of the Matlab fzero function. Finds zeros in univariate functions.
 *
 * @author Philip DeCamp
 */
public class FZero
	{

		public static final int CODE_OK = 0;

		private static final double SQRT2 = Math.sqrt(2.0);
		private static final double TOL = 2.2204E-15;


		/**
		 * @param func A function to search.
		 * @param x0   Start value
		 * @return Value of x near a sign change.
		 * @throws MathException If fails to find a sign change.
		 */
		public static double findZeroNear(Function11 func, double x0) throws MathException
			{
				return findZeroNear(func, x0, TOL);
			}

		/**
		 * @param func A function to search.
		 * @param x    Start value
		 * @param tol  Tolerance of result
		 * @return Value of x near sign change.
		 * @throws MathException If fails to find a sign change.
		 */
		public static double findZeroNear(Function11 func, double x, double tol) throws MathException
			{
				double fx = applyFunc(func, x);
				if (fx == 0)
					{
						return x;
					}

				double dx;

				if (x == 0.0)
					{
						dx = 1.0 / 50.0;
					} else
					{
						dx = x / 50.0;
					}

				// Find change of sign.
				double a = x;
				double fa = fx;
				double b = x;
				double fb = fx;

				while ((fa > 0.0) == (fb > 0.0))
					{
						dx *= SQRT2;
						a = x - dx;
						fa = applyFunc(func, a);

						// Check for different sign.
						if ((fa > 0.0) != (fb > 0.0))
							{
								break;
							}

						b = x + dx;
						fb = applyFunc(func, b);
					}

				return refine(func, a, fa, b, fb, tol);
			}

		/**
		 * @param func A function to search
		 * @param x0   Start of interval to search
		 * @param x1   End of interval to search
		 * @return Value of x within interval [x0,x1] where sign changes.
		 * @throws MathException If fails to find a sign change, or if <code>f(x0)</code> has same sign as <code>f(x1)</code>.
		 */
		public static double findZeroIn(Function11 func, double x0, double x1) throws MathException
			{
				return findZeroIn(func, x0, x1, TOL);
			}

		/**
		 * @param func A function to search
		 * @param x0   Start of interval to search
		 * @param x1   End of interval to search
		 * @param tol  Tolerance of answer.
		 * @return Value of x within interval [x0,x1] where sign changes.
		 * @throws MathException If fails to find a sign change, or if <code>f(x0)</code> has same sign as <code>f(x1)</code>.
		 */
		public static double findZeroIn(Function11 func, double x0, double x1, double tol) throws MathException
			{
				double a = x0;
				double b = x1;

				double fa = applyFunc(func, x0);
				double fb = applyFunc(func, x1);

				if (fa == 0.0)
					{
						return a;
					}
				if (fb == 0.0)
					{
						return b;
					}
				if ((fa > 0) == (fb > 0))
					{
						throw new MathException("The function values at the interval endpoints must differ in sign.");
					}

				return refine(func, a, fa, b, fb, tol);
			}


		private static double refine(Function11 func, double a, double fa, double b, double fb, double tol) throws MathException
			{
				double c = a;
				double fc = fa;
				double d = b - a;
				double e = d;

				while (fb != 0.0 && a != b)
					{
						// Ensure that b is the best result so far, a is the previous value of b,
						// and c is on the opposite side of the zero from b.
						if ((fb > 0.0) == (fc > 0.0))
							{
								c = a;
								fc = fa;
								d = b - a;
								e = d;
							}

						if (Math.abs(fc) < Math.abs(fb))
							{
								a = b;
								b = c;
								c = a;
								fa = fb;
								fb = fc;
								fc = fa;
							}

						// Convergence test and possible exit.
						double m = 0.5 * (c - b);
						double toler = 2.0 * tol * Math.max(Math.abs(b), 1.0);

						if (Math.abs(m) <= toler || fb == 0.0)
							{
								break;
							}

						// Choose bisection or interpolation
						if ((Math.abs(e) < toler) || (Math.abs(fa) <= Math.abs(fb)))
							{
								// Bisection.
								d = m;
								e = m;

							} else
							{
								// Interpolation
								double s = fb / fa;
								double p, q;

								if (a == c)
									{
										// Linear interpolation.
										p = 2.0 * m * s;
										q = 1.0 - s;

									} else
									{
										// Inverse quadratic interpolation
										double r = fb / fc;
										q = fa / fc;
										p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
										q = (q - 1.0) * (r - 1.0) * (s - 1.0);
									}

								if (p > 0.0)
									{
										q = -q;
									} else
									{
										p = -p;
									}

								// Is interpolated point acceptable
								if ((2.0 * p < 3.0 * m * q - Math.abs(toler * q)) && (p < Math.abs(0.5 * e * q)))
									{
										e = d;
										d = p / q;
									} else
									{
										// Bisect
										d = m;
										e = m;
									}
							}

						// Next point.
						a = b;
						fa = fb;

						if (Math.abs(d) > toler)
							{
								b = b + d;
							} else if (b > c)
							{
								b = b - toler;
							} else
							{
								b = b + toler;
							}

						fb = applyFunc(func, b);
					}

				return b;
			}


		private static double applyFunc(Function11 func, double x) throws MathException
			{
				double y = func.apply(x);

				if (Double.isInfinite(x))
					{
						throw new MathException("Zero not found.");
					}

				if (Double.isInfinite(y))
					{
						throw new MathException("User function returned infinite value.");
					}

				if (Double.isNaN(y))
					{
						throw new MathException("User function returned NaN.");
					}

				return y;
			}

	}
