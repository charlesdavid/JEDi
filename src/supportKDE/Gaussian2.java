/*
 * Copyright (c) 2014, Massachusetts Institute of Technology Released under the BSD 2-Clause License http://opensource.org/licenses/BSD-2-Clause
 */

package supportKDE;

/**
 * @author Philip DeCamp
 */
public class Gaussian2 implements Function21
	{


		public static Gaussian2 fromCovariance(double exx, double eyy, double exy)
			{
				double sigX = Math.sqrt(exx);
				double sigY = Math.sqrt(eyy);

				return new Gaussian2(0f, 0f, sigX, sigY, exy / (sigX * sigY), 1f);
			}


		public static Gaussian2 fromCovariance(double meanX, double meanY, double exx, double eyy, double exy)
			{
				double sigX = Math.sqrt(exx);
				double sigY = Math.sqrt(eyy);
				double p = exy / (sigX * sigY);

				return new Gaussian2(meanX, meanY, sigX, sigY, p, 1.0f);
			}


		public static Gaussian2 fromCovariance(double meanX, double meanY, double exx, double eyy, double exy, double integral)
			{
				double sigX = Math.sqrt(exx);
				double sigY = Math.sqrt(eyy);
				double p = exy / (sigX * sigY);

				return new Gaussian2(meanX, meanY, sigX, sigY, p, integral);
			}


		public static Gaussian2 fromSigma(double sigX, double sigY, double corr)
			{
				return new Gaussian2(0f, 0f, sigX, sigY, corr, 1f);
			}


		public static Gaussian2 fromSigma(double meanX, double meanY, double sigX, double sigY, double corr)
			{
				return new Gaussian2(meanX, meanY, sigX, sigY, corr, 1.0f);
			}


		public static Gaussian2 fromSigma(double meanX, double meanY, double sigX, double sigY, double corr, double integral)
			{
				return new Gaussian2(meanX, meanY, sigX, sigY, corr, integral);
			}


		private final double mMeanX;
		private final double mMeanY;
		private final double mSigmaX;
		private final double mSigmaY;
		private final double mCorr;
		private final double mIntegral;

		private final double mCoefVolume;
		private final double mCoefXX;
		private final double mCoefYY;
		private final double mCoefXY;


		private Gaussian2(double meanX, double meanY, double sigmaX, double sigmaY, double corr, double integral)
			{
				mMeanX = meanX;
				mMeanY = meanY;
				mSigmaX = sigmaX;
				mSigmaY = sigmaY;
				mCorr = corr;
				mIntegral = integral;

				final double ompp = 1 - corr * corr;

				mCoefVolume = (integral / (2.0 * Math.PI * sigmaX * sigmaY * Math.sqrt(ompp)));
				mCoefXX = -1.0f / (2.0f * sigmaX * sigmaX * ompp);
				mCoefYY = -1.0f / (2.0f * sigmaY * sigmaY * ompp);
				mCoefXY = corr / (sigmaX * sigmaY * ompp);
			}


		public double meanX()
			{
				return mMeanX;
			}

		public double meanY()
			{
				return mMeanY;
			}

		public double exx()
			{
				return mSigmaX * mSigmaX;
			}

		public double eyy()
			{
				return mSigmaY * mSigmaY;
			}

		public double exy()
			{
				return mCorr * mSigmaX * mSigmaY;
			}

		public double sigmaX()
			{
				return mSigmaX;
			}

		public double sigmaY()
			{
				return mSigmaY;
			}

		public double correlation()
			{
				return mCorr;
			}

		public double integral()
			{
				return mIntegral;
			}

		public void contour(double val, double[] outX, double[] outY)
			{
				val = Math.log(val) - Math.log(mCoefVolume);

				for (int i = 0; i < outX.length; i++)
					{
						double ang = i * Math.PI * 2.0 / outX.length;
						double cos = Math.cos(ang);
						double sin = Math.sin(ang);

						double d = Math.sqrt(val / (mCoefXX * cos * cos + mCoefYY * sin * sin + mCoefXY * cos * sin));

						outX[i] = d * cos + mMeanX;
						outY[i] = d * sin + mMeanY;
					}
			}


		public Gaussian2 setMean(double x, double y)
			{
				return new Gaussian2(x, y, mSigmaX, mSigmaY, mCorr, mIntegral);
			}

		public Gaussian2 setCovariance(double exx, double eyy, double exy)
			{
				double sigmaX = Math.sqrt(exx);
				double sigmaY = Math.sqrt(eyy);
				double p = exy / (sigmaX * sigmaY);

				return new Gaussian2(mMeanX, mMeanY, sigmaX, sigmaY, p, mIntegral);
			}

		public Gaussian2 setSigma(double sigmaX, double sigmaY, double p)
			{
				return new Gaussian2(mMeanX, mMeanY, sigmaX, sigmaY, p, mIntegral);
			}

		public Gaussian2 setIntegral(double integral)
			{
				return new Gaussian2(mMeanX, mMeanY, mSigmaX, mSigmaY, mCorr, integral);
			}


		public Gaussian2 translate(double dx, double dy)
			{
				return new Gaussian2(mMeanX + dx, mMeanY + dy, mSigmaX, mSigmaY, mCorr, mIntegral);
			}

		public Gaussian2 scaleDomain(double scale)
			{
				return new Gaussian2(mMeanX * scale, mMeanY * scale, mSigmaX * Math.abs(scale), mSigmaY * Math.abs(scale), mCorr, mIntegral);
			}

		public Gaussian2 scaleDomain(double scaleX, double scaleY)
			{
				return new Gaussian2(mMeanX * scaleX, mMeanY * scaleY, mSigmaX * Math.abs(scaleX), mSigmaY * Math.abs(scaleY), mCorr, mIntegral);
			}

		public Gaussian2 scaleSigma(double scale)
			{
				return new Gaussian2(mMeanX, mMeanY, mSigmaX * scale, mSigmaY * scale, mCorr, mIntegral);
			}

		public Gaussian2 scaleSigma(double scaleX, double scaleY)
			{
				return new Gaussian2(mMeanX, mMeanY, mSigmaX * scaleX, mSigmaY * scaleY, mCorr, mIntegral);
			}

		public Gaussian2 scaleRange(double scale)
			{
				return new Gaussian2(mMeanX, mMeanY, mSigmaX, mSigmaY, mCorr, mIntegral * scale);
			}


		@Override
		public double apply(double x, double y)
			{
				x -= mMeanX;
				y -= mMeanY;
				return mCoefVolume * Math.exp(mCoefXX * x * x + mCoefYY * y * y + mCoefXY * x * y);
			}


		@Override
		public String toString()
			{
				return String.format("Gaussian [ux: %f, uy: %f, sx: %f, sy: %f, p: %f, vol: %f]", mMeanX, mMeanY, mSigmaX, mSigmaY, mCorr, mIntegral);
			}

		@Override
		public boolean equals(Object obj)
			{
				if (!(obj instanceof Gaussian2))
					{
						return false;
					}

				Gaussian2 g = (Gaussian2) obj;

				return this == g || mMeanX == g.mMeanX && mMeanY == g.mMeanY && mSigmaX == g.mSigmaX && mSigmaY == g.mSigmaY && mCorr == g.mCorr && mIntegral == g.mIntegral;
			}

		@Override
		public int hashCode()
			{
				return (int) Double.doubleToLongBits(mMeanX + mMeanY + mSigmaX + mSigmaY);

			}

	}
