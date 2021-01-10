/*
 * Copyright (c) 2014, Massachusetts Institute of Technology Released under the BSD 2-Clause License http://opensource.org/licenses/BSD-2-Clause
 */
package supportKDE;

/**
 * @author Philip DeCamp
 */
public class BilinearSampler implements Function21
	{

		private final double[] mTable;
		private final int mCols;
		private final int mRows;

		private final double mMinX;
		private final double mMinY;
		private final double mMaxX;
		private final double mMaxY;

		private final double mScaleX;
		private final double mScaleY;

		private double mOobValue = 0.0;


		public BilinearSampler(double[] tableRef, int cols, int rows, double minX, double minY, double maxX, double maxY)
			{
				mTable = tableRef;
				mCols = cols;
				mRows = rows;

				mMinX = minX;
				mMinY = minY;
				mMaxX = maxX;
				mMaxY = maxY;

				mScaleX = (cols - 1) / (maxX - minX);
				mScaleY = (rows - 1) / (maxY - minY);
			}


		/**
		 * Computes value at a given point using a bilinear interpolation within the table.
		 */
		@Override
		public double apply(double x, double y)
			{
				double px = (x - mMinX) * mScaleX;
				double py = (y - mMinY) * mScaleY;

				int xx = (int) px;
				int yy = (int) py;

				if (xx < 0 || yy < 0 || xx >= mCols - 1 || yy >= mRows - 1)
					{
						return mOobValue;
					}

				final double u = px - xx;
				final double v = py - yy;
				final int i = xx + yy * mCols;

				double ret;

				ret = mTable[i] * (1.0 - u) * (1.0 - v);
				ret += mTable[i + 1] * u * (1.0 - v);
				ret += mTable[i + mCols] * (1.0 - u) * v;
				ret += mTable[i + 1 + mCols] * u * v;

				return ret;
			}

		/**
		 * @return width of the cells used in the sample table.
		 */
		public double cellWidth()
			{
				return (mMaxX - mMinX) / mCols;
			}

		/**
		 * @return height of the cells used in the sample table.
		 */
		public double cellHeight()
			{
				return (mMaxY - mMinY) / mRows;
			}

		/**
		 * @return number of cells in one column of the sample table
		 */
		public int tableCols()
			{
				return mCols;
			}

		/**
		 * @return number of cells in one row of the sample table
		 */
		public int tableRows()
			{
				return mRows;
			}

		/**
		 * @return bounds of the sample table [minX, minY, maxX, maxY]
		 */
		public double[] bounds()
			{
				return new double[] { mMinX, mMinY, mMaxX, mMaxY };
			}

		/**
		 * @return min x value represented by the sample table
		 */
		public double minX()
			{
				return mMinX;
			}

		/**
		 * @return min y value represented by the sample table
		 */
		public double minY()
			{
				return mMinY;
			}

		/**
		 * @return max x value represented by the sample table
		 */
		public double maxX()
			{
				return mMaxX;
			}

		/**
		 * @return max y value represented by the sample table
		 */
		public double maxY()
			{
				return mMaxY;
			}

		/**
		 * @return value that <code>apply(x,y)</code> returns if [x,y] is out of bounds. Default is 0.
		 */
		public double oobValue()
			{
				return mOobValue;
			}

		/**
		 * @param value value that <code>apply(x,y)</code> should return if [x,y] is out of bounds. Default is 0.
		 */
		public void oobValue(double value)
			{
				mOobValue = value;
			}


		public double[] tableRef()
			{
				return mTable;
			}

		/**
		 * Normalizes the volume of the table to specified value.
		 */
		public void normalizeVolume(double val)
			{
				double ww = (mMaxX - mMinX) / (mCols - 1);
				double hh = (mMaxY - mMinY) / (mRows - 1);
				double cellScale = ww * hh * 0.25;
				double sum = 0.0;

				for (int y = 0; y < mRows - 1; y++)
					{
						for (int x = 0; x < mCols - 1; x++)
							{
								final int ii = x + y * mCols;

								double v00 = mTable[ii];
								double v10 = mTable[ii + 1];
								double v01 = mTable[ii + mCols];
								double v11 = mTable[ii + 1 + mCols];

								// Add four post values, divide by four, and multiple by area of cell (ww * hh).
								sum += (v00 + v10 + v01 + v11) * cellScale;
							}
					}

				double scale = val / sum;

				for (int i = 0; i < mRows * mCols; i++)
					{
						mTable[i] *= scale;
					}
			}


		public double colToX(int col)
			{
				return mMinX + col * (mMaxX - mMinX) / (mCols - 1.0);
			}


		public double rowToY(int row)
			{
				return mMinY + row * (mMaxY - mMinY) / (mRows - 1.0);
			}

	}
