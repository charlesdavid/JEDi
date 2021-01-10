package supportPlot;

import java.awt.Color;
import java.io.File;

import com.orsoncharts.Chart3D;
import com.orsoncharts.Chart3DFactory;
import com.orsoncharts.Colors;
import com.orsoncharts.axis.ValueAxis3D;
import com.orsoncharts.data.xyz.XYZSeries;
import com.orsoncharts.data.xyz.XYZSeriesCollection;
import com.orsoncharts.graphics3d.Dimension3D;
import com.orsoncharts.graphics3d.ExportUtils;
import com.orsoncharts.graphics3d.ViewPoint3D;
import com.orsoncharts.legend.LegendAnchor;
import com.orsoncharts.plot.XYZPlot;
import com.orsoncharts.renderer.xyz.ScatterXYZRenderer;
import com.orsoncharts.util.Orientation;

import Jama.Matrix;


public class FES_Plot
{
	String out_dir, chart_title, chart_subtitle, xAxis, yAxis, zAxis, view_Point;
	Matrix data, col1, col2, col3;
	int ROWS;
	double maxX, minX, maxZ, minZ;
	XYZSeriesCollection<Double> FES_Data;
	XYZSeries<Double> FES_Series1, FES_Series2, FES_Series3, FES_Series4;
	Chart3D chart;

	public FES_Plot(String output_dir, String title, String subtitle, String X, String Y, String Z, Matrix input_data, String viewPoint)
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.chart_subtitle = subtitle;
		this.xAxis = X;
		this.yAxis = Y;
		this.zAxis = Z;
		this.data = input_data;
		this.view_Point = viewPoint;

		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);
		this.col3 = data.getMatrix(0, ROWS - 1, 2, 2);

		createDataset(col1, col2, col3);
		createChart();

		try
			{
				ExportUtils.writeAsPNG(chart, 1600, 1600, new File(out_dir + File.separatorChar + chart_title + ".png"));
			}
		catch (Exception e)
			{
				System.err.println("Error writing the file:" + out_dir + File.separatorChar + chart_title + ".png");
				e.printStackTrace();
			}
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void createDataset(Matrix m1, Matrix m2, Matrix m3)
	{
		FES_Data = new XYZSeriesCollection();

		FES_Series1 = new XYZSeries("First Quarter");
		FES_Series2 = new XYZSeries("Second Quarter");
		FES_Series3 = new XYZSeries("Third Quarter");
		FES_Series4 = new XYZSeries("Fourth Quarter");

		maxX = 0;
		minX = 0;
		maxZ = 0;
		minZ = 0;

		int quarter = (ROWS / 4);

		for (int i = 0; i < m1.getRowDimension(); i++)
			{
				double x = m1.get(i, 0);
				if (x > maxX) maxX = x;
				if (x < minX) minX = x;

				double z = m2.get(i, 0);
				if (z > maxX) maxZ = z;
				if (z < minZ) minZ = z;

				double y = m3.get(i, 0);

				if (i < quarter) FES_Series1.add(x, y, z);
				if (i < 2 * quarter & i >= quarter - 1) FES_Series2.add(x, y, z);
				if (i < 3 * quarter & i >= 2 * quarter - 1) FES_Series3.add(x, y, z);
				if (i < 4 * quarter & i >= 3 * quarter - 1) FES_Series4.add(x, y, z);
			}

		FES_Data.add(FES_Series1);
		FES_Data.add(FES_Series2);
		FES_Data.add(FES_Series3);
		FES_Data.add(FES_Series4);
	}


	public Chart3D createChart()
	{
		chart = Chart3DFactory.createScatterChart(chart_title, chart_subtitle, FES_Data, xAxis, yAxis, zAxis);
		chart.setLegendPosition(LegendAnchor.BOTTOM_RIGHT, Orientation.VERTICAL);
		chart.setAntiAlias(true);

		XYZPlot plot = (XYZPlot) chart.getPlot();
		plot.setDataset(FES_Data);

		ValueAxis3D xAxis = plot.getXAxis();
		xAxis.setRange(minX, maxX);

		ValueAxis3D zAxis = plot.getZAxis();
		zAxis.setRange(minZ, maxZ);

		plot.setDimensions(new Dimension3D(20, 20, 20));

		ScatterXYZRenderer renderer = (ScatterXYZRenderer) plot.getRenderer();
		renderer.setSize(0.15);
		renderer.setColors(Colors.createIntenseColors());
		renderer.setColors(Color.RED, Color.BLUE, Color.GREEN, Color.ORANGE);

		plot.setRenderer(renderer);

		if (view_Point.contentEquals("LEFT")) chart.setViewPoint(ViewPoint3D.createAboveLeftViewPoint(40));
		if (view_Point.contentEquals("RIGHT")) chart.setViewPoint(ViewPoint3D.createAboveRightViewPoint(40));

		return chart;
	}
}
