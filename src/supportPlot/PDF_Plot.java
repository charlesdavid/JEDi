package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class PDF_Plot
{
	public static void createChart(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension(); // The input matrix is expected to have 2 columns: x and PDF(x)

		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(title);

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series.add(x, y);
			}
		result.addSeries(series);
		JFreeChart chart = ChartFactory.createXYLineChart(title, "x", "PDF(x)", result, PlotOrientation.VERTICAL, true, true, false);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void createChart2Series(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension(); // The input matrix is expected to have 3 columns: x, PDF(x), and CDF(x)

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();

		XYSeries series1 = new XYSeries("PDF(x)");
		XYSeries series2 = new XYSeries("CDF(x)");

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series1.add(x, y);
				double z = input_data.get(i, 2);
				series2.add(x, z);
			}

		series1.setDescription("PDF(x)");
		series2.setDescription("CDF(x)");

		dataset1.addSeries(series1);
		dataset2.addSeries(series2);

		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());// use default fill paint for first series
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);
		plot.setRangeAxis(0, new NumberAxis("Density"));
		plot.setRangeAxis(1, new NumberAxis("Cumulative Density"));
		plot.setDomainAxis(new NumberAxis("x"));

		plot.getRangeAxis(1).setRange(0, 1);

		// Map the data to the appropriate axis
		plot.mapDatasetToRangeAxis(0, 0);
		plot.mapDatasetToRangeAxis(1, 1);

		// generate the chart
		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}
}