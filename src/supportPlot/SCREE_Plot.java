package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class SCREE_Plot
{
	/* ********************************* METHODS ********************************* */

	public static void createChart(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension(); // The input matrix is expected to have 2 columns

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Variance");
		XYSeries series2 = new XYSeries("Cumulative Variance");

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series1.add((i + 1), x);
				series2.add((i + 1), y);
			}
		dataset1.addSeries(series1);
		dataset2.addSeries(series2);

		// construct the plot
		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());// use default fill paint for first series
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);
		plot.setRangeAxis(0, new NumberAxis("Variance"));
		plot.setRangeAxis(1, new NumberAxis("Cumulative Variance"));
		plot.setDomainAxis(new NumberAxis("Mode Number"));

		plot.getDomainAxis().setRange(1, ROWS);
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
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