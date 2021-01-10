package supportPlot;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
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

public class PC_Plot
{
	public static void createChart(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension(); // The input matrix is expected to have 2 columns: PC1 and PC2

		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(title);

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series.add(x, y);
			}
		result.addSeries(series);

		JFreeChart chart = ChartFactory.createScatterPlot(title, "PC1", "PC2", result, PlotOrientation.VERTICAL, true, true, false);
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

	@SuppressWarnings("deprecation")
	public static void createChart4Series(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension(); // The input matrix is expected to have 2 columns: PC1 and PC2
		int quarter = (ROWS / 4);

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeriesCollection dataset3 = new XYSeriesCollection();
		XYSeriesCollection dataset4 = new XYSeriesCollection();

		XYSeries series1 = new XYSeries("First Quarter");
		XYSeries series2 = new XYSeries("Second Quarter");
		XYSeries series3 = new XYSeries("Third Quarter");
		XYSeries series4 = new XYSeries("Fourth Quarter");

		double minX = Double.POSITIVE_INFINITY, maxX = 0, minY = Double.POSITIVE_INFINITY, maxY = 0;

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				if (x < minX) minX = x;
				if (x > maxX) maxX = x;

				double y = input_data.get(i, 1);
				if (y < minY) minY = y;
				if (y > maxY) maxY = y;

				if (i < quarter) series1.add(x, y);
				if (i < 2 * quarter & i >= quarter - 1) series2.add(x, y);
				if (i < 3 * quarter & i >= 2 * quarter - 1) series3.add(x, y);
				if (i < 4 * quarter & i >= 3 * quarter - 1) series4.add(x, y);
			}

		dataset1.addSeries(series1);
		dataset2.addSeries(series2);
		dataset3.addSeries(series3);
		dataset4.addSeries(series4);

		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);
		plot.setDataset(2, dataset3);
		plot.setDataset(3, dataset4);

		// customize the plot with renderers
		XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesFillPaint(0, Color.RED);
		renderer0.setLinesVisible(false);
		plot.setRenderer(0, renderer0);

		XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesFillPaint(1, Color.BLUE);
		renderer1.setLinesVisible(false);
		plot.setRenderer(1, renderer1);

		XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
		renderer2.setSeriesFillPaint(2, Color.GREEN);
		renderer2.setLinesVisible(false);
		plot.setRenderer(2, renderer2);

		Rectangle rect = new Rectangle();
		Dimension dd = new Dimension();
		dd.setSize(6, 6);
		rect.setSize(dd);


		XYLineAndShapeRenderer renderer3 = new XYLineAndShapeRenderer();
		renderer3.setSeriesShape(3, rect);
		renderer3.setSeriesFillPaint(3, Color.ORANGE);
		renderer3.setLinesVisible(false);
		plot.setRenderer(3, renderer3);

		// Assign the axes and set the ranges
		plot.setRangeAxis(0, new NumberAxis("PC2"));
		plot.setDomainAxis(new NumberAxis("PC1"));

		plot.getDomainAxis().setRange(minX, maxX);
		plot.getRangeAxis().setRange(minY, maxY);

		plot.getRendererForDataset(plot.getDataset(0)).setSeriesPaint(0, Color.RED);
		plot.getRendererForDataset(plot.getDataset(1)).setSeriesPaint(0, Color.BLUE);
		plot.getRendererForDataset(plot.getDataset(2)).setSeriesPaint(0, Color.GREEN);
		plot.getRendererForDataset(plot.getDataset(3)).setSeriesPaint(0, Color.ORANGE);

		plot.getRendererForDataset(plot.getDataset(3)).setShape(rect);

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