package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class COLLECTIVITY_Plot
{
	@SuppressWarnings("deprecation")
	public static void createPlot(String out_dir, String title, Matrix data)
	{
		// JFreeChart chart = ChartFactory.createXYLineChart(title, "Mode Number", "Collectivity", createDataset(data), PlotOrientation.VERTICAL, true, true, false);
		JFreeChart chart = ChartFactory.createLineChart(title, "Mode Number", "Collectivity", createDataset(data), PlotOrientation.VERTICAL, true, false, false);
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		LineAndShapeRenderer renderer = (LineAndShapeRenderer) plot.getRenderer();
		renderer.setShapesVisible(true);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		ValueAxis range = plot.getRangeAxis();
		range.setLowerBound(0);
		range.setUpperBound(1.00);
		plot.setRangeAxis(range);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	private static DefaultCategoryDataset createDataset(Matrix m)
	{
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < m.getRowDimension(); i++)
			{
				String index = Integer.toString(i + 1);
				dataset.addValue(m.get(i, 0), "collectivity", index);
			}
		return dataset;
	}

	@SuppressWarnings("unused")
	private static XYDataset createDataset2(Matrix m1)
	{
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Eigenvector Collectivity Plot");

		for (int i = 0; i < m1.getRowDimension(); i++)
			{
				double x = i + 1;
				double y = m1.get(i, 0);
				series.add(x, y);
			}
		result.addSeries(series);
		return result;
	}
}