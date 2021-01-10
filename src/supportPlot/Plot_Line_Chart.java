package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;

import Jama.Matrix;

public class Plot_Line_Chart
{
	@SuppressWarnings("deprecation")
	public static void createChart_One_Series(String out_dir, String title, String x_axis_label, String y_axis_label, List<Double> input_data)
	{
		// Plots a single series from a list
		int ROWS = input_data.size();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				dataset.addValue(input_data.get(i), y_axis_label, x_axis_label);
			}
		JFreeChart chart = ChartFactory.createLineChart(title, x_axis_label, y_axis_label, dataset, PlotOrientation.VERTICAL, true, true, false);

		// customize the plot with renderers and axis
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		LineAndShapeRenderer renderer = (LineAndShapeRenderer) plot.getRenderer();
		renderer.setShapesVisible(true);
		CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 2.0));

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	@SuppressWarnings("deprecation")
	public static void createChart_One_Series(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data)
	{
		// Plots a single series from a matrix
		int ROWS = input_data.getRowDimension();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				// dataset.addValue(input_data.get(i, 0), y_axis_label, "Mode " + (i + 1));
				dataset.addValue(input_data.get(i, 0), y_axis_label, "" + (i + 1));
			}
		JFreeChart chart = ChartFactory.createLineChart(title, x_axis_label, y_axis_label, dataset, PlotOrientation.VERTICAL, true, false, false);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		// customize the plot with renderers and axis
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		LineAndShapeRenderer renderer = (LineAndShapeRenderer) plot.getRenderer();
		renderer.setShapesVisible(true);
		CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 2.0));

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	@SuppressWarnings("deprecation")
	public static void createChart_One_Series(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data, int column)
	{
		// Plots a single series from a matrix using specified column
		int ROWS = input_data.getRowDimension();
		int COL = column - 1;

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				dataset.addValue(input_data.get(i, COL), y_axis_label, "" + (i + 1));
			}
		JFreeChart chart = ChartFactory.createLineChart(title, x_axis_label, y_axis_label, dataset, PlotOrientation.VERTICAL, true, false, false);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		CategoryPlot plot = (CategoryPlot) chart.getPlot();

		ValueAxis range = plot.getRangeAxis();
		range.setLowerBound(0);
		range.setUpperBound(1.00);
		plot.setRangeAxis(range);

		// customize the plot with renderers and axis
		LineAndShapeRenderer renderer = (LineAndShapeRenderer) plot.getRenderer();
		renderer.setShapesVisible(true);
		CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 2.0));

		// Label Control and Clarity:
		if (ROWS >= 50 & ROWS < 100)
			{
				for (int i = 0; i < ROWS; i++)
					{
						if ((i + 1) % 2 == 0)
							{
								String index = "" + (i + 1);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (ROWS >= 100 & ROWS < 200)
			{
				for (int i = 0; i < ROWS; i++)
					{
						if (!((i) % 3 == 0))
							{
								String index = "" + (i + 1);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (ROWS >= 300)
			{
				for (int i = 0; i < ROWS; i++)
					{
						if (!((i) % 4 == 0))
							{
								String index = "" + (i + 1);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void createChart_Two_Series(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data)
	{
		// Plots two series from a matrix
		int ROWS = input_data.getRowDimension();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				dataset.addValue(input_data.get(i, 0), "MSD", "Residue " + (i + 1));
				dataset.addValue(input_data.get(i, 1), "MSD", "Residue " + (i + 1));
			}
		JFreeChart chart = ChartFactory.createLineChart(title, x_axis_label, y_axis_label, dataset, PlotOrientation.VERTICAL, true, true, false);
		chart.setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void createChart_Three_Series(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data)
	{
		// Plots three series from a matrix
		int ROWS = input_data.getRowDimension();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				dataset.addValue(input_data.get(i, 0), y_axis_label, "Mode " + (i + 1));
				dataset.addValue(input_data.get(i, 1), "MSD", "Residue " + (i + 1));
				dataset.addValue(input_data.get(i, 3), "MSD", "Residue " + (i + 1));
			}
		JFreeChart chart = ChartFactory.createLineChart(title, x_axis_label, y_axis_label, dataset, PlotOrientation.VERTICAL, true, true, false);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void displayChart(JFreeChart chart)
	{
		String title = chart.getTitle().toString();
		ChartFrame frame = new ChartFrame(title, chart);
		frame.pack();
		frame.setVisible(true);
	}

	public static void saveChart(JFreeChart chart, String path)
	{
		try
			{
				String title = chart.getTitle().toString();
				ChartUtilities.saveChartAsPNG(new File(path + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}
}