package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.DefaultCategoryItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;
import support.Atom;

public class STATS_Plot
{
	public static void create_Variables_Stat_Line_Chart(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data, List<Atom> atoms, int col)
	{
		// Plots three series of length N from a matrix with 3N rows: The X, Y, and Z variables stats
		// 'col' should be a natural number, 1,2,3,etc...so we subtract 1 for java array numbering
		DefaultCategoryDataset dataset1 = new DefaultCategoryDataset();
		DefaultCategoryDataset dataset2 = new DefaultCategoryDataset();
		DefaultCategoryDataset dataset3 = new DefaultCategoryDataset();

		int width = (1600);
		int height = (900);

		int num_atoms = atoms.size();
		for (int i = 0; i < num_atoms; i++)
			{
				int atom_number = atoms.get(i).getAtom_number();
				String index = Integer.toString(atom_number);
				double varX = input_data.get(i, col - 1);
				double varY = input_data.get(i + num_atoms, col - 1);
				double varZ = input_data.get(i + 2 * num_atoms, col - 1);
				dataset1.addValue(varX, "X-variables", index);
				dataset2.addValue(varY, "Y-variables", index);
				dataset3.addValue(varZ, "Z-variables", index);
			}

		// construct the plot
		CategoryPlot plot = new CategoryPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);
		plot.setDataset(2, dataset3);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new DefaultCategoryItemRenderer());// use default fill paint for first series

		DefaultCategoryItemRenderer renderer = new DefaultCategoryItemRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);

		renderer = new DefaultCategoryItemRenderer();
		renderer.setSeriesFillPaint(0, Color.GREEN);
		plot.setRenderer(2, renderer);


		plot.setRangeAxis(0, new NumberAxis(y_axis_label));
		plot.setDomainAxis(new CategoryAxis(x_axis_label));
		// plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());

		CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 2.0));
		domainAxis.setTickMarksVisible(true);
		domainAxis.setLowerMargin(0);
		domainAxis.setUpperMargin(0);

		// Label Control and Clarity: "Turn off" labels based on number of points to plot
		if (num_atoms >= 50 & num_atoms < 100)
			{
				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 2 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (num_atoms >= 100 & num_atoms < 250)
			{
				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 5 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (num_atoms >= 250 & num_atoms < 500)
			{
				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 10 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (num_atoms >= 500 & num_atoms < 1000)
			{
				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 15 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (num_atoms >= 1000 & num_atoms <= 2000)
			{
				width = 2400;

				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 20 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		if (num_atoms > 2000)
			{
				width = 3200;

				for (int i = 0; i < num_atoms; i++)
					{
						if (!(i % 40 == 0))
							{
								int atom_number = atoms.get(i).getAtom_number();
								String index = Integer.toString(atom_number);
								domainAxis.setTickLabelPaint(index, Color.white);
							}
					}
			}

		plot.configureDomainAxes();

		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void create_Variables_Stat_XY_Chart(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data, List<Atom> atoms, int col)
	{
		// Plots three series of length N from a matrix with 3N rows: The X, Y, and Z variables stats
		// The value of 'col' should be a natural number, 1,2,3,etc...so we subtract 1 for java array numbering
		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeriesCollection dataset3 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("X-variables", false);
		XYSeries series2 = new XYSeries("Y-variables", false);
		XYSeries series3 = new XYSeries("Z-variables", false);

		int width = (1600);
		int height = (900);

		int num_atoms = atoms.size();
		if (num_atoms > 1000) width = 2400;
		if (num_atoms > 2000) width = 3200;

		for (int i = 0; i < num_atoms; i++)
			{
				int atom_number = (i + 1);
				double varX = input_data.get(i, col - 1);
				double varY = input_data.get(i + num_atoms, col - 1);
				double varZ = input_data.get(i + 2 * num_atoms, col - 1);
				series1.add(atom_number, varX);
				series2.add(atom_number, varY);
				series3.add(atom_number, varZ);
			}
		dataset1.addSeries(series1);
		dataset2.addSeries(series2);
		dataset3.addSeries(series3);

		// construct the plot
		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);
		plot.setDataset(2, dataset3);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());// use default fill paint for first series

		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);

		renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.GREEN);
		plot.setRenderer(2, renderer);
		plot.setRangeAxis(0, new NumberAxis(y_axis_label));

		plot.setDomainAxis(new NumberAxis(x_axis_label));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setVerticalTickLabels(true);
		plot.getDomainAxis().setLowerMargin(0);
		plot.getDomainAxis().setUpperMargin(0);

		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void create_Variables_Stat_XY_Chart(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data, int col)
	{
		// Plots three series of length N from a matrix with 3N rows: The X, Y, and Z variables stats
		// The value of 'col' should be a natural number, 1,2,3,etc...so we subtract 1 for java array numbering
		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeriesCollection dataset3 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("X-variables", false);
		XYSeries series2 = new XYSeries("Y-variables", false);
		XYSeries series3 = new XYSeries("Z-variables", false);

		int width = (1600);
		int height = (900);

		int num_atoms = input_data.getRowDimension() / 3;
		if (num_atoms > 1000) width = 2400;
		if (num_atoms > 2000) width = 3200;

		for (int i = 0; i < num_atoms; i++)
			{
				int atom_number = (i + 1);
				double varX = input_data.get(i, col - 1);
				double varY = input_data.get(i + num_atoms, col - 1);
				double varZ = input_data.get(i + 2 * num_atoms, col - 1);
				series1.add(atom_number, varX);
				series2.add(atom_number, varY);
				series3.add(atom_number, varZ);
			}
		dataset1.addSeries(series1);
		dataset2.addSeries(series2);
		dataset3.addSeries(series3);

		// construct the plot
		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);
		plot.setDataset(2, dataset3);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());// use default fill paint for first series

		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);

		renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesFillPaint(0, Color.GREEN);
		plot.setRenderer(2, renderer);
		plot.setRangeAxis(0, new NumberAxis(y_axis_label));

		plot.setDomainAxis(new NumberAxis(x_axis_label));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setVerticalTickLabels(true);
		plot.getDomainAxis().setLowerMargin(0);
		plot.getDomainAxis().setUpperMargin(0);

		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);
		chart.getPlot().setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}
}