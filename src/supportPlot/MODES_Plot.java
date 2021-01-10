package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.DefaultCategoryItemRenderer;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;
import support.Atom;

public class MODES_Plot
{
	public static void create_Line_Chart_2Series(String out_dir, String title, Matrix input_data, List<Atom> atoms)
	{
		DefaultCategoryDataset dataset1 = new DefaultCategoryDataset();
		DefaultCategoryDataset dataset2 = new DefaultCategoryDataset();

		int num_atoms = atoms.size();
		for (int i = 0; i < num_atoms; i++)
			{
				int atom_number = atoms.get(i).getAtom_number();
				String index = Integer.toString(atom_number);
				double X = input_data.get(i, 0);
				double Y = input_data.get(i, 1);
				dataset1.addValue(X, "Mode 1", index);
				dataset2.addValue(Y, "Mode 2", index);
			}

		// construct the plot
		CategoryPlot plot = new CategoryPlot();
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new DefaultCategoryItemRenderer()); // use default fill paint for first series

		DefaultCategoryItemRenderer renderer = new DefaultCategoryItemRenderer();
		renderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, renderer);

		plot.setRangeAxis(0, new NumberAxis("Variance (Angstrom^2)"));
		plot.setDomainAxis(new CategoryAxis("Atom Number"));

		CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 2.0));
		domainAxis.setTickMarksVisible(true);
		// domainAxis.setLabelFont(new Font("Calibri", Font.PLAIN, 2));
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
				ChartUtilities.saveChartAsPNG(new File(out_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	@SuppressWarnings("deprecation")
	public static void createLineChart2Series(String output_dir, String title, Matrix input_data, List<Atom> atoms)
	{
		int ROWS = input_data.getRowDimension();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				String index = Integer.toString(atoms.get(i).atom_number);
				dataset.addValue(x, "Mode 1", index);
				dataset.addValue(y, "Mode 2", index);
			}

		// generate the chart
		JFreeChart chart = ChartFactory.createLineChart(title, "Atom Index", "Variance (Angstrom^2)", dataset, PlotOrientation.VERTICAL, true, false, false);
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
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	@SuppressWarnings("deprecation")
	public static void createLineChart1SeriesDist(String output_dir, String title, String x_axis_label, String y_axis_label, Matrix input_data, List<Integer> atoms1,
			List<Integer> atoms2, int col)
	{
		int ROWS = input_data.getRowDimension();
		int COL = col - 1;

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		for (int i = 0; i < ROWS; i++)
			{
				double val = input_data.get(i, COL);
				String index = Integer.toString(atoms1.get(i)) + "---" + Integer.toString(atoms2.get(i));
				dataset.addValue(val, "Mode 1", index);
			}

		// generate the chart
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
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	@SuppressWarnings("deprecation")
	public static void createLineChart2SeriesDist(String output_dir, String title, Matrix input_data, List<Atom> atoms1, List<Atom> atoms2)
	{
		int ROWS = input_data.getRowDimension();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				String index = Integer.toString(atoms1.get(i).atom_number) + "---" + Integer.toString(atoms2.get(i).atom_number);
				dataset.addValue(x, "Mode 1", index);
				dataset.addValue(y, "Mode 2", index);
			}

		// generate the chart
		JFreeChart chart = ChartFactory.createLineChart(title, "Atom Index", "Variance (Angstrom^2)", dataset, PlotOrientation.VERTICAL, true, false, false);
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
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void create_XY_Chart_2Series(String output_dir, String title, Matrix input_data, List<Atom> atoms)
	{
		int ROWS = input_data.getRowDimension();

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Mode 1");
		XYSeries series2 = new XYSeries("Mode 2");

		for (int i = 0; i < ROWS; i++)
			{
				int atom_number = atoms.get(i).getAtom_number();
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series1.add(atom_number, x);
				series2.add(atom_number, y);
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
		plot.setRangeAxis(0, new NumberAxis("Variance Mode 1"));
		plot.setRangeAxis(1, new NumberAxis("Variance Mode 2"));
		plot.setDomainAxis(new NumberAxis("Atom Number"));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setRange(atoms.get(0).getAtom_number(), atoms.get(ROWS - 1).getAtom_number());

		// Map the data to the appropriate axis
		plot.mapDatasetToRangeAxis(0, 0);
		plot.mapDatasetToRangeAxis(1, 1);

		// generate the chart
		JFreeChart chart = new JFreeChart(title, plot);
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

	public static void create_XY_Chart_2Series(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension();

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Mode 1");
		XYSeries series2 = new XYSeries("Mode 2");

		for (int i = 0; i < ROWS; i++)
			{
				int atom_number = i + 1;
				double x = input_data.get(i, 0);
				double y = input_data.get(i, 1);
				series1.add(atom_number, x);
				series2.add(atom_number, y);
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
		plot.setRangeAxis(0, new NumberAxis("Variance Mode 1"));
		plot.setRangeAxis(1, new NumberAxis("Variance Mode 2"));
		plot.setDomainAxis(new NumberAxis("Atom Index"));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setRange(0, ROWS);

		// Map the data to the appropriate axis
		plot.mapDatasetToRangeAxis(0, 0);
		plot.mapDatasetToRangeAxis(1, 1);

		// generate the chart
		JFreeChart chart = new JFreeChart(title, plot);
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
}