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

import support.Atom;

public class RMSD_Plot
{
	/* ********************************* METHODS ********************************* */

	public static void create_RMSF_XY_Chart(String output_dir, String title, String x_axis_label, List<Double> input_data)
	{
		int width = (1600);
		int height = (900);

		int ROWS = input_data.size();
		if (ROWS > 1000) width = 2400;
		if (ROWS > 2000) width = 3200;

		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Atomic RMSF", false);

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i);
				series1.add((i + 1), x);
			}
		dataset.addSeries(series1);

		// construct the plot
		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());
		plot.setRangeAxis(0, new NumberAxis("RMSF (Angstrom)"));
		plot.setDomainAxis(new NumberAxis(x_axis_label));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setLowerBound(0);
		plot.getDomainAxis().setUpperBound(ROWS);
		plot.getDomainAxis().setLowerMargin(0);
		plot.getDomainAxis().setUpperMargin(0);
		// plot.getDomainAxis().setVerticalTickLabels(true);

		// generate the chart
		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void create_RMSD_XY_Chart(String output_dir, String title, String x_axis_label, List<Double> input_data)
	{
		int width = (1600);
		int height = (900);

		int ROWS = input_data.size();
		if (ROWS > 1000) width = 2400;
		if (ROWS > 2000) width = 3200;

		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("RMSD", false);

		for (int i = 0; i < ROWS; i++)
			{
				double x = input_data.get(i);
				series1.add((i + 1), x);
			}
		dataset.addSeries(series1);

		// construct the plot
		XYPlot plot = new XYPlot();
		plot.setDataset(0, dataset);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new XYLineAndShapeRenderer());
		plot.setRangeAxis(0, new NumberAxis("RMSD"));
		plot.setDomainAxis(new NumberAxis(x_axis_label));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setLowerBound(0);
		plot.getDomainAxis().setUpperBound(ROWS);
		plot.getDomainAxis().setLowerMargin(0);
		plot.getDomainAxis().setUpperMargin(0);
		// plot.getDomainAxis().setVerticalTickLabels(true);

		// generate the chart
		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.WHITE);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}

	public static void create_RMSF_Line_Chart(String output_dir, String title, String x_axis_label, List<Double> input_data, List<Atom> atoms)
	{
		int width = (1600);
		int height = (900);

		int num_atoms = atoms.size();
		if (num_atoms > 1000) width = 2400;
		if (num_atoms > 2000) width = 3200;

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		for (int i = 0; i < num_atoms; i++)
			{
				double x = input_data.get(i);
				String index = Integer.toString(atoms.get(i).atom_number);
				dataset.addValue(x, "Atomic RMSF", index);
			}
		// construct the plot
		CategoryPlot plot = new CategoryPlot();
		plot.setDataset(0, dataset);

		// customize the plot with renderers and axis
		plot.setRenderer(0, new DefaultCategoryItemRenderer());// use default fill paint for first series

		plot.setRangeAxis(0, new NumberAxis("RMSF (Angstroms)"));
		plot.setDomainAxis(new CategoryAxis(x_axis_label));

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
				ChartUtilities.saveChartAsPNG(new File(output_dir + title + ".png"), chart, width, height);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}
	}
}