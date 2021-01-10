package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;
import support.Atom;

public class Modes_XYPlot
{
	/* ********************************* METHODS ********************************* */

	public static void create_XY_Chart_2Series(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension();

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Mode 1");
		XYSeries series2 = new XYSeries("Mode 2");

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
		plot.setRangeAxis(0, new NumberAxis("Variance Mode 1"));
		plot.setRangeAxis(1, new NumberAxis("Variance Mode 2"));
		plot.setDomainAxis(new NumberAxis("Atom Index"));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());

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
}