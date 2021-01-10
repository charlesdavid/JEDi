package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class Plot_XY_Line_2_Series
{
	int ROWS;
	Matrix data, col1, col2;
	String out_dir, chart_title, xAxis, yAxis;

	/* ********************************* CONSTRUCTORS ********************************* */

	public Plot_XY_Line_2_Series(String output_dir, String title, Matrix input_data)
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.xAxis = "x";
		this.yAxis = "PDF(x)";
		this.data = input_data;
		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);

		saveChart(createChart(chart_title, xAxis, yAxis));
	}

	public Plot_XY_Line_2_Series(String output_dir, String title, String x_ax, String y_ax, Matrix input_data)
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.xAxis = x_ax;
		this.yAxis = y_ax;
		this.data = input_data;
		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);

		saveChart(createChart(chart_title, xAxis, yAxis));
	}

	/* ********************************* METHODS ********************************* */

	public JFreeChart createChart(String title, String x_axis_label, String y_axis_label)
	{
		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("Variance");
		XYSeries series2 = new XYSeries("Cumulative Variance");

		for (int i = 0; i < ROWS; i++)
			{
				double x = col1.get(i, 0);
				double y = col2.get(i, 0);
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
		plot.setRenderer(0, new XYSplineRenderer());// use default fill paint for first series
		XYSplineRenderer splinerenderer = new XYSplineRenderer();
		splinerenderer.setSeriesFillPaint(0, Color.BLUE);
		plot.setRenderer(1, splinerenderer);
		plot.setRangeAxis(0, new NumberAxis("Variance"));
		plot.setRangeAxis(1, new NumberAxis("Cumulative Variance"));
		plot.setDomainAxis(new NumberAxis("Mode Index"));

		// Map the data to the appropriate axis
		plot.mapDatasetToRangeAxis(0, 0);
		plot.mapDatasetToRangeAxis(1, 1);

		// generate the chart
		JFreeChart chart = new JFreeChart(chart_title, plot);
		chart.setBackgroundPaint(Color.WHITE);
		// JPanel jpanel = new ChartPanel(chart);


		// NEW PART THAT MAKES IT WORK
		// ChartPanel chartPanel = new ChartPanel(chart);
		// chartPanel.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));
		// chartPanel.setBackground(Color.white);
		// add(chartPanel);


		return chart;
	}

	public void displayChart(JFreeChart chart)
	{
		ChartFrame frame = new ChartFrame(chart_title, chart);
		frame.pack();
		frame.setVisible(true);
	}

	public void saveChart(JFreeChart chart)
	{
		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_dir + chart_title + ".png"), chart, 1600, 900);
			}
		catch (IOException ex)
			{
				ex.printStackTrace();
			}

	}
}