package supportPlot;

import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class Plot_XY_Scatter
{
	int ROWS;
	Matrix data, col1, col2;
	String out_dir, chart_title, xAxis, yAxis;

	/* ********************************* CONSTRUCTORS ********************************* */

	public Plot_XY_Scatter(String output_dir, String title, Matrix input_data)
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.xAxis = "PC1";
		this.yAxis = "PC2";
		this.data = input_data;
		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);

		saveChart(createChart(chart_title, xAxis, yAxis));
	}

	public Plot_XY_Scatter(String output_dir, String title, String x_ax, String y_ax, Matrix input_data)
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

	public XYDataset createDataset(Matrix m1, Matrix m2)
	{
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("PC-Plot");

		for (int i = 0; i < m1.getRowDimension(); i++)
			{
				double x = m1.get(i, 0);
				double y = m2.get(i, 0);
				series.add(x, y);
			}
		result.addSeries(series);
		return result;
	}

	public JFreeChart createChart(String title, String x_axis_label, String y_axis_label)
	{
		JFreeChart chart = ChartFactory.createScatterPlot(title, x_axis_label, y_axis_label, createDataset(col1, col2), PlotOrientation.VERTICAL, true, true, false);

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