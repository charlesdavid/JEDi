package supportPlot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.chart.renderer.category.StatisticalLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import Jama.Matrix;

public class RMSIP_Plot
{
	/* ********************************* METHODS ********************************* */

	public static void create_RMSIP_Chart_ErrorBars(String out_dir, String title, String x_axis_label, String y_axis_label, Matrix data)
	{
		int ROWS = data.getRowDimension();
		DefaultStatisticalCategoryDataset dataset1 = new DefaultStatisticalCategoryDataset();
		DefaultCategoryDataset dataset2 = new DefaultCategoryDataset();
		DefaultCategoryDataset dataset3 = new DefaultCategoryDataset();
		for (int i = 0; i < ROWS; i++)
			{
				String index = String.valueOf((i + 1));
				dataset1.add(data.get(i, 1), data.get(i, 2), "Mean RMSIP", index);
				dataset2.addValue(data.get(i, 0), "RMSIP Score", index);
				dataset3.addValue(data.get(i, 3), "Z-Scores", index);
			}

		CategoryAxis domain = new CategoryAxis();
		domain.setLabel(x_axis_label);

		ValueAxis range = new NumberAxis();
		range.setLabel(y_axis_label);
		range.setRange(0, 1);

		StatisticalLineAndShapeRenderer renderer1 = new StatisticalLineAndShapeRenderer(true, true);
		LineAndShapeRenderer renderer2 = new LineAndShapeRenderer();
		LineAndShapeRenderer renderer3 = new LineAndShapeRenderer();

		renderer1.setSeriesItemLabelsVisible(0, true);
		renderer2.setSeriesItemLabelsVisible(1, true);
		renderer3.setSeriesItemLabelsVisible(2, true);

		CategoryPlot plot = new CategoryPlot(dataset1, domain, range, renderer1);
		plot.setDataset(0, dataset1);
		plot.setDataset(1, dataset2);
		plot.setDataset(2, dataset3);

		plot.setRenderer(0, renderer1);
		plot.setRenderer(1, renderer2);
		plot.setRenderer(2, renderer3);

		plot.getRendererForDataset(plot.getDataset(0)).setSeriesPaint(0, Color.blue);
		plot.getRendererForDataset(plot.getDataset(1)).setSeriesPaint(0, Color.red);
		plot.getRendererForDataset(plot.getDataset(2)).setSeriesPaint(0, Color.green);

		plot.setRangeAxis(0, range);
		plot.setRangeAxis(1, new NumberAxis("Z-SCORE"));

		plot.mapDatasetToRangeAxis(0, 0);
		plot.mapDatasetToRangeAxis(1, 0);
		plot.mapDatasetToRangeAxis(2, 1);

		LegendTitle lt = new LegendTitle(plot);

		JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, plot, false);
		chart.setAntiAlias(true);
		chart.addLegend(lt);
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

	public static void createChart2Series(String output_dir, String title, Matrix input_data)
	{
		int ROWS = input_data.getRowDimension();

		XYSeriesCollection dataset1 = new XYSeriesCollection();
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("RMSIP");
		XYSeries series2 = new XYSeries("Z-SCORE");

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
		plot.setRangeAxis(0, new NumberAxis("RMSIP"));
		plot.setRangeAxis(1, new NumberAxis("Z-SCORE"));
		plot.setDomainAxis(new NumberAxis("Subspace Dimension"));
		plot.getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getRangeAxis(0).setRange(0, 1);

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