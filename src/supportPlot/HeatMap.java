package supportPlot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.util.Arrays;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.RectangleEdge;

import Jama.Matrix;
import support.Descriptive_Stats;

public class HeatMap
{
	final String title;
	final String out_directory;
	final String x_axis_title;
	final String y_axis_title;
	final int number_of_residues, tickUnit;
	final double MIN, MAX, RANGE, MEAN, STD_DEV;
	final Matrix coupling_scores;

	public HeatMap(Matrix data, String out_dir, String chart_title, String xAxis, String yAxis)
	{
		this.coupling_scores = data;
		this.out_directory = out_dir;
		this.title = chart_title;
		this.x_axis_title = xAxis;
		this.y_axis_title = yAxis;
		this.number_of_residues = coupling_scores.getColumnDimension();
		this.tickUnit = (number_of_residues / 50) + 1;

		double[] SCORES = new double[number_of_residues];
		for (int i = 0; i < number_of_residues; i++)
			{
				for (int j = 1; j < number_of_residues; j++)
					{
						double score = coupling_scores.get(i, j);
						SCORES[i] = score;
					}
			}
		Descriptive_Stats ds = new Descriptive_Stats();
		double mean = ds.get_mean(SCORES);
		STD_DEV = ds.get_standard_deviation(SCORES, mean);

		Arrays.sort(SCORES);
		double min = SCORES[0];
		double max = SCORES[SCORES.length - 1];
		double val = Math.max(Math.abs(min), Math.abs(max));

		if (min == max) val = 1.00;

		/* For visualizing COV and CORR matrices, we use a symmetric range centered on zero */
		MIN = -val;
		MAX = val;
		RANGE = (MAX - MIN);
		MEAN = 0;
	}

	public JFreeChart createCoupingScorePlot()
	{
		XYZDataset dataset = createDataset();

		Font tick_label_font = new Font("Dialog", Font.PLAIN, 25);
		Font axis_label_font = new Font("Dialog", Font.PLAIN, 36);
		Font title_font = new Font("Dialog", Font.PLAIN, 72);

		NumberAxis xAxis = new NumberAxis(x_axis_title);
		xAxis.setTickUnit(new NumberTickUnit(tickUnit));
		xAxis.setTickMarksVisible(true);
		xAxis.setRange(0, number_of_residues + 1);
		xAxis.setLowerMargin(.25);
		xAxis.setUpperMargin(.25);
		xAxis.setLabelFont(axis_label_font);
		xAxis.setTickLabelFont(tick_label_font);
		xAxis.setVerticalTickLabels(true);

		NumberAxis yAxis = new NumberAxis(y_axis_title);
		yAxis.setTickUnit(new NumberTickUnit(tickUnit));
		yAxis.setTickMarksVisible(true);
		yAxis.setRange(0, number_of_residues + 1);
		yAxis.setLowerMargin(.25);
		yAxis.setUpperMargin(.25);
		yAxis.setLabelFont(axis_label_font);
		yAxis.setTickLabelFont(tick_label_font);

		// LookupPaintScale paintScale = get_Coupling_Scores_Paint_Scale_Z();


		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(0, 100, Color.black);
		paintScale.add(0.0, Color.black);
		// paintScale.add(25, Color.darkGray);
		// paintScale.add(50, Color.darkGray);
		// paintScale.add(75, Color.lightGray);
		// paintScale.add(85, Color.white);
		paintScale.add(95, Color.blue);
		// paintScale.add(95, Color.green);
		// paintScale.add(95, Color.green);
		// paintScale.add(97.5, Color.orange);
		paintScale.add(99, Color.red);



		PaintScaleLegend psl = new PaintScaleLegend(paintScale, new NumberAxis());
		psl.setPosition(RectangleEdge.RIGHT);
		psl.setAxisLocation(AxisLocation.TOP_OR_RIGHT);
		psl.setMargin(100.0, 40.0, 160.0, 50.0);
		psl.setStripOutlineStroke(new BasicStroke(10.0f));
		psl.setWidth(200);
		psl.setStripWidth(100);


		// a renderer and a plot
		XYBlockRenderer renderer = new XYBlockRenderer();
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
		renderer.setPaintScale(paintScale);
		plot.setDomainAxis(xAxis);
		plot.setRangeAxis(yAxis);
		plot.setDomainGridlinesVisible(true);
		plot.setRangeGridlinesVisible(true);
		plot.setDomainGridlinePaint(Color.GRAY);
		plot.setRangeGridlinePaint(Color.GRAY);

		plot.setDomainGridlineStroke(new BasicStroke(3.0f));
		plot.setRangeGridlineStroke(new BasicStroke(3.0f));

		// PaintScaleLegend psl = get_Paint_Scale_Legend(paintScale);

		JFreeChart chart = new JFreeChart(title, title_font, plot, false);
		chart.addSubtitle(psl);
		chart.setAntiAlias(true);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_directory + title + ".png"), chart, 3200, 3200);
			}
		catch (Exception e)
			{
				System.err.println("Sorry, could not find the specified file." + out_directory + File.separatorChar + title + ".png");
				e.printStackTrace();
			}
		return chart;
	}

	public JFreeChart createMatrixHeatmap()
	{
		XYZDataset dataset = createDataset();

		Font tick_label_font = new Font("Dialog", Font.PLAIN, 25);
		Font axis_label_font = new Font("Dialog", Font.PLAIN, 36);
		Font title_font = new Font("Dialog", Font.PLAIN, 72);

		NumberAxis xAxis = new NumberAxis(x_axis_title);
		xAxis.setTickUnit(new NumberTickUnit(tickUnit));
		xAxis.setTickMarksVisible(true);
		xAxis.setRange(0, number_of_residues + 1);
		xAxis.setLowerMargin(.25);
		xAxis.setUpperMargin(.25);
		xAxis.setLabelFont(axis_label_font);
		xAxis.setTickLabelFont(tick_label_font);
		xAxis.setVerticalTickLabels(true);

		NumberAxis yAxis = new NumberAxis(y_axis_title);
		yAxis.setTickUnit(new NumberTickUnit(tickUnit));
		yAxis.setTickMarksVisible(true);
		yAxis.setRange(0, number_of_residues + 1);
		yAxis.setLowerMargin(.25);
		yAxis.setUpperMargin(.25);
		yAxis.setLabelFont(axis_label_font);
		yAxis.setTickLabelFont(tick_label_font);

		LookupPaintScale paintScale = get_HeatMap_Paint_Scale_Z4();
		// LookupPaintScale paintScale = get_HeatMap_Paint_Scale2();

		// a renderer and a plot
		XYBlockRenderer renderer = new XYBlockRenderer();
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
		renderer.setPaintScale(paintScale);
		plot.setDomainAxis(xAxis);
		plot.setRangeAxis(yAxis);
		plot.setDomainGridlinesVisible(true);
		plot.setRangeGridlinesVisible(true);
		plot.setDomainGridlinePaint(Color.GRAY);
		plot.setRangeGridlinePaint(Color.GRAY);

		plot.setDomainGridlineStroke(new BasicStroke(3.0f));
		plot.setRangeGridlineStroke(new BasicStroke(3.0f));

		PaintScaleLegend psl = get_Paint_Scale_Legend(paintScale);

		JFreeChart chart = new JFreeChart(title, title_font, plot, false);
		chart.addSubtitle(psl);
		chart.setAntiAlias(true);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_directory + title + ".png"), chart, 3200, 3200);
			}
		catch (Exception e)
			{
				System.err.println("Sorry, could not find the specified file." + out_directory + File.separatorChar + title + ".png");
				e.printStackTrace();
			}
		return chart;
	}

	public JFreeChart createAdjacencyHeatmap()
	{
		XYZDataset dataset = createDataset();

		Font font = new Font("Dialog", Font.PLAIN, 25);
		Font font3 = new Font("Dialog", Font.PLAIN, 36);

		NumberAxis xAxis = new NumberAxis(x_axis_title);
		xAxis.setTickUnit(new NumberTickUnit(tickUnit));
		xAxis.setTickMarksVisible(true);
		xAxis.setRange(0, number_of_residues + 1);
		xAxis.setLowerMargin(.25);
		xAxis.setUpperMargin(.25);
		xAxis.setLabelFont(font3);
		xAxis.setTickLabelFont(font);
		xAxis.setVerticalTickLabels(true);

		NumberAxis yAxis = new NumberAxis(y_axis_title);
		yAxis.setTickUnit(new NumberTickUnit(tickUnit));
		yAxis.setTickMarksVisible(true);
		yAxis.setRange(0, number_of_residues + 1);
		yAxis.setLowerMargin(.25);
		yAxis.setUpperMargin(.25);
		yAxis.setLabelFont(font3);
		yAxis.setTickLabelFont(font);

		LookupPaintScale paintScale = get_HeatMap_Paint_Scale_Adjacency();

		// a renderer and a plot
		XYBlockRenderer renderer = new XYBlockRenderer();
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
		renderer.setPaintScale(paintScale);
		plot.setDomainAxis(xAxis);
		plot.setRangeAxis(yAxis);
		plot.setDomainGridlinesVisible(true);
		plot.setRangeGridlinesVisible(true);
		plot.setDomainGridlinePaint(Color.GRAY);
		plot.setRangeGridlinePaint(Color.GRAY);

		plot.setDomainGridlineStroke(new BasicStroke(3.0f));
		plot.setRangeGridlineStroke(new BasicStroke(3.0f));

		PaintScaleLegend psl = get_Paint_Scale_Legend(paintScale);

		JFreeChart chart = new JFreeChart(null, null, plot, false);
		chart.addSubtitle(psl);
		chart.setAntiAlias(true);

		try
			{
				ChartUtilities.saveChartAsPNG(new File(out_directory + title + ".png"), chart, 3200, 3200);
			}
		catch (Exception e)
			{
				System.err.println("Sorry, could not find the specified file." + out_directory + File.separatorChar + title + ".png");
				e.printStackTrace();
			}
		return chart;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_HeatMap_Paint_Scale()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN * 1.00, Color.CYAN);
		paintScale.add(MIN * 0.75, Color.CYAN);
		paintScale.add(MIN * 0.50, Color.blue);
		paintScale.add(MIN * 0.25, Color.green);
		paintScale.add(MIN * 0.10, Color.white);
		paintScale.add(MIN * 0.01, Color.lightGray);
		paintScale.add(MIN * 0.001, Color.black);
		paintScale.add(MAX * 0.001, Color.black);
		paintScale.add(MAX * 0.01, Color.lightGray);
		paintScale.add(MAX * 0.10, Color.white);
		paintScale.add(MAX * 0.25, Color.yellow);
		paintScale.add(MAX * 0.50, Color.orange);
		paintScale.add(MAX * 0.75, Color.red);
		paintScale.add(MAX * 1.00, Color.red);

		return paintScale;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_HeatMap_Paint_Scale2()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.CYAN);
		paintScale.add(MIN + (1 * (RANGE / 10)), Color.CYAN);
		paintScale.add(MIN + (2 * (RANGE / 10)), Color.blue);
		paintScale.add(MIN + (3 * (RANGE / 10)), Color.green);
		paintScale.add(MIN + (4 * (RANGE / 10)), Color.lightGray);
		paintScale.add(MIN + (4.95 * (RANGE / 10)), Color.black);
		paintScale.add(MIN + (5.05 * (RANGE / 10)), Color.black);
		paintScale.add(MIN + (6 * (RANGE / 10)), Color.white);
		paintScale.add(MIN + (7 * (RANGE / 10)), Color.yellow);
		paintScale.add(MIN + (8 * (RANGE / 10)), Color.orange);
		paintScale.add(MIN + (9 * (RANGE / 10)), Color.red);
		paintScale.add(0.95 * MAX, Color.magenta);
		paintScale.add(MAX, Color.magenta);

		return paintScale;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_HeatMap_Paint_Scale_Z()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.CYAN);
		paintScale.add(MEAN - 3.00 * STD_DEV, Color.CYAN);
		paintScale.add(MEAN - 2.00 * STD_DEV, Color.blue);
		paintScale.add(MEAN - 1.00 * STD_DEV, Color.green);
		paintScale.add(MEAN - 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN - 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN - 0.10 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.10 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN + 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN + 1.00 * STD_DEV, Color.yellow);
		paintScale.add(MEAN + 2.00 * STD_DEV, Color.orange);
		paintScale.add(MEAN + 3.00 * STD_DEV, Color.red);
		paintScale.add(MAX, Color.red);

		return paintScale;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_HeatMap_Paint_Scale_Z2()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.CYAN);
		paintScale.add(MEAN - 2.00 * STD_DEV, Color.CYAN);
		paintScale.add(MEAN - 1.50 * STD_DEV, Color.blue);
		paintScale.add(MEAN - 1.00 * STD_DEV, Color.green);
		paintScale.add(MEAN - 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN - 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN - 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN - 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN + 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN + 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN + 1.00 * STD_DEV, Color.yellow);
		paintScale.add(MEAN + 1.50 * STD_DEV, Color.orange);
		paintScale.add(MEAN + 2.00 * STD_DEV, Color.red);
		paintScale.add(MAX, Color.red);

		return paintScale;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_HeatMap_Paint_Scale_Z3()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.cyan);
		paintScale.add(MEAN - 3.00 * STD_DEV, Color.cyan);
		paintScale.add(MEAN - 2.00 * STD_DEV, Color.blue);
		paintScale.add(MEAN - 1.00 * STD_DEV, Color.green);
		paintScale.add(MEAN - 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN - 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN - 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN - 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN, Color.black);
		paintScale.add(MEAN + 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN + 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN + 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN + 1.00 * STD_DEV, Color.orange);
		paintScale.add(MEAN + 2.00 * STD_DEV, Color.red);
		paintScale.add(MEAN + 3.00 * STD_DEV, Color.magenta);
		paintScale.add(MAX, Color.magenta);

		return paintScale;
	}

	private LookupPaintScale get_HeatMap_Paint_Scale_Z4()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.blue);
		paintScale.add(MEAN - 2.00 * STD_DEV, Color.blue);
		paintScale.add(MEAN - 1.00 * STD_DEV, Color.green);
		paintScale.add(MEAN - 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN - 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN - 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN - 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.01 * STD_DEV, Color.black);
		paintScale.add(MEAN + 0.10 * STD_DEV, Color.darkGray);
		paintScale.add(MEAN + 0.25 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN + 0.50 * STD_DEV, Color.white);
		paintScale.add(MEAN + 1.00 * STD_DEV, Color.orange);
		paintScale.add(MEAN + 2.00 * STD_DEV, Color.red);
		paintScale.add(MAX, Color.red);

		return paintScale;
	}

	private LookupPaintScale get_HeatMap_Paint_Scale_Adjacency()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(-1, 1, Color.white);
		paintScale.add(-1, Color.blue);
		paintScale.add(-0.5, Color.blue);
		paintScale.add(-0.49, Color.white);
		paintScale.add(0.49, Color.white);
		paintScale.add(0.5, Color.red);
		paintScale.add(1, Color.red);

		return paintScale;
	}

	@SuppressWarnings("unused")
	private LookupPaintScale get_Coupling_Scores_Paint_Scale_Z()
	{
		// create a paint-scale and a legend showing it
		LookupPaintScale paintScale = new LookupPaintScale(MIN, MAX, Color.black);
		paintScale.add(MIN, Color.black);
		paintScale.add(MEAN, Color.darkGray);
		paintScale.add(MEAN + 1.00 * STD_DEV, Color.lightGray);
		paintScale.add(MEAN + 2.00 * STD_DEV, Color.white);
		paintScale.add(MEAN + 3.00 * STD_DEV, Color.orange);
		paintScale.add(MEAN + 4.00 * STD_DEV, Color.red);
		paintScale.add(MAX, Color.red);

		return paintScale;
	}

	private PaintScaleLegend get_Paint_Scale_Legend(LookupPaintScale paintScale)
	{
		PaintScaleLegend psl = new PaintScaleLegend(paintScale, new NumberAxis());
		psl.setPosition(RectangleEdge.RIGHT);
		psl.setAxisLocation(AxisLocation.TOP_OR_RIGHT);
		psl.setMargin(100.0, 40.0, 160.0, 50.0);
		psl.setStripOutlineStroke(new BasicStroke(10.0f));
		psl.setWidth(200);
		psl.setStripWidth(100);

		return psl;
	}

	private XYZDataset createDataset()
	{
		double[] xvalues = new double[number_of_residues * number_of_residues];
		double[] yvalues = new double[number_of_residues * number_of_residues];
		double[] zvalues = new double[number_of_residues * number_of_residues];

		for (int i = 1; i < number_of_residues + 1; i++)
			{
				for (int j = 1; j < number_of_residues + 1; j++)
					{
						final int idx = (i - 1) * number_of_residues + (j - 1);
						xvalues[idx] = i;
						yvalues[idx] = j;
						zvalues[idx] = coupling_scores.get(i - 1, j - 1);
					}
			}

		DefaultXYZDataset dataset = new DefaultXYZDataset();
		dataset.addSeries("Coupling Scores", new double[][] { xvalues, yvalues, zvalues });

		return dataset;
	}
}