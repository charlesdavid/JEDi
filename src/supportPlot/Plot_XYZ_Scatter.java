package supportPlot;

import java.awt.Color;
import java.io.File;

import com.orsoncharts.Chart3D;
import com.orsoncharts.Chart3DFactory;
import com.orsoncharts.Colors;
import com.orsoncharts.axis.ValueAxis3D;
import com.orsoncharts.data.xyz.XYZDataset;
import com.orsoncharts.data.xyz.XYZSeries;
import com.orsoncharts.data.xyz.XYZSeriesCollection;
import com.orsoncharts.graphics3d.Dimension3D;
import com.orsoncharts.graphics3d.ExportUtils;
import com.orsoncharts.graphics3d.ViewPoint3D;
import com.orsoncharts.legend.LegendAnchor;
import com.orsoncharts.plot.XYZPlot;
import com.orsoncharts.renderer.xyz.ScatterXYZRenderer;
import com.orsoncharts.util.Orientation;

import Jama.Matrix;
import supportIO.Matrix_IO;


public class Plot_XYZ_Scatter
{
	String out_dir, chart_title, chart_subtitle, xAxis, yAxis, zAxis;
	Matrix data, col1, col2, col3;
	int ROWS;
	double maxX, minX, maxZ, minZ;
	XYZSeriesCollection<Double> FES_Data;
	XYZSeries<Double> FES_Series;
	Chart3D chart;

	public Plot_XYZ_Scatter(String output_dir, String title, String subtitle, String X, String Y, String Z, Matrix input_data)
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.chart_subtitle = subtitle;
		this.xAxis = X;
		this.yAxis = Y;
		this.zAxis = Z;
		this.data = input_data;

		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);
		this.col3 = data.getMatrix(0, ROWS - 1, 2, 2);

		createDataset(col1, col2, col3);
		createChart();

		try
			{
				ExportUtils.writeAsPNG(chart, 1600, 1600, new File(out_dir + File.separatorChar + chart_title + ".png"));
			}
		catch (Exception e)
			{
				System.err.println("Sorry, there are problems writing the file." + out_dir + File.separatorChar + chart_title + ".png");
				e.printStackTrace();
			}
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public XYZDataset createDataset(Matrix m1, Matrix m2, Matrix m3)
	{
		FES_Data = new XYZSeriesCollection();
		FES_Series = new XYZSeries("FES-Plot");

		maxX = 0;
		minX = 0;
		maxZ = 0;
		minZ = 0;

		for (int i = 0; i < m1.getRowDimension(); i++)
			{
				double x = m1.get(i, 0);
				if (x > maxX) maxX = x;
				if (x < minX) minX = x;

				double z = m2.get(i, 0);
				if (z > maxX) maxZ = z;
				if (z < minZ) minZ = z;

				double y = m3.get(i, 0);

				// System.out.println(x + "\t" + y + "\t" + z);

				FES_Series.add(x, y, z);
			}
		FES_Data.add(FES_Series);
		return FES_Data;
	}


	public Chart3D createChart()
	{
		chart = Chart3DFactory.createScatterChart(chart_title, "Free-Energy Landscape", FES_Data, "X", "Y", "Z");
		chart.setLegendPosition(LegendAnchor.BOTTOM_RIGHT, Orientation.VERTICAL);
		chart.setAntiAlias(true);

		XYZPlot plot = (XYZPlot) chart.getPlot();

		ValueAxis3D xAxis = plot.getXAxis();
		xAxis.setRange(minX, maxX);

		ValueAxis3D zAxis = plot.getZAxis();
		zAxis.setRange(minZ, maxZ);

		plot.setDimensions(new Dimension3D(20, 20, 20));

		ScatterXYZRenderer renderer = (ScatterXYZRenderer) plot.getRenderer();
		renderer.setSize(0.15);
		renderer.setColors(Colors.createIntenseColors());
		renderer.setColors(new Color(0, 0, 255), new Color(0, 255, 0));
		// renderer.setItemLabelBackgroundColor(Color.WHITE);
		// renderer.setItemLabelColor(Color.RED);
		chart.setViewPoint(ViewPoint3D.createAboveLeftViewPoint(40));

		return chart;
	}

	public static void main(String[] args)
	{
		String path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release2\\test\\BL\\MD\\1LHY\\TEM52\\JEDi_RESULTS_residues1\\FES\\Cartesian_All_Atom_PCA\\CORR\\residues1_FES_1_2.txt";
		Matrix input_data = Matrix_IO.read_Matrix(path);

		new Plot_XYZ_Scatter(System.getProperty("user.dir"), "Test-3D-Plot", "Free Energy Landscape", "X", "Y", "Z", input_data);
	}
}
