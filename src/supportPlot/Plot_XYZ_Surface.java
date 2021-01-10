package supportPlot;

import java.io.File;

import com.orsoncharts.Chart3D;
import com.orsoncharts.Chart3DFactory;
import com.orsoncharts.Range;
import com.orsoncharts.axis.ValueAxis3D;
import com.orsoncharts.data.function.Function3D;
import com.orsoncharts.data.xyz.XYZDataset;
import com.orsoncharts.data.xyz.XYZSeries;
import com.orsoncharts.data.xyz.XYZSeriesCollection;
import com.orsoncharts.graphics3d.Dimension3D;
import com.orsoncharts.graphics3d.ExportUtils;
import com.orsoncharts.legend.LegendAnchor;
import com.orsoncharts.plot.XYZPlot;
import com.orsoncharts.renderer.ColorScale;
import com.orsoncharts.renderer.RainbowScale;
import com.orsoncharts.renderer.xyz.SurfaceRenderer;
import com.orsoncharts.util.Orientation;

import Jama.Matrix;
import supportIO.Matrix_IO;

public class Plot_XYZ_Surface
{
	final long serialVersionUID = 1L;
	String out_dir, chart_title, xAxis, yAxis, zAxis;
	Matrix data, col1, col2, col3;
	int ROWS;
	XYZSeriesCollection<Double> FES_Data;
	XYZSeries<Double> FES_Series;
	Chart3D chart;

	public Plot_XYZ_Surface(String output_dir, String title, Matrix input_data) throws Exception
	{
		this.out_dir = output_dir;
		this.chart_title = title;
		this.xAxis = "OP1";
		this.yAxis = "FE";
		this.zAxis = "OP2";
		this.data = input_data;
		this.ROWS = data.getRowDimension();
		this.col1 = data.getMatrix(0, ROWS - 1, 0, 0);
		this.col2 = data.getMatrix(0, ROWS - 1, 1, 1);
		this.col3 = data.getMatrix(0, ROWS - 1, 2, 2);

		createDataset(col1, col2, col3);
		createChart();

		ExportUtils.writeAsPNG(chart, 1600, 1600, new File(out_dir + File.separatorChar + chart_title + ".png"));
		// ExportUtils.writeAsPDF(chart, 1200, 800, new File(out_dir + File.separatorChar + chart_title + ".pdf"));
	}


	@SuppressWarnings({ "unchecked", "rawtypes" })
	public XYZDataset<Double> createDataset(Matrix m1, Matrix m2, Matrix m3)
	{
		FES_Data = new XYZSeriesCollection<Double>();
		FES_Series = new XYZSeries("FES-Plot");

		for (int i = 0; i < m1.getRowDimension(); i++)
			{
				double x = m1.get(i, 0);
				double z = m2.get(i, 0);
				double y = m3.get(i, 0);

				System.out.println(x + "\t" + y + "\t" + z);

				FES_Series.add(x, y, z);
			}
		FES_Data.add(FES_Series);

		return FES_Data;
	}


	@SuppressWarnings("serial")
	public Chart3D createChart()
	{
		Function3D function = new Function3D()
		{
			@Override
			public double getValue(double x, double z)
			{

				// double xx = x * 100;
				// double zz = z * 100;

				// double fractionalPartX = xx % 1;
				// double integralPartX = xx - fractionalPartX;
				// double fractionalPartZ = zz % 1;
				// double integralPartY = zz - fractionalPartZ;

				return Math.sin(x) * Math.sin(z);
			}
		};


		chart = Chart3DFactory.createSurfaceChart(chart_title, "Free-Energy Landscape", function, "X", "Y", "Z");
		chart.setLegendPosition(LegendAnchor.BOTTOM_RIGHT, Orientation.VERTICAL);

		XYZPlot plot = (XYZPlot) chart.getPlot();
		plot.setDimensions(new Dimension3D(25, 10, 25));

		ValueAxis3D xAxis = plot.getXAxis();
		xAxis.setRange(-6, 6);

		ValueAxis3D zAxis = plot.getZAxis();
		zAxis.setRange(-6, 6);

		Range r = new Range(-1, 1);
		Range BLUE_TO_RED_RANGE = new Range(0.0, 0.6666);
		ColorScale cs = new RainbowScale(r, 256, BLUE_TO_RED_RANGE);

		SurfaceRenderer renderer = (SurfaceRenderer) plot.getRenderer();

		renderer.setColorScale(new RainbowScale(new Range(-1.0, 1.0)));
		renderer.setColorScale(cs);
		renderer.setDrawFaceOutlines(false);

		return chart;
	}

	public static void main(String[] args) throws Exception
	{
		String path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release2\\test\\BL\\MD\\1LHY\\TEM52\\JEDi_RESULTS_BL_MD_TEM52_Mechanistic_Site\\FES\\Cartesian_All_Atom_PCA\\CORR\\BL_MD_TEM52_Mechanistic_Site_FES_1_2.txt";
		Matrix input_data = Matrix_IO.read_Matrix(path);
		new Plot_XYZ_Surface(System.getProperty("user.dir"), "3D-Plot", input_data);
	}
}
