/*
 * 
 */
package supportStats;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import Jama.Matrix;
import supportPlot.PDF_Plot;

/**
 *
 * @author jenny
 */
public class EstimatePDF
{
	// static Plot_XY_Line plot;
	// static Matrix plotMatix;
	// static String outDir, plotTitle;

	final String directory, title;
	final Matrix sampleData;
	final ArrayList<Double> arrayData;

	public EstimatePDF(String directory, String title, ArrayList<Double> sample)
	{
		super();
		this.directory = directory;
		this.title = title;
		this.sampleData = null;
		this.arrayData = sample;
	}

	public EstimatePDF(String directory, String title, Matrix sample)
	{
		super();
		this.directory = directory;
		this.title = title;
		this.sampleData = sample;
		arrayData = new ArrayList<Double>();
		for (int i = 0; i < sampleData.getColumnDimension(); i++)
			{
				arrayData.add(sampleData.get(0, i));
			}
	}

	public Matrix estimate_MAT()
	{
		InputParameters input = new InputParameters();
		input.debug = false;
		InputData data = new InputData(input);
		data.setData(arrayData);

		Score score = new ScoreQZ(input.SURDTarget, input.SURDMinimum, input.SURDMaximum);
		data.out.debug = input.debug;

		MinimizeScore minimumPDF = new MinimizeScore();
		minimumPDF.out.debug = input.debug;
		Results results = new Results(); // moved out of the if clause
		boolean failed = minimumPDF.minimize(input, data, score);
		if (!failed)
			{
				results.out.debug = input.debug;
				results.createSolution(input, data, minimumPDF, score);
				for (int j = 0; j < results.x.size(); j++)
					{
						// writeFile(fileName, results.x, results.PDF, results.CDF);
						plotOutput(directory, title, results.x, results.PDF, results.CDF);
					}
			}
		Matrix margProbs = new Matrix(1, results.px.size());
		for (int i = 0; i < results.px.size(); i++)
			{
				margProbs.set(0, i, results.px.get(i));
			}
		return margProbs;
	}

	public ArrayList<Double> estimate_ALD()
	{
		InputParameters input = new InputParameters();
		input.debug = false;
		InputData data = new InputData(input);
		data.setData(arrayData);

		Score score = new ScoreQZ(input.SURDTarget, input.SURDMinimum, input.SURDMaximum);
		data.out.debug = input.debug;

		MinimizeScore minimumPDF = new MinimizeScore();
		minimumPDF.out.debug = input.debug;
		Results results = new Results(); // moved out of the if clause
		boolean failed = minimumPDF.minimize(input, data, score);
		if (!failed)
			{
				results.out.debug = input.debug;
				results.createSolution(input, data, minimumPDF, score);
				plotOutput(directory, title, results.x, results.PDF, results.CDF);

//				for (int j = 0; j < results.x.size(); j++)
					{
						// writeFile(fileName, results.x, results.PDF, results.CDF);
						// plotOutput(directory, title, results.x, results.PDF, results.CDF);
					}
			}
		return results.px;
	}

	public static void main(String[] args)
	{
		for (int i = 0; i < args.length; i++)
			{

				String fileName = args[i];
				ArrayList<Double> sample = readFile(fileName);
				InputParameters input = new InputParameters();
				input.debug = true;
				InputData data = new InputData(input);
				data.setData(sample);

				Score score = new ScoreQZ(input.SURDTarget, input.SURDMinimum, input.SURDMaximum);
				data.out.debug = input.debug;

				MinimizeScore minimumPDF = new MinimizeScore();
				minimumPDF.out.debug = input.debug;
				boolean failed = minimumPDF.minimize(input, data, score);
				if (!failed)
					{
						Results results = new Results();
						results.out.debug = input.debug;
						results.createSolution(input, data, minimumPDF, score);
						for (int j = 0; j < results.x.size(); j++)
							{
								// writeFile(fileName, results.x, results.PDF, results.CDF);
							}
					}
			}
	}

	static ArrayList<Double> readFile(String fileName)
	{
		FileReader inputFile;

		try
			{
				inputFile = new FileReader(fileName);
			}
		catch (IOException e)
			{
				System.out.println("Error opening input file " + fileName + ": " + e.toString());
				return null;
			}

		ArrayList<Double> sample = new ArrayList<Double>();

		BufferedReader buff = new BufferedReader(inputFile);

		String line;
		while (true)
			{
				try
					{
						line = buff.readLine();
						if (line == null) break;
						sample.add(Double.valueOf(line));
					}
				catch (IOException e)
					{
						System.out.println("" + e);
					}
			}
		try
			{
				buff.close();
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
		return sample;
	}

	void writeFile(String fileName, ArrayList<Double> x, ArrayList<Double> pdf, ArrayList<Double> cdf)
	{
		try
			{
				FileWriter file = new FileWriter("Out_" + fileName);
				BufferedWriter buffer = new BufferedWriter(file);

				for (int i = 0; i < x.size(); i++)
					{
						String line = String.valueOf(x.get(i)) + "  " + String.valueOf(pdf.get(i)) + "  " + String.valueOf(cdf.get(i));
						buffer.write(line);
						buffer.newLine();
					}
				buffer.close();
				file.close();

			}
		catch (IOException e)
			{
				System.err.println("Error writing to outputfile: " + e.toString());
			}
	}

	void plotOutput(String out_dir, String title, ArrayList<Double> x, ArrayList<Double> pdf, ArrayList<Double> cdf)
	{
		Matrix plotMatrix = new Matrix(x.size(), 3);
		for (int i = 0; i < x.size(); i++)
			{
				plotMatrix.set(i, 0, x.get(i));
				plotMatrix.set(i, 1, pdf.get(i));
				plotMatrix.set(i, 2, cdf.get(i));
			}
		PDF_Plot.createChart2Series(out_dir, title, plotMatrix);
	}
}
