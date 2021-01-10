/*
 * 
 */
package supportStats;

/**
 *
 * @author jenny
 */
public class OutputControl
{
	boolean debug;

	public void print(String output)
	{
		if (debug)
			{
				System.out.println(output);
			}
	}

	public void print(String output, int value)
	{
		if (debug)
			{
				System.out.println(output + ": " + value);
			}
	}

	public void print(String output, double value)
	{
		if (debug)
			{
				System.out.println(output + ": " + value);
			}
	}

	public void error(String output)
	{
		System.err.println(output);
	}

	public void error(String output, int value)
	{
		System.err.println(output + ": " + value);
	}

	public void error(String output, double value)
	{
		System.err.println(output + ": " + value);
	}
}
