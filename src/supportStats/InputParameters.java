/*
 * 
 */
package supportStats;

/**
 *
 * @author jenny
 */
public class InputParameters
{

	boolean debug = false;

	float lowerBound;
	float upperBound;
	boolean lowerBoundSpecified = false;
	boolean upperBoundSpecified = false;

	String scoreType = "QZ";
	double SURDMinimum = 5;
	double SURDTarget = 40;
	double SURDMaximum = 100;
	int initPartitionSize = 1025;
	int integrationPoints = -1;
	int maxLagrange = 100;
	int minLagrange = 1;
	int nLagrangeAdd = 5;
	double outlierCutoff = 7.0;

	double fractionLagrangeAdd = 0.1;
	double initSigma = 0.1;
	double finalSigma = 0.001;
	double decayFactor = Math.sqrt(2.0);
	int loopMax = 30;
}
