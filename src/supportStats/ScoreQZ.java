/*
 * To change this license header, choose License Headers in Project Properties. To change this template file, choose Tools | Templates and open the template in the editor.
 */
package supportStats;

//import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author jenny
 */
public class ScoreQZ extends Score
{
	ScoreQZ(double confidenceTarget, double confidenceMin, double confidenceMax)
	{
		getValues();
		targetScore = getTargetScore(confidenceTarget);
		minimumScore = getTargetScore(confidenceMin);
		maximumScore = getTargetScore(confidenceMax);
	}

	@Override
	public ArrayList<Integer> getIndices(int N, int p)
	{
		ArrayList<Integer> indicesLocal = new ArrayList<Integer>();
		if (N == p)
			{
				for (int i = 0; i < p; i++)
					{
						indicesLocal.add(i, i);
					}
			}
		else
			{
				indicesLocal.add(0, 0);
				double div = N * 1.0 / (p - 1);
				double index = div - 1;
				for (int i = 1; i < p; i++)
					{
						indicesLocal.add(i, (int) Math.round(index));
						index += div;
					}
				indicesLocal.add(p - 1, N - 1);
			}
		return indicesLocal;
	}

	@Override
	public ArrayList<Integer> setIndices(int N, int p)
	{
		indices = getIndices(N, p);
		return indices;
	}

	@Override
	public double calculateScore(double r[], int N, int p)
	{
		double mu;
		double QZ = 0;
		QZVariance = 0;

		for (int i = 0; i < p; i++)
			{
				mu = (indices.get(i) + 1.0) / (N + 1);
				double difference = mu - r[i];
				double qz = difference * difference * (N + 2) / (mu * (1 - mu));
				QZ += qz;
				QZVariance += qz * qz;
			}

		QZ /= p;
		QZVariance /= p;
		likelihood = -QZ;
		SURD = -QZ;

		return -QZVariance;
	}

	public void getValues()
	{
		double arrayScore[] = { -6.474, -6.4649, -6.4557, -6.4466, -6.4375, -6.4283, -6.4192, -6.4101, -6.4009, -6.3918, -6.3827, -6.3735, -6.3644, -6.3553, -6.3462, -6.337,
				-6.3279, -6.3188, -6.3096, -6.3005, -6.2914, -6.2822, -6.2731, -6.264, -6.2548, -6.2457, -6.2366, -6.2274, -6.2183, -6.2092, -6.2, -6.1909, -6.1818, -6.1727,
				-6.1635, -6.1544, -6.1453, -6.1361, -6.127, -6.1179, -6.1087, -6.0996, -6.0905, -6.0813, -6.0722, -6.0631, -6.0539, -6.0448, -6.0357, -6.0265, -6.0174, -6.0083,
				-5.9991, -5.99, -5.9809, -5.9718, -5.9626, -5.9535, -5.9444, -5.9352, -5.9261, -5.917, -5.9078, -5.8987, -5.8896, -5.8804, -5.8713, -5.8622, -5.853, -5.8439,
				-5.8348, -5.8256, -5.8165, -5.8074, -5.7983, -5.7891, -5.78, -5.7709, -5.7617, -5.7526, -5.7435, -5.7343, -5.7252, -5.7161, -5.7069, -5.6978, -5.6887, -5.6795,
				-5.6704, -5.6613, -5.6521, -5.643, -5.6339, -5.6248, -5.6156, -5.6065, -5.5974, -5.5882, -5.5791, -5.57, -5.5608, -5.5517, -5.5426, -5.5334, -5.5243, -5.5152,
				-5.506, -5.4969, -5.4878, -5.4786, -5.4695, -5.4604, -5.4513, -5.4421, -5.433, -5.4239, -5.4147, -5.4056, -5.3965, -5.3873, -5.3782, -5.3691, -5.3599, -5.3508,
				-5.3417, -5.3325, -5.3234, -5.3143, -5.3051, -5.296, -5.2869, -5.2777, -5.2686, -5.2595, -5.2504, -5.2412, -5.2321, -5.223, -5.2138, -5.2047, -5.1956, -5.1864,
				-5.1773, -5.1682, -5.159, -5.1499, -5.1408, -5.1316, -5.1225, -5.1134, -5.1042, -5.0951, -5.086, -5.0769, -5.0677, -5.0586, -5.0495, -5.0403, -5.0312, -5.0221,
				-5.0129, -5.0038, -4.9947, -4.9855, -4.9764, -4.9673, -4.9581, -4.949, -4.9399, -4.9307, -4.9216, -4.9125, -4.9034, -4.8942, -4.8851, -4.876, -4.8668, -4.8577,
				-4.8486, -4.8394, -4.8303, -4.8212, -4.812, -4.8029, -4.7938, -4.7846, -4.7755, -4.7664, -4.7572, -4.7481, -4.739, -4.7299, -4.7207, -4.7116, -4.7025, -4.6933,
				-4.6842, -4.6751, -4.6659, -4.6568, -4.6477, -4.6385, -4.6294, -4.6203, -4.6111, -4.602, -4.5929, -4.5837, -4.5746, -4.5655, -4.5563, -4.5472, -4.5381, -4.529,
				-4.5198, -4.5107, -4.5016, -4.4924, -4.4833, -4.4742, -4.465, -4.4559, -4.4468, -4.4376, -4.4285, -4.4194, -4.4102, -4.4011, -4.392, -4.3828, -4.3737, -4.3646,
				-4.3555, -4.3463, -4.3372, -4.3281, -4.3189, -4.3098, -4.3007, -4.2915, -4.2824, -4.2733, -4.2641, -4.255, -4.2459, -4.2367, -4.2276, -4.2185, -4.2093, -4.2002,
				-4.1911, -4.182, -4.1728, -4.1637, -4.1546, -4.1454, -4.1363, -4.1272, -4.118, -4.1089, -4.0998, -4.0906, -4.0815, -4.0724, -4.0632, -4.0541, -4.045, -4.0358,
				-4.0267, -4.0176, -4.0085, -3.9993, -3.9902, -3.9811, -3.9719, -3.9628, -3.9537, -3.9445, -3.9354, -3.9263, -3.9171, -3.908, -3.8989, -3.8897, -3.8806, -3.8715,
				-3.8623, -3.8532, -3.8441, -3.835, -3.8258, -3.8167, -3.8076, -3.7984, -3.7893, -3.7802, -3.771, -3.7619, -3.7528, -3.7436, -3.7345, -3.7254, -3.7162, -3.7071,
				-3.698, -3.6888, -3.6797, -3.6706, -3.6614, -3.6523, -3.6432, -3.6341, -3.6249, -3.6158, -3.6067, -3.5975, -3.5884, -3.5793, -3.5701, -3.561, -3.5519, -3.5427,
				-3.5336, -3.5245, -3.5153, -3.5062, -3.4971, -3.4879, -3.4788, -3.4697, -3.4606, -3.4514, -3.4423, -3.4332, -3.424, -3.4149, -3.4058, -3.3966, -3.3875, -3.3784,
				-3.3692, -3.3601, -3.351, -3.3418, -3.3327, -3.3236, -3.3144, -3.3053, -3.2962, -3.2871, -3.2779, -3.2688, -3.2597, -3.2505, -3.2414, -3.2323, -3.2231, -3.214,
				-3.2049, -3.1957, -3.1866, -3.1775, -3.1683, -3.1592, -3.1501, -3.1409, -3.1318, -3.1227, -3.1136, -3.1044, -3.0953, -3.0862, -3.077, -3.0679, -3.0588, -3.0496,
				-3.0405, -3.0314, -3.0222, -3.0131, -3.004, -2.9948, -2.9857, -2.9766, -2.9674, -2.9583, -2.9492, -2.94, -2.9309, -2.9218, -2.9127, -2.9035, -2.8944, -2.8853,
				-2.8761, -2.867, -2.8579, -2.8487, -2.8396, -2.8305, -2.8213, -2.8122, -2.8031, -2.7939, -2.7848, -2.7757, -2.7665, -2.7574, -2.7483, -2.7392, -2.73, -2.7209,
				-2.7118, -2.7026, -2.6935, -2.6844, -2.6752, -2.6661, -2.657, -2.6478, -2.6387, -2.6296, -2.6204, -2.6113, -2.6022, -2.593, -2.5839, -2.5748, -2.5657, -2.5565,
				-2.5474, -2.5383, -2.5291, -2.52, -2.5109, -2.5017, -2.4926, -2.4835, -2.4743, -2.4652, -2.4561, -2.4469, -2.4378, -2.4287, -2.4195, -2.4104, -2.4013, -2.3922,
				-2.383, -2.3739, -2.3648, -2.3556, -2.3465, -2.3374, -2.3282, -2.3191, -2.31, -2.3008, -2.2917, -2.2826, -2.2734, -2.2643, -2.2552, -2.246, -2.2369, -2.2278,
				-2.2187, -2.2095, -2.2004, -2.1913, -2.1821, -2.173, -2.1639, -2.1547, -2.1456, -2.1365, -2.1273, -2.1182, -2.1091, -2.0999, -2.0908, -2.0817, -2.0725, -2.0634,
				-2.0543, -2.0451, -2.036, -2.0269, -2.0178, -2.0086, -1.9995, -1.9904, -1.9812, -1.9721, -1.963, -1.9538, -1.9447, -1.9356, -1.9264, -1.9173, -1.9082, -1.899,
				-1.8899, -1.8808, -1.8716, -1.8625, -1.8534, -1.8443, -1.8351, -1.826, -1.8169, -1.8077, -1.7986, -1.7895, -1.7803, -1.7712, -1.7621, -1.7529, -1.7438, -1.7347,
				-1.7255, -1.7164, -1.7073, -1.6981, -1.689, -1.6799, -1.6708, -1.6616, -1.6525, -1.6434, -1.6342, -1.6251, -1.616, -1.6068, -1.5977, -1.5886, -1.5794, -1.5703,
				-1.5612, -1.552, -1.5429, -1.5338, -1.5246, -1.5155, -1.5064, -1.4973, -1.4881, -1.479, -1.4699, -1.4607, -1.4516, -1.4425, -1.4333, -1.4242, -1.4151, -1.4059,
				-1.3968, -1.3877, -1.3785, -1.3694, -1.3603, -1.3511, -1.342, -1.3329, -1.3237, -1.3146, -1.3055, -1.2964, -1.2872, -1.2781, -1.269, -1.2598, -1.2507, -1.2416,
				-1.2324, -1.2233, -1.2142, -1.205, -1.1959, -1.1868, -1.1776, -1.1685, -1.1594, -1.1502, -1.1411, -1.132, -1.1229, -1.1137, -1.1046, -1.0955, -1.0863, -1.0772,
				-1.0681, -1.0589, -1.0498, -1.0407, -1.0315, -1.0224, -1.0133, -1.0041, -0.99501, -0.98588, -0.97675, -0.96762, -0.95848, -0.94935, -0.94022, -0.93109, -0.92196,
				-0.91283, -0.90369, -0.89456, -0.88543, -0.8763, -0.86717, -0.85804, -0.8489, -0.83977, -0.83064, -0.82151, -0.81238, -0.80325, -0.79411, -0.78498, -0.77585,
				-0.76672, -0.75759, -0.74846, -0.73932, -0.73019, -0.72106, -0.71193, -0.7028, -0.69367, -0.68454, -0.6754, -0.66627, -0.65714, -0.64801, -0.63888, -0.62975,
				-0.62061, -0.61148, -0.60235, -0.59322, -0.58409, -0.57496, -0.56582, -0.55669, -0.54756, -0.53843, -0.5293, -0.52017, -0.51103, -0.5019, -0.49277, -0.48364,
				-0.47451, -0.46538, -0.45624, -0.44711, -0.43798, -0.42885, -0.41972, -0.41059, -0.40145, -0.39232, -0.38319, -0.37406, -0.36493, -0.3558, -0.34666, -0.33753,
				-0.3284, -0.31927, -0.31014, -0.30101, -0.29187, -0.28274, -0.27361, -0.26448, -0.25535, -0.24622, -0.23708, -0.22795, -0.21882, -0.20969, -0.20056, -0.19143,
				-0.18229, -0.17316, -0.16403, -0.1549, -0.14577, -0.13664, -0.12751, -0.11837, -0.10924, -0.10011, -0.090979, -0.081847 };
		double arraySurd[] = { 1.7764e-15, 1.761e-06, 3.7079e-06, 5.8413e-06, 8.1782e-06, 1.0748e-05, 1.3589e-05, 1.6745e-05, 2.0266e-05, 2.4202e-05, 2.8601e-05, 3.3513e-05,
				3.8979e-05, 4.5038e-05, 5.1718e-05, 5.9042e-05, 6.7022e-05, 7.5664e-05, 8.4965e-05, 9.4915e-05, 0.0001055, 0.0001167, 0.0001285, 0.00014087, 0.00015379, 0.00016722,
				0.00018116, 0.00019556, 0.00021042, 0.0002257, 0.00024139, 0.00025746, 0.00027389, 0.00029066, 0.00030775, 0.00032514, 0.00034281, 0.00036074, 0.00037891,
				0.00039729, 0.00041586, 0.00043462, 0.00045352, 0.00047256, 0.00049171, 0.00051095, 0.00053026, 0.00054963, 0.00056903, 0.00058845, 0.00060788, 0.00062728,
				0.00064667, 0.00066601, 0.00068529, 0.00070452, 0.00072367, 0.00074275, 0.00076173, 0.00078063, 0.00079942, 0.00081812, 0.00083671, 0.0008552, 0.00087359,
				0.00089188, 0.00091006, 0.00092814, 0.00094613, 0.00096403, 0.00098184, 0.00099957, 0.0010172, 0.0010348, 0.0010523, 0.0010698, 0.0010872, 0.0011046, 0.0011219,
				0.0011392, 0.0011565, 0.0011738, 0.001191, 0.0012083, 0.0012255, 0.0012428, 0.00126, 0.0012773, 0.0012946, 0.001312, 0.0013293, 0.0013467, 0.0013641, 0.0013816,
				0.0013991, 0.0014166, 0.0014341, 0.0014517, 0.0014692, 0.0014868, 0.0015044, 0.0015221, 0.0015397, 0.0015573, 0.001575, 0.0015926, 0.0016102, 0.0016278, 0.0016453,
				0.0016628, 0.0016803, 0.0016977, 0.0017151, 0.0017325, 0.0017497, 0.0017669, 0.0017841, 0.0018011, 0.0018181, 0.001835, 0.0018519, 0.0018686, 0.0018853, 0.0019019,
				0.0019184, 0.0019349, 0.0019513, 0.0019676, 0.0019838, 0.002, 0.0020161, 0.0020321, 0.0020481, 0.002064, 0.0020799, 0.0020957, 0.0021116, 0.0021273, 0.0021431,
				0.0021588, 0.0021745, 0.0021902, 0.0022059, 0.0022216, 0.0022373, 0.0022531, 0.0022688, 0.0022845, 0.0023003, 0.0023161, 0.0023319, 0.0023478, 0.0023637, 0.0023797,
				0.0023957, 0.0024117, 0.0024278, 0.0024439, 0.0024602, 0.0024764, 0.0024927, 0.0025091, 0.0025256, 0.0025421, 0.0025586, 0.0025753, 0.002592, 0.0026088, 0.0026256,
				0.0026426, 0.0026596, 0.0026766, 0.0026938, 0.002711, 0.0027283, 0.0027457, 0.0027632, 0.0027808, 0.0027984, 0.0028162, 0.0028341, 0.002852, 0.0028701, 0.0028882,
				0.0029065, 0.0029249, 0.0029435, 0.0029621, 0.0029809, 0.0029998, 0.0030189, 0.0030381, 0.0030575, 0.003077, 0.0030967, 0.0031166, 0.0031367, 0.003157, 0.0031774,
				0.0031981, 0.003219, 0.0032401, 0.0032615, 0.0032831, 0.003305, 0.0033271, 0.0033495, 0.0033722, 0.0033951, 0.0034184, 0.003442, 0.0034659, 0.0034901, 0.0035146,
				0.0035396, 0.0035648, 0.0035905, 0.0036165, 0.0036429, 0.0036697, 0.0036969, 0.0037246, 0.0037526, 0.0037812, 0.0038101, 0.0038396, 0.0038695, 0.0038999, 0.0039309,
				0.0039623, 0.0039943, 0.0040268, 0.0040599, 0.0040935, 0.0041278, 0.0041626, 0.0041981, 0.0042341, 0.0042709, 0.0043083, 0.0043464, 0.0043852, 0.0044248, 0.0044651,
				0.0045061, 0.004548, 0.0045907, 0.0046342, 0.0046785, 0.0047238, 0.0047699, 0.004817, 0.0048651, 0.0049142, 0.0049642, 0.0050154, 0.0050676, 0.0051209, 0.0051753,
				0.005231, 0.0052878, 0.0053459, 0.0054052, 0.0054659, 0.0055278, 0.0055912, 0.005656, 0.0057222, 0.0057898, 0.005859, 0.0059298, 0.0060021, 0.006076, 0.0061516,
				0.0062288, 0.0063078, 0.0063885, 0.006471, 0.0065553, 0.0066414, 0.0067294, 0.0068193, 0.0069111, 0.0070048, 0.0071006, 0.0071983, 0.007298, 0.0073998, 0.0075037,
				0.0076097, 0.0077177, 0.0078279, 0.0079402, 0.0080547, 0.0081713, 0.0082901, 0.0084111, 0.0085343, 0.0086597, 0.0087874, 0.0089173, 0.0090494, 0.0091837, 0.0093204,
				0.0094592, 0.0096004, 0.0097438, 0.0098895, 0.010037, 0.010188, 0.01034, 0.010495, 0.010652, 0.010811, 0.010973, 0.011137, 0.011303, 0.011472, 0.011643, 0.011816,
				0.011991, 0.012169, 0.012349, 0.012531, 0.012716, 0.012902, 0.013092, 0.013283, 0.013477, 0.013673, 0.013871, 0.014072, 0.014274, 0.01448, 0.014687, 0.014896,
				0.015108, 0.015322, 0.015539, 0.015757, 0.015978, 0.016201, 0.016426, 0.016654, 0.016883, 0.017115, 0.017349, 0.017585, 0.017823, 0.018064, 0.018306, 0.018551,
				0.018797, 0.019046, 0.019297, 0.01955, 0.019805, 0.020062, 0.020321, 0.020582, 0.020845, 0.02111, 0.021378, 0.021647, 0.021918, 0.022192, 0.022467, 0.022745,
				0.023025, 0.023306, 0.02359, 0.023876, 0.024165, 0.024455, 0.024748, 0.025043, 0.02534, 0.025639, 0.025941, 0.026245, 0.026552, 0.026861, 0.027173, 0.027487,
				0.027804, 0.028124, 0.028446, 0.028771, 0.029099, 0.02943, 0.029764, 0.030101, 0.030441, 0.030784, 0.031131, 0.031481, 0.031834, 0.032191, 0.032551, 0.032915,
				0.033282, 0.033654, 0.034029, 0.034408, 0.034791, 0.035178, 0.035569, 0.035964, 0.036364, 0.036768, 0.037176, 0.037589, 0.038006, 0.038428, 0.038855, 0.039286,
				0.039722, 0.040163, 0.040609, 0.041059, 0.041515, 0.041976, 0.042442, 0.042914, 0.04339, 0.043872, 0.04436, 0.044853, 0.045351, 0.045855, 0.046364, 0.04688,
				0.047401, 0.047927, 0.04846, 0.048999, 0.049543, 0.050094, 0.050651, 0.051214, 0.051783, 0.052359, 0.052941, 0.05353, 0.054125, 0.054727, 0.055336, 0.055951,
				0.056574, 0.057203, 0.05784, 0.058484, 0.059135, 0.059794, 0.060461, 0.061135, 0.061817, 0.062507, 0.063205, 0.063911, 0.064626, 0.065349, 0.066081, 0.066821,
				0.067571, 0.068329, 0.069097, 0.069874, 0.070661, 0.071458, 0.072265, 0.073081, 0.073908, 0.074746, 0.075594, 0.076453, 0.077323, 0.078204, 0.079097, 0.080001,
				0.080918, 0.081846, 0.082786, 0.083739, 0.084705, 0.085683, 0.086674, 0.087679, 0.088697, 0.089729, 0.090774, 0.091834, 0.092907, 0.093995, 0.095098, 0.096215,
				0.097348, 0.098495, 0.099658, 0.10084, 0.10203, 0.10324, 0.10447, 0.10571, 0.10697, 0.10824, 0.10953, 0.11084, 0.11216, 0.1135, 0.11486, 0.11624, 0.11763, 0.11904,
				0.12047, 0.12192, 0.12338, 0.12486, 0.12636, 0.12788, 0.12942, 0.13098, 0.13255, 0.13415, 0.13576, 0.1374, 0.13905, 0.14073, 0.14243, 0.14414, 0.14588, 0.14765,
				0.14943, 0.15124, 0.15307, 0.15493, 0.15681, 0.15872, 0.16065, 0.16261, 0.1646, 0.16661, 0.16866, 0.17073, 0.17284, 0.17498, 0.17715, 0.17935, 0.18159, 0.18386,
				0.18617, 0.18852, 0.1909, 0.19333, 0.19579, 0.19829, 0.20084, 0.20342, 0.20605, 0.20873, 0.21144, 0.21421, 0.21702, 0.21987, 0.22277, 0.22572, 0.22871, 0.23175,
				0.23484, 0.23798, 0.24116, 0.2444, 0.24768, 0.251, 0.25438, 0.2578, 0.26127, 0.26479, 0.26836, 0.27197, 0.27563, 0.27934, 0.28309, 0.28689, 0.29074, 0.29464,
				0.29859, 0.30258, 0.30663, 0.31073, 0.31487, 0.31908, 0.32333, 0.32764, 0.33201, 0.33644, 0.34093, 0.34548, 0.35009, 0.35478, 0.35953, 0.36435, 0.36925, 0.37423,
				0.37928, 0.38442, 0.38964, 0.39494, 0.40034, 0.40582, 0.4114, 0.41707, 0.42283, 0.42869, 0.43465, 0.44071, 0.44686, 0.45311, 0.45945, 0.46589, 0.47243, 0.47906,
				0.48579, 0.4926, 0.49951, 0.5065, 0.51357, 0.52073, 0.52798, 0.5353, 0.5427, 0.55018, 0.55774, 0.56537, 0.57309, 0.58087, 0.58874, 0.59669, 0.60471, 0.61282,
				0.62101, 0.62929, 0.63765, 0.6461, 0.65464, 0.66328, 0.672, 0.68082, 0.68972, 0.69872, 0.70779, 0.71695, 0.72617, 0.73546, 0.7448, 0.75419, 0.7636, 0.77303,
				0.78245, 0.79186, 0.80124, 0.81056, 0.81982, 0.829, 0.83809, 0.84706, 0.8559, 0.86462, 0.87319, 0.8816, 0.88986, 0.89795, 0.90585, 0.91357, 0.92109, 0.9284,
				0.93547, 0.94229, 0.94882, 0.95505, 0.96093, 0.96642, 0.97149, 0.97611, 0.98025, 0.9839, 0.98704, 0.98971, 0.99192, 0.99371, 0.99513, 0.99625, 0.99711, 0.99777,
				0.99828, 0.99867, 0.99898, 0.99923, 0.99945, 0.99964, 0.99982 };
		for (int i = 0; i < 701; i++)
			{
				scores.add(arrayScore[i]);
				SURDs.add(arraySurd[i]);
			}
	}
}
