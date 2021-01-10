package support;

import java.util.ArrayList;

import Jama.Matrix;

public class SSA
{
	int iterations;
	double RMSIP, max_angle;
	Matrix Avg_RMSIP_Score_Matrix, RMSIP_Std_Dev_Matrix, Avg_CO_Score_Matrix, Avg_PA_Matrix, projections, projections_abs, projections_squared, CO_matrix_1_2, CO_matrix_2_1, PAs;
	ArrayList<Double> projections_squared_CO_k, projections_squared_RMSIP, cumulative_overlaps_1_2, cumulative_overlaps_2_1, principle_angles_svd, cosine_products,
			vectorial_sum_of_angles, RMSIPs;

	// ********************************** CONSTRUCTORS *********************************************************** //

	public SSA()
	{
		super();
	}

	// ************************* GETTERS & SETTERS ****************************************************** //

	public int getIterations()
	{
		return iterations;
	}

	public void setIterations(int iterations)
	{
		this.iterations = iterations;
	}

	public double getRMSIP()
	{
		return RMSIP;
	}

	public void setRMSIP(double rMSIP)
	{
		RMSIP = rMSIP;
	}

	public double getMax_angle()
	{
		return max_angle;
	}

	public void setMax_angle(double max_angle)
	{
		this.max_angle = max_angle;
	}

	public Matrix getAvg_RMSIP_Score_Matrix()
	{
		return Avg_RMSIP_Score_Matrix;
	}

	public void setAvg_RMSIP_Score_Matrix(Matrix avg_RMSIP_Score_Matrix)
	{
		Avg_RMSIP_Score_Matrix = avg_RMSIP_Score_Matrix;
	}

	public Matrix getRMSIP_Std_Dev_Matrix()
	{
		return RMSIP_Std_Dev_Matrix;
	}

	public void setRMSIP_Std_Dev_Matrix(Matrix rMSIP_Std_Dev_Matrix)
	{
		RMSIP_Std_Dev_Matrix = rMSIP_Std_Dev_Matrix;
	}

	public Matrix getAvg_CO_Score_Matrix()
	{
		return Avg_CO_Score_Matrix;
	}

	public void setAvg_CO_Score_Matrix(Matrix avg_CO_Score_Matrix)
	{
		Avg_CO_Score_Matrix = avg_CO_Score_Matrix;
	}

	public Matrix getAvg_PA_Matrix()
	{
		return Avg_PA_Matrix;
	}

	public void setAvg_PA_Matrix(Matrix avg_PA_Matrix)
	{
		Avg_PA_Matrix = avg_PA_Matrix;
	}

	public Matrix getProjections()
	{
		return projections;
	}

	public void setProjections(Matrix projections)
	{
		this.projections = projections;
	}

	public Matrix getProjections_abs()
	{
		return projections_abs;
	}

	public void setProjections_abs(Matrix projections_abs)
	{
		this.projections_abs = projections_abs;
	}

	public Matrix getProjections_squared()
	{
		return projections_squared;
	}

	public void setProjections_squared(Matrix projections_squared)
	{
		this.projections_squared = projections_squared;
	}

	public Matrix getCO_matrix_1_2()
	{
		return CO_matrix_1_2;
	}

	public void setCO_matrix_1_2(Matrix cO_matrix_1_2)
	{
		CO_matrix_1_2 = cO_matrix_1_2;
	}

	public Matrix getCO_matrix_2_1()
	{
		return CO_matrix_2_1;
	}

	public void setCO_matrix_2_1(Matrix cO_matrix_2_1)
	{
		CO_matrix_2_1 = cO_matrix_2_1;
	}

	public Matrix getPAs()
	{
		return PAs;
	}

	public void setPAs(Matrix pAs)
	{
		PAs = pAs;
	}

	public ArrayList<Double> getProjections_squared_CO_k()
	{
		return projections_squared_CO_k;
	}

	public void setProjections_squared_CO_k(ArrayList<Double> projections_squared_CO_k)
	{
		this.projections_squared_CO_k = projections_squared_CO_k;
	}

	public ArrayList<Double> getProjections_squared_RMSIP()
	{
		return projections_squared_RMSIP;
	}

	public void setProjections_squared_RMSIP(ArrayList<Double> projections_squared_RMSIP)
	{
		this.projections_squared_RMSIP = projections_squared_RMSIP;
	}

	public ArrayList<Double> getCumulative_overlaps_1_2()
	{
		return cumulative_overlaps_1_2;
	}

	public void setCumulative_overlaps_1_2(ArrayList<Double> cumulative_overlaps_1_2)
	{
		this.cumulative_overlaps_1_2 = cumulative_overlaps_1_2;
	}

	public ArrayList<Double> getCumulative_overlaps_2_1()
	{
		return cumulative_overlaps_2_1;
	}

	public void setCumulative_overlaps_2_1(ArrayList<Double> cumulative_overlaps_2_1)
	{
		this.cumulative_overlaps_2_1 = cumulative_overlaps_2_1;
	}

	public ArrayList<Double> getPrinciple_angles_svd()
	{
		return principle_angles_svd;
	}

	public void setPrinciple_angles_svd(ArrayList<Double> principle_angles_svd)
	{
		this.principle_angles_svd = principle_angles_svd;
	}

	public ArrayList<Double> getCosine_products()
	{
		return cosine_products;
	}

	public void setCosine_products(ArrayList<Double> cosine_products)
	{
		this.cosine_products = cosine_products;
	}

	public ArrayList<Double> getVectorial_sum_of_angles()
	{
		return vectorial_sum_of_angles;
	}

	public void setVectorial_sum_of_angles(ArrayList<Double> vectorial_sum_of_angles)
	{
		this.vectorial_sum_of_angles = vectorial_sum_of_angles;
	}

	public ArrayList<Double> getRMSIPs()
	{
		return RMSIPs;
	}

	public void setRMSIPs(ArrayList<Double> rMSIPs)
	{
		RMSIPs = rMSIPs;
	}
}
