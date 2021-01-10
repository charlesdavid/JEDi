package jedi;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import support.Outlier_Processing;
import support.Transform_Coordinates;
import supportIO.Input_Parameters;

/**
 * JED class JED_Get_Transformed_Coordinates: Optimally aligns all frames of a trajectory to a reference frame using Quaternion operations. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class JEDi_Get_Transformed_Coordinates
{
	boolean exist, success;
	final int ROWS, COLS, number_of_atoms;
	double z_cut_off, mad_cut_off;
	List<Double> original_conformation_rmsds, aligned_conformation_rmsds, aligned_atomic_rmsf_list, Z_Scores;
	Matrix aligned_subset_REF_PDB_coordinates, aligned_subset_PDB_coordinates, adjusted_PDB_coordinates_rows, coordinates_Outliers_REMOVED, coordinates_Outliers_SELECTED,
			adjustments_per_variable_REMOVE_MAD, adjustments_per_variable_REMOVE_Z, adjustments_per_variable_SELECT_MAD, adjustments_per_variable_SELECT_Z, coordinate_stats;
	final Matrix subset_REF_PDB_coordinates, subset_PDB_coordinates;
	Transform_Coordinates tfc;
	Outlier_Processing op;

	/* ********************************************* CONSTRUCTORS ************************************************************************** */

	/**
	 * Constructor for aligning the PDB coordinates using quaternion operations and calculating the conformation RMSDs and atomic RMSFs.
	 * 
	 * Outliers are detected using a Z-cutoff threshold, and their values replaced with the variable mean value.
	 * 
	 * This class is used for Cartesian and Displacement type analyses, but not the distance analysis, which uses internal coordinates (distances)
	 * 
	 * @param coords     The PDB coordinates
	 * @param ref_coords The reference PDB coordinates to use for the transformation
	 */
	public JEDi_Get_Transformed_Coordinates(Matrix coords, Matrix ref_coords)
	{
		this.subset_PDB_coordinates = coords;
		this.subset_REF_PDB_coordinates = ref_coords;

		this.ROWS = subset_PDB_coordinates.getRowDimension();
		this.COLS = subset_PDB_coordinates.getColumnDimension();
		this.number_of_atoms = (ROWS / 3);

		this.mad_cut_off = Input_Parameters.MAD_SCORE_CUTOFF;
		this.z_cut_off = Input_Parameters.Z_SCORE_CUTOFF;

		aligned_atomic_rmsf_list = new ArrayList<>(ROWS);
		original_conformation_rmsds = new ArrayList<>(COLS);
		aligned_conformation_rmsds = new ArrayList<>(COLS);
		aligned_subset_PDB_coordinates = new Matrix(ROWS, COLS);

		tfc = new Transform_Coordinates(subset_REF_PDB_coordinates, subset_REF_PDB_coordinates);
		aligned_subset_REF_PDB_coordinates = tfc.get_transformed_coords();
	}

	/* ************************************** SETTERS ******************************************************************************** */

	public void set_MAD_Score_Cutoff(double mad_score)
	{
		this.mad_cut_off = mad_score;
	}

	public void set_Z_Score_Cutoff(double zscore)
	{
		this.z_cut_off = zscore;
	}

	/* ************************************** METHODS ******************************************************************************** */

	public Matrix get_SS_Transformed_coords()
	{
		for (int i = 0; i < COLS; i++)
			{
				Matrix fc = subset_PDB_coordinates.getMatrix(0, ROWS - 1, i, i);
				tfc = new Transform_Coordinates(aligned_subset_REF_PDB_coordinates, fc);
				Matrix transformed_col_vector = tfc.get_transformed_coords();
				aligned_subset_PDB_coordinates.setMatrix(0, ROWS - 1, i, i, transformed_col_vector);
				tfc = null;
			}
		return aligned_subset_PDB_coordinates;
	}

	public Matrix get_SS_coordinate_STATS()
	{
		op = new Outlier_Processing(aligned_subset_PDB_coordinates);
		op.set_Mad_threshold(mad_cut_off);
		op.set_Z_threshold(z_cut_off);

		return coordinate_stats = op.getCoordinateStats();
	}

	public Matrix get_SS_transformed_coordinates_OUTLIERS_REMOVED()
	{
		if (mad_cut_off > 0)
			{
				coordinates_Outliers_REMOVED = op.remove_outliers_MAD();
				adjustments_per_variable_REMOVE_MAD = op.getCounts_remove_MAD();
			}
		if (z_cut_off > 0)
			{
				coordinates_Outliers_REMOVED = op.remove_outliers_Z();
				adjustments_per_variable_REMOVE_Z = op.getCounts_remove_Z();
			}

		return coordinates_Outliers_REMOVED;
	}

	public Matrix get_SS_transformed_coordinates_OUTLIERS_SELECTED()
	{
		if (mad_cut_off > 0)
			{
				coordinates_Outliers_SELECTED = op.select_outliers_MAD();
				adjustments_per_variable_SELECT_MAD = op.getCounts_select_MAD();
			}
		if (z_cut_off > 0)
			{
				coordinates_Outliers_SELECTED = op.select_outliers_Z();
				adjustments_per_variable_SELECT_Z = op.getCounts_select_Z();
			}

		return coordinates_Outliers_SELECTED;
	}

	public List<Double> get_SS_Conformation_RMSDs()
	{
		for (int i = 0; i < COLS; i++)
			{
				Matrix fc_O = subset_PDB_coordinates.getMatrix(0, ROWS - 1, i, i);
				Matrix fc_T = aligned_subset_PDB_coordinates.getMatrix(0, ROWS - 1, i, i);

				JEDi_Get_Conformation_RMSD gRMSD_O = new JEDi_Get_Conformation_RMSD(subset_REF_PDB_coordinates, fc_O);
				JEDi_Get_Conformation_RMSD gRMSD_T = new JEDi_Get_Conformation_RMSD(aligned_subset_REF_PDB_coordinates, fc_T);

				double rmsd_O = gRMSD_O.get_RMSD();
				original_conformation_rmsds.add(rmsd_O);
				double rmsd_T = gRMSD_T.get_RMSD();
				aligned_conformation_rmsds.add(rmsd_T);
			}
		return aligned_conformation_rmsds;
	}

	public List<Double> get_SS_RMSF()
	{
		JEDi_Get_Atomic_RMSFs rrmsd = new JEDi_Get_Atomic_RMSFs(aligned_subset_PDB_coordinates);
		aligned_atomic_rmsf_list = rrmsd.get_residue_rmsfs();
		return aligned_atomic_rmsf_list;
	}

	/* ************************************** GETTERS ******************************************************************************** */

	public Matrix get_Original_reference_coordinates()
	{
		return subset_REF_PDB_coordinates;
	}

	public Matrix get_Aligned_Reference_Coordinates()
	{
		return aligned_subset_REF_PDB_coordinates;
	}

	public Matrix getStats()
	{
		return coordinate_stats;
	}

	public Matrix getCoordinates_Outliers_REMOVED()
	{
		return coordinates_Outliers_REMOVED;
	}

	public Matrix getCoordinates_Outliers_SELECTED()
	{
		return coordinates_Outliers_SELECTED;
	}

	public Matrix getAdjustments_per_variable_REMOVE_MAD()
	{
		return adjustments_per_variable_REMOVE_MAD;
	}

	public Matrix getAdjustments_per_variable_REMOVE_Z()
	{
		return adjustments_per_variable_REMOVE_Z;
	}

	public Matrix getAdjustments_per_variable_SELECT_MAD()
	{
		return adjustments_per_variable_SELECT_MAD;
	}

	public Matrix getAdjustments_per_variable_SELECT_Z()
	{
		return adjustments_per_variable_SELECT_Z;
	}

	public int getNumber_of_atoms()
	{
		return number_of_atoms;
	}
}
