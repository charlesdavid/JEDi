package jedi;

import java.math.RoundingMode;
import java.text.NumberFormat;
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

public class JEDi_Get_Transformed_Residue_Coordinates
{
	String res_index;
	final int ROWS, COLS, number_of_atoms;
	List<Double> original_conformation_rmsds, aligned_conformation_rmsds, transformed_residue_rmsd_list, Z_Scores;
	Matrix aligned_subset_REF_PDB_coordinates, aligned_subset_PDB_coordinates, residue_stats, coordinates_Outliers_REMOVED, coordinates_Outliers_SELECTED,
			adjustments_per_variable_REMOVE_Outliers, adjustments_per_variable_SELECT_Outliers;
	final Matrix subset_REF_PDB_coordinates, subset_PDB_coordinates;
	final NumberFormat nf;
	final RoundingMode rm;
	Transform_Coordinates tf_coords;
	Outlier_Processing op;

	/* ************************************** CONSTRUCTORS ******************************************************************************** */

	/**
	 * Constructor for transforming the original coordinates using quaternion operations and calculating the conformation and residue RMSDs. If specified, the coordinates outliers
	 * will be handled This class is used for Cartesian analysis, but not the distance analyses which use internal coordinates (distances)
	 * 
	 * @param coords     The original PDB coordinates
	 * @param ref_coords The reference coordinates to use for the transformation
	 * @param dir        The working directory
	 * @param des        The job description
	 */
	JEDi_Get_Transformed_Residue_Coordinates(Matrix coords, Matrix ref_coords)
	{

		this.nf = NumberFormat.getInstance();
		this.rm = RoundingMode.HALF_UP;
		this.nf.setRoundingMode(rm);
		this.nf.setMaximumFractionDigits(3);
		this.nf.setMinimumFractionDigits(3);

		this.subset_PDB_coordinates = coords;
		this.subset_REF_PDB_coordinates = ref_coords;

		this.ROWS = subset_PDB_coordinates.getRowDimension();
		this.COLS = subset_PDB_coordinates.getColumnDimension();
		this.number_of_atoms = (ROWS / 3);

		transformed_residue_rmsd_list = new ArrayList<>();
		original_conformation_rmsds = new ArrayList<>();
		aligned_conformation_rmsds = new ArrayList<>();
		aligned_subset_PDB_coordinates = new Matrix(ROWS, COLS);
		tf_coords = new Transform_Coordinates(subset_REF_PDB_coordinates, subset_REF_PDB_coordinates);
		aligned_subset_REF_PDB_coordinates = tf_coords.get_transformed_coords();
	}

	/* ************************************** METHODS ******************************************************************************** */

	/**
	 * @return The Transformed coordinates for the specified subset
	 */
	public Matrix get_SS_Transformed_coords()
	{
		for (int i = 0; i < COLS; i++)
			{
				Matrix fc = subset_PDB_coordinates.getMatrix(0, ROWS - 1, i, i);
				tf_coords = new Transform_Coordinates(aligned_subset_REF_PDB_coordinates, fc);
				Matrix transformed_col_vector = tf_coords.get_transformed_coords();
				aligned_subset_PDB_coordinates.setMatrix(0, ROWS - 1, i, i, transformed_col_vector);
				tf_coords = null;
			}
		return aligned_subset_PDB_coordinates;
	}

	/**
	 * @return The transformed conformation RMSDs
	 */
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

	public Matrix get_SS_coordinate_STATS()
	{
		op = new Outlier_Processing(aligned_subset_PDB_coordinates);

		return residue_stats = op.getCoordinateStats();
	}

	public Matrix get_SS_transformed_coordinates_OUTLIERS_REMOVED()
	{
		if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
			{
				coordinates_Outliers_REMOVED = op.remove_outliers_MAD();
				adjustments_per_variable_REMOVE_Outliers = op.getCounts_remove_MAD();
			}
		if (Input_Parameters.Z_SCORE_CUTOFF > 0)
			{
				coordinates_Outliers_REMOVED = op.remove_outliers_Z();
				adjustments_per_variable_REMOVE_Outliers = op.getCounts_remove_Z();
			}

		return coordinates_Outliers_REMOVED;
	}

	public Matrix get_SS_transformed_coordinates_OUTLIERS_SELECTED()
	{
		if (Input_Parameters.MAD_SCORE_CUTOFF > 0)
			{
				coordinates_Outliers_SELECTED = op.select_outliers_MAD();
				adjustments_per_variable_SELECT_Outliers = op.getCounts_select_MAD();
			}
		if (Input_Parameters.Z_SCORE_CUTOFF > 0)
			{
				coordinates_Outliers_SELECTED = op.select_outliers_Z();
				adjustments_per_variable_SELECT_Outliers = op.getCounts_select_Z();
			}

		return coordinates_Outliers_SELECTED;
	}


	/**
	 * @return The Residue RMSDs (RMSFs)
	 */
	public List<Double> get_SS_Residue_RMSFs()
	{
		JEDi_Get_Atomic_RMSFs rrmsd = new JEDi_Get_Atomic_RMSFs(aligned_subset_PDB_coordinates);
		transformed_residue_rmsd_list = rrmsd.get_residue_rmsfs();
		return transformed_residue_rmsd_list;
	}

	/* ************************************** GETTERS ******************************************************************************** */

	public Matrix get_Residue_stats()
	{
		return residue_stats;
	}

	public Matrix get_Transformed_reference_coordinates()
	{
		return aligned_subset_REF_PDB_coordinates;
	}

	public Matrix getCoordinates_Outliers_REMOVED()
	{
		return coordinates_Outliers_REMOVED;
	}

	public Matrix getCoordinates_Outliers_SELECTED()
	{
		return coordinates_Outliers_SELECTED;
	}

	public Matrix getAdjustments_per_variable_REMOVE_Outliers()
	{
		return adjustments_per_variable_REMOVE_Outliers;
	}

	public Matrix getAdjustments_per_variable_SELECT_Outliers()
	{
		return adjustments_per_variable_SELECT_Outliers;
	}
}
