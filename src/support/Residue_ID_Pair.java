package support;

import java.io.Serializable;

/**
 * JEDi class Residue_ID_Pair: Data structure for pairing a residue Chain ID with its corresponding Residue Number. Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Residue_ID_Pair implements Serializable, Comparable<Residue_ID_Pair>
{
	private static final long serialVersionUID = 2911455749954226544L;
	String chain_ID;
	int residue_Number;

	/**
	 * Constructs a residue ID pair using the specified chain ID and residue number
	 * 
	 * @param chainID
	 *            The chain ID
	 * @param resNum
	 *            The corresponding residue number
	 */
	public Residue_ID_Pair(String chainID, int resNum)
	{
		this.chain_ID = chainID;
		this.residue_Number = resNum;
	}

	@Override
	public int compareTo(Residue_ID_Pair aThat)
	{
		final int EQUAL = 0;
		final int NOT_EQUAL = 1;

		if (this == aThat) return EQUAL;

		if (this.chain_ID != aThat.chain_ID) return NOT_EQUAL;

		if (this.residue_Number != aThat.residue_Number) return NOT_EQUAL;

		return EQUAL;
	}

	@Override
	public boolean equals(Object o)
	{
		return o instanceof Residue_ID_Pair && this.getChain_ID().equals(((Residue_ID_Pair) o).getChain_ID())
				&& this.getResidue_Number().equals(((Residue_ID_Pair) o).getResidue_Number());
	}

	@Override
	public int hashCode()
	{
		Integer val = residue_Number;
		return val.hashCode() + chain_ID.hashCode();
	}

	@Override
	public String toString()
	{
		return String.format("%-10s%-3s%-15s%-3s", "ChainID = ", chain_ID, "Residue Number = ", residue_Number);
		// return String.format("Residue_ID_Pair: " + chain_ID + "\t" + residue_Number);
	}

	public String toShortString()
	{
		return String.format("%-5s%-10s", chain_ID, residue_Number);
		// return String.format("Residue_ID_Pair: " + chain_ID + "\t" + residue_Number);
	}

	// ******************************** GETTERS ******************************************************* //

	public String getChain_ID()
	{
		return chain_ID;
	}

	public Integer getResidue_Number()
	{
		return residue_Number;
	}

	// ******************************** SETTERS ******************************************************* //

	public void setChain_ID(String chain_ID)
	{
		this.chain_ID = chain_ID;
	}

	public void setResidue_Number(Integer residue_Number)
	{
		this.residue_Number = residue_Number;
	}
}
