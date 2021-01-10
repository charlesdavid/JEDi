package support;

import java.io.Serializable;

/**
 * JEDi class Atom_ID_Pair: Data structure for pairing a atom Chain ID with its corresponding atom Number.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * @author Dr. Charles David
 */

public class Atom_ID_Pair implements Serializable, Comparable<Atom_ID_Pair>
{
	private static final long serialVersionUID = 6216348688379341181L;
	String chain_ID;
	Integer atom_Number;

	/**
	 * Constructs a atom ID pair using the specified chain ID and residue number
	 * 
	 * @param chainID
	 *            The chain ID
	 * @param atmNum
	 *            The corresponding atom number
	 */
	public Atom_ID_Pair(String chainID, Integer atmNum)
	{
		this.chain_ID = chainID;
		this.atom_Number = atmNum;
	}

	@Override
	public int compareTo(Atom_ID_Pair aThat)
	{
		final int EQUAL = 0;
		final int NOT_EQUAL = 1;

		if (this == aThat) return EQUAL;

		if (this.chain_ID != aThat.chain_ID) return NOT_EQUAL;

		if (this.atom_Number != aThat.atom_Number) return NOT_EQUAL;

		return EQUAL;
	}

	@Override
	public boolean equals(Object o)
	{
		return o instanceof Atom_ID_Pair && this.getChain_ID().equals(((Atom_ID_Pair) o).getChain_ID()) && this.getAtom_Number().equals(((Atom_ID_Pair) o).getAtom_Number());
	}

	@Override
	public int hashCode()
	{
		return atom_Number.hashCode() + chain_ID.hashCode();
	}

	@Override
	public String toString()
	{
		return String.format("%-9s%-3s%-14s%-3s", "Chain = ", chain_ID, "Atom Number = ", atom_Number);
		// return String.format("Atom_ID_Pair: [chain_ID=%s, atom_Number=%s]", chain_ID, atom_Number);
	}

	// ******************************** GETTERS ******************************************************* //

	public String getChain_ID()
	{
		return chain_ID;
	}

	public Integer getAtom_Number()
	{
		return atom_Number;
	}

	// ******************************** SETTERS ******************************************************* //

	public void setChain_ID(String chain_ID)
	{
		this.chain_ID = chain_ID;
	}

	public void setAtom_Number(Integer atom_Number)
	{
		this.atom_Number = atom_Number;
	}
}
