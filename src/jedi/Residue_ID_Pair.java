package jedi;

import java.io.Serializable;

/**
 * JED class Residue_ID_Pair: Data structure for pairing a residue Chain ID with its corresponding Residue Number. Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

public class Residue_ID_Pair implements Serializable, Comparable<Residue_ID_Pair>
	{
		private static final long serialVersionUID = 2911455749954226544L;
		String chain_ID;
		Integer residue_Number;

		/**
		 * Constructs a residue ID pair using the specified chain ID and residue number
		 * 
		 * @param chainID The chain ID
		 * @param resNum  The corresponding residue number
		 */
		public Residue_ID_Pair(String chainID, Integer resNum)
			{
				this.chain_ID = chainID;
				this.residue_Number = resNum;
			}

		/**
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Residue_ID_Pair aThat)
			{
				final int EQUAL = 0;
				final int NOT_EQUAL = 1;

				if (this == aThat)
					return EQUAL;

				if (this.chain_ID != aThat.chain_ID)
					return NOT_EQUAL;

				if (this.residue_Number != aThat.residue_Number)
					return NOT_EQUAL;

				return EQUAL;
			}

		/**
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object aThat)
			{
				if (this == aThat)
					return true;

				if (!(aThat instanceof Residue_ID_Pair))
					return false;

				Residue_ID_Pair that = (Residue_ID_Pair) aThat;

				return (this.residue_Number == that.residue_Number) && (this.chain_ID == that.chain_ID);
			}

		/**
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode()
			{
				return residue_Number.hashCode() + chain_ID.hashCode();
			}

		/**
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString()
			{
				return chain_ID + "\t" + residue_Number;
			}
	}
