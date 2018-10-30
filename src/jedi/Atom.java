package jedi;

/**
 * JED class Atom: The data object for atoms in a PDB file. Copyright (C) 2012 Dr. Charles David
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

public class Atom
	{

		public int res_number, atom_number, charge;
		public String header, symbol, chainID, res_type, code;
		public double x, y, z, occupancy, b_factor;

		/*
		 * **************************** CONSTRUCTORS ****************************************
		 */

		/**
		 * Constructs an empty atom with no data
		 */
		public Atom()
			{

			}

		/**
		 * Constructs an atom with all of its data fields specified.
		 * 
		 * @param res_number
		 * @param atom_number
		 * @param charge
		 * @param header
		 * @param symbol
		 * @param chainID
		 * @param res_type
		 * @param code
		 * @param x
		 * @param y
		 * @param z
		 * @param occupancy
		 * @param b_factor
		 */
		public Atom(int res_number, int atom_number, int charge, String header, String symbol, String chainID, String res_type, String code, double x, double y, double z,
				double occupancy, double b_factor)
			{
				this.res_number = res_number;
				this.atom_number = atom_number;
				this.charge = charge;
				this.header = header;
				this.symbol = symbol;
				this.chainID = chainID;
				this.res_type = res_type;
				this.code = code;
				this.x = x;
				this.y = y;
				this.z = z;
				this.occupancy = occupancy;
				this.b_factor = b_factor;
			}

		/*
		 * ******************************* METHODS *******************************************
		 */

		/**
		 * Provides a meaningful string representation of an atom.
		 * 
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString()
			{
				return String.format("Atom [header=%s, atom_number=%s, symbol=%s, res_type=%s, chainID=%s, res_number=%s, x=%s, y=%s, z=%s, b_factor=%s, code=%s]", header,
						atom_number, symbol, res_type, chainID, res_number, x, y, z, b_factor, code);
			}

		/**
		 * Ensures that the hashcode is unique.
		 * 
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode()
			{
				final int prime = 31;
				int result = 1;
				result = prime * result + atom_number;
				result = prime * result + ((chainID == null) ? 0 : chainID.hashCode());
				result = prime * result + ((header == null) ? 0 : header.hashCode());
				result = prime * result + res_number;
				result = prime * result + ((res_type == null) ? 0 : res_type.hashCode());
				result = prime * result + ((symbol == null) ? 0 : symbol.hashCode());
				long temp;
				temp = Double.doubleToLongBits(x);
				result = prime * result + (int) (temp ^ (temp >>> 32));
				temp = Double.doubleToLongBits(y);
				result = prime * result + (int) (temp ^ (temp >>> 32));
				temp = Double.doubleToLongBits(z);
				result = prime * result + (int) (temp ^ (temp >>> 32));
				return result;
			}

		/**
		 * Provides a way to determine if two atoms are the same.
		 * 
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj)
			{
				if (this == obj)
					return true;
				if (obj == null)
					return false;
				if (getClass() != obj.getClass())
					return false;
				Atom other = (Atom) obj;
				if (atom_number != other.atom_number)
					return false;
				if (chainID == null)
					{
						if (other.chainID != null)
							return false;
					} else if (!chainID.equals(other.chainID))
					return false;
				if (header == null)
					{
						if (other.header != null)
							return false;
					} else if (!header.equals(other.header))
					return false;
				if (res_number != other.res_number)
					return false;
				if (res_type == null)
					{
						if (other.res_type != null)
							return false;
					} else if (!res_type.equals(other.res_type))
					return false;
				if (symbol == null)
					{
						if (other.symbol != null)
							return false;
					} else if (!symbol.equals(other.symbol))
					return false;
				if (Double.doubleToLongBits(x) != Double.doubleToLongBits(other.x))
					return false;
				if (Double.doubleToLongBits(y) != Double.doubleToLongBits(other.y))
					return false;
				if (Double.doubleToLongBits(z) != Double.doubleToLongBits(other.z))
					return false;
				return true;
			}

		/*
		 * ******************************* SETTERS *******************************************
		 */
		/**
		 * @param res_number the res_number to set
		 */
		public void setRes_number(int res_number)
			{
				this.res_number = res_number;
			}

		/**
		 * @param atom_number the atom_number to set
		 */
		public void setAtom_number(int atom_number)
			{
				this.atom_number = atom_number;
			}

		/**
		 * @param charge the charge to set
		 */
		public void setCharge(int charge)
			{
				this.charge = charge;
			}

		/**
		 * @param header the header to set
		 */
		public void setHeader(String header)
			{
				this.header = header;
			}

		/**
		 * @param symbol the symbol to set
		 */
		public void setSymbol(String symbol)
			{
				this.symbol = symbol;
			}

		/**
		 * @param chainID the chainID to set
		 */
		public void setChainID(String chainID)
			{
				this.chainID = chainID;
			}

		/**
		 * @param res_type the res_type to set
		 */
		public void setRes_type(String res_type)
			{
				this.res_type = res_type;
			}

		/**
		 * @param code the code to set
		 */
		public void setCode(String code)
			{
				this.code = code;
			}

		/**
		 * @param x the x to set
		 */
		public void setX(double x)
			{
				this.x = x;
			}

		/**
		 * @param y the y to set
		 */
		public void setY(double y)
			{
				this.y = y;
			}

		/**
		 * @param z the z to set
		 */
		public void setZ(double z)
			{
				this.z = z;
			}

		/**
		 * @param occupancy the occupancy to set
		 */
		public void setOccupancy(double occupancy)
			{
				this.occupancy = occupancy;
			}

		/**
		 * @param b_factor the b_factor to set
		 */
		public void setB_factor(double b_factor)
			{
				this.b_factor = b_factor;
			}

		/*
		 * ******************************* GETTERS *******************************************
		 */
		/**
		 * @return the res_number
		 */
		public int getRes_number()
			{
				return res_number;
			}

		/**
		 * @return the atom_number
		 */
		public int getAtom_number()
			{
				return atom_number;
			}

		/**
		 * @return the charge
		 */
		public int getCharge()
			{
				return charge;
			}

		/**
		 * @return the header
		 */
		public String getHeader()
			{
				return header;
			}

		/**
		 * @return the symbol
		 */
		public String getSymbol()
			{
				return symbol;
			}

		/**
		 * @return the chainID
		 */
		public String getChainID()
			{
				return chainID;
			}

		/**
		 * @return the res_type
		 */
		public String getRes_type()
			{
				return res_type;
			}

		/**
		 * @return the code
		 */
		public String getCode()
			{
				return code;
			}

		/**
		 * @return the x
		 */
		public double getX()
			{
				return x;
			}

		/**
		 * @return the y
		 */
		public double getY()
			{
				return y;
			}

		/**
		 * @return the z
		 */
		public double getZ()
			{
				return z;
			}

		/**
		 * @return the occupancy
		 */
		public double getOccupancy()
			{
				return occupancy;
			}

		/**
		 * @return the b_factor
		 */
		public double getB_factor()
			{
				return b_factor;
			}
	}
