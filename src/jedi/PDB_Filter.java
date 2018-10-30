package jedi;

import java.io.File;
import java.io.FilenameFilter;

class PDB_Filter implements FilenameFilter
	{
		String filter_String = ".pdb";

		public PDB_Filter(String filter_text)
			{
				this.filter_String = filter_text;
			}

		@Override
		public boolean accept(File dir, String name)
			{
				return (name.endsWith(filter_String));
			}
	}
