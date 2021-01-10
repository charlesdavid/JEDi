package supportIO;

import java.io.File;
import java.io.FilenameFilter;

public class PDB_Filter implements FilenameFilter
{
	boolean BZIP, GZIP, ZIP, TAR;
	String filter_String = ".pdb";

	public PDB_Filter(String filter_text)
	{
		this.filter_String = filter_text;
	}

	@Override
	public boolean accept(File dir, String name)
	{
		if (name.endsWith("bz2")) BZIP = true;
		if (name.endsWith("gz")) GZIP = true;
		if (name.endsWith("zip")) ZIP = true;
		if (name.endsWith("tar.bz2") || name.endsWith("tar.gz")) TAR = true;

		if (BZIP) return (name.contains(filter_String));
		else if (GZIP) return (name.contains(filter_String));
		else if (ZIP) return (name.contains(filter_String));
		else if (TAR) return (name.contains(filter_String));
		else
			return (name.endsWith(filter_String));
	}
}
