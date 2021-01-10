package supportIO;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

public class ZIP_IO
{
	public ZIP_IO()
	{

	}

	public void write_compressed_file(String inPath, String outPath) throws java.io.IOException
	{
		FileInputStream in = new FileInputStream(inPath);

		GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(outPath));

		byte[] buf = new byte[1024];
		int len = 0;
		while ((len = in.read(buf)) > 0)
			{
				out.write(buf, 0, len);
			}
		in.close();
		out.close();
	}

	public void read_Compressed_File(String inPath) throws java.io.IOException
	{
		GZIPInputStream in = new GZIPInputStream(new FileInputStream(inPath));

		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0)
			{
				System.out.write(buf, 0, len);
			}
		in.close();
	}

	/* ********************************************************************************************************************* */

	public List<String> read_Lines_from_BZIP2_File(String inPath) // 'inPath' is the path to the '.bz2' file.
	{
		List<String> lines = new ArrayList<String>();
		String line;
		BufferedReader reader;
		BZip2CompressorInputStream bzip;

		try
			{
				bzip = new BZip2CompressorInputStream(new FileInputStream(inPath));
				reader = new BufferedReader(new InputStreamReader(bzip));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						// System.out.println(line);
					}
				reader.close();
				bzip.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return lines;
	}

	public List<String> read_Lines_from_GZIP_File(String inPath) // 'inPath' is the path to the '.gz' file.
	{
		List<String> lines = new ArrayList<String>();
		String line;
		BufferedReader reader;
		GzipCompressorInputStream gzip;

		try
			{
				gzip = new GzipCompressorInputStream(new FileInputStream(inPath));
				reader = new BufferedReader(new InputStreamReader(gzip));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						// System.out.println(line);
					}
				reader.close();
				gzip.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return lines;
	}

	public List<String> read_Lines_from_TAR_File(String inPath)  // 'inPath' is the path to the '.tar' file.
	{
		List<String> lines = new ArrayList<String>();
		String line;
		TarArchiveInputStream tarInput = null;
		BufferedReader reader = null;
		try
			{
				// Define the input streams based on archive type: tar.bz2, tar.gz, or just .tar
				if (inPath.endsWith("bz2")) tarInput = new TarArchiveInputStream(new BZip2CompressorInputStream(new FileInputStream(inPath)));
				else if (inPath.endsWith("gz")) tarInput = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(inPath)));
				else
					tarInput = new TarArchiveInputStream(new FileInputStream(inPath));

				TarArchiveEntry currentEntry = tarInput.getNextTarEntry();
				// System.out.println("\tProcessing file: " + currentEntry.getName() + "\n");
				if (currentEntry.isDirectory())
					{
						System.err.println(" Warning: The first entry is a DIRECTORY... Trying the second...");
						currentEntry = tarInput.getNextTarEntry();
						// System.out.println("\tProcessing file: " + currentEntry.getName() + "\n");
					}

				reader = new BufferedReader(new InputStreamReader(tarInput));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						// System.out.println(line);
					}
				TarArchiveEntry nextEntry = tarInput.getNextTarEntry();
				if (!(nextEntry == null || nextEntry == currentEntry))
					{
						System.err.println(" Warning: The TAR Archive contains more than one file.");
						System.err.println("\tThe next file is: " + nextEntry.getName());
					}
				reader.close();
				tarInput.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return lines;
	}

	public List<String> read_Lines_from_ZIP_File(String inPath) // 'inPath' is the path to the '.zip' file.
	{
		List<String> lines = new ArrayList<String>();
		String line;
		ZipFile zipFile;
		BufferedReader reader;
		InputStream zip;
		try
			{
				zipFile = new ZipFile(inPath);
				Enumeration<? extends ZipEntry> entries = zipFile.entries();
				ZipEntry entry = entries.nextElement();
				zip = zipFile.getInputStream(entry);
				// System.out.println("\tProcessing file: " + entry.getName() + "\n");
				reader = new BufferedReader(new InputStreamReader(zip, "UTF-8"));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						// System.out.println(line);
					}
				if (entries.hasMoreElements()) System.err.println("\tWarning: ZIP file is NOT a SINGLE FILE ARCHIVE.\n");
				reader.close();
				zip.close();
				zipFile.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return lines;
	}

	/**
	 * This method reads text files compressed with GZIP or BZIP2 and returns their contents as an array of strings.
	 * 
	 * @param inPath
	 * @return List<String> lines: The non-empty lines of the file
	 */
	public List<String> read_Lines_from_Compressed_File(String inPath)  // 'inPath' is the path to the 'bbz2' or '.gz' file.
	{
		List<String> lines = new ArrayList<String>();
		String line;
		BufferedReader reader;
		try
			{
				FileInputStream fin = new FileInputStream(inPath);
				BufferedInputStream bis = new BufferedInputStream(fin);
				CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
				reader = new BufferedReader(new InputStreamReader(input));
				while ((line = reader.readLine()) != null && line.length() >= 1)
					{
						lines.add(line);
						// System.out.println(line);
					}
				reader.close();
				fin.close();
				bis.close();
				input.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return lines;
	}

	/* ********************************************************************************************************************* */

	public List<ZipEntry> read_Entries_In_ZIP_Archive(String inPath) // 'inPath' is the path to the '.zip' archive.
	{
		List<ZipEntry> zip_entries = new ArrayList<ZipEntry>();
		ZipFile zipFile;
		try
			{
				zipFile = new ZipFile(inPath);
				Enumeration<? extends ZipEntry> entries = zipFile.entries();
				ZipEntry entry;
				while (entries.hasMoreElements())
					{
						entry = entries.nextElement();
						zip_entries.add(entry);
						// System.out.println(entry.getName());
					}
				zipFile.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return zip_entries;
	}

	public List<TarArchiveEntry> read_Entries_In_TAR_Archive(String inPath)
	{
		List<TarArchiveEntry> entries = new ArrayList<TarArchiveEntry>();
		TarArchiveInputStream tarInput = null;
		try
			{
				// Define the input streams based on archive type: tar.bz2, tar.gz, or just .tar
				if (inPath.endsWith("bz2")) tarInput = new TarArchiveInputStream(new BZip2CompressorInputStream(new FileInputStream(inPath)));
				else if (inPath.endsWith("gz")) tarInput = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(inPath)));
				else
					tarInput = new TarArchiveInputStream(new FileInputStream(inPath));

				TarArchiveEntry currentEntry;
				while ((currentEntry = tarInput.getNextTarEntry()) != null)
					{
						entries.add(currentEntry);
						// System.out.println(currentEntry.getName());
					}
				tarInput.close();
			}
		catch (Exception e)
			{
				e.printStackTrace();
			}
		return entries;
	}

	public void read_Large_Compressed_File(String inPath)
	{
		final int BUFFER_SIZE = (1024 * 64);
		Path textFile = Paths.get(".... .txt");
		Path gzFile = textFile.resolveSibling(textFile.getFileName().toString() + ".gz");

		try
			{
				OutputStream out = new GZIPOutputStream(Files.newOutputStream(gzFile), BUFFER_SIZE);
				Files.copy(textFile, out);


				InputStream in = new GZIPInputStream(Files.newInputStream(gzFile), BUFFER_SIZE);
				Files.copy(in, textFile);
			}
		catch (IOException e)
			{
				e.printStackTrace();
			}
	}

	public static void driver_Read_Entries()
	{
		ZIP_IO zio = new ZIP_IO();
		String path;

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\TAR_GZ\\1A6N_froda_00000001.pdb.tar.gz";
		System.out.println("\nReading entries in a TAR.GZ ARCHIVE  " + path + "\n");
		List<TarArchiveEntry> entries_tar = zio.read_Entries_In_TAR_Archive(path);
		for (ArchiveEntry z : entries_tar)
			{
				System.out.println(z.getName() + "\n");
			}

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\TAR_BZ2\\1A6N_froda_00000001.pdb.tar.bz2";
		System.out.println("\nReading entries in a TAR.BZ2 ARCHIVE  " + path + "\n");
		entries_tar = zio.read_Entries_In_TAR_Archive(path);
		for (ArchiveEntry z : entries_tar)
			{
				System.out.println(z.getName() + "\n");
			}

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\ZIP_ARCHIVE\\archive.zip";
		System.out.println("\nReading entries in a ZIP ARCHIVE  " + path + "\n");
		List<ZipEntry> entries = zio.read_Entries_In_ZIP_Archive(path);
		for (ZipEntry z : entries)
			{
				System.out.println(z.getName() + "\n");
			}
	}

	public static void driver_Read_Lines()
	{
		ZIP_IO zio = new ZIP_IO();
		String path;

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\ZIP_FILE\\1A6N.pdb.zip";
		System.out.println("\nReading ZIP file  " + path + "\n");
		zio.read_Lines_from_ZIP_File(path);

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\GZIP\\1A6N.pdb.gz";
		System.out.println("\nReading GZIP file  " + path + "\n");
		zio.read_Lines_from_Compressed_File(path);

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\BZIP\\1A6N.pdb.bz2";
		System.out.println("\nReading a BZIP file  " + path + "\n");
		zio.read_Lines_from_Compressed_File(path);
	}

	public static void driver_Read_Lines2()
	{
		ZIP_IO zio = new ZIP_IO();
		String path;

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\TAR_BZ2\\1A6N_froda_00000001.pdb.tar.bz2";
		System.out.println("\nReading a SINGLE FILE TAR.BZ2 ARCHIVE  " + path + "\n");
		zio.read_Lines_from_TAR_File(path);

		path = "C:\\Users\\cflcyd\\eclipse-workspace\\JEDi_Release\\test\\compression\\TAR_GZ\\1A6N_froda_00000001.pdb.tar.gz";
		System.out.println("\nReading a SINGLE FILE TAR.GZ ARCHIVE  " + path + "\n");
		zio.read_Lines_from_TAR_File(path);
	}

	public static void main(String[] args) throws IOException
	{
		// driver_Read_Entries();
		// driver_Read_Lines();
		// driver_Read_Lines2();
	}
}
