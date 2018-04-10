//package fpgen;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.GraphOnlyFingerprinter;
import org.openscience.cdk.fingerprint.KlekotaRothFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.EStateFingerprinter;
import org.openscience.cdk.fingerprint.HybridizationFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;

import org.openscience.cdk.interfaces.IAtomContainer;

import org.openscience.cdk.smiles.SmilesParser;

public class FingerprintGenerator {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		//int i = 0;
		String compoundList = args[0];
		//String compoundList = "/Users/nsadawi/Downloads/Compound_000000001_000025000.smi.txt";
		List<List<String>> input = readFile(compoundList, false);
		String outputDir = args[1];
		long startTime = System.currentTimeMillis();
		for(List<String> line : input){
			//i++;
			String smiles = line.get(0);
			String ID  = line.get(1);
			System.out.println(smiles + " " + ID);
			//if(i>10) break;
			
			PrintWriter out = new PrintWriter(outputDir+"/"+ID+".fp");

			SmilesParser sp=new SmilesParser(DefaultChemObjectBuilder.getInstance());
			//String smiles = "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2";

			IAtomContainer molecule1;
			try {
				molecule1 = sp.parseSmiles(smiles);

				KlekotaRothFingerprinter krFingerprinter = new KlekotaRothFingerprinter();
				BitSet bs = krFingerprinter.getFingerprint(molecule1);
				printFP(bs, krFingerprinter.getSize(),out);

				MACCSFingerprinter maccsFingerprinter = new MACCSFingerprinter(); 
				bs = maccsFingerprinter.getFingerprint(molecule1);
				printFP(bs, maccsFingerprinter.getSize(),out);

				EStateFingerprinter esFingerprinter = new EStateFingerprinter(); 
				bs = esFingerprinter.getFingerprint(molecule1);
				printFP(bs, esFingerprinter.getSize(),out);

				HybridizationFingerprinter hFingerprinter = new HybridizationFingerprinter();//can have size and depth
				bs = hFingerprinter.getFingerprint(molecule1);
				printFP(bs, hFingerprinter.getSize(),out);

				PubchemFingerprinter pchemFingerprinter = new PubchemFingerprinter(); 
				bs = pchemFingerprinter.getFingerprint(molecule1);
				printFP(bs, pchemFingerprinter.getSize(),out);

				SubstructureFingerprinter subFingerprinter = new SubstructureFingerprinter();//accepts user defined set of fragments.
				bs = subFingerprinter.getFingerprint(molecule1);
				printFP(bs, subFingerprinter.getSize(),out);

				ExtendedFingerprinter eFingerprinter = new ExtendedFingerprinter(2048,4); //can have size and depth
				bs = eFingerprinter.getFingerprint(molecule1);
				printFP(bs, eFingerprinter.getSize(),out);

				GraphOnlyFingerprinter graphFingerprinter = new GraphOnlyFingerprinter(); //can have size and depth
				bs = graphFingerprinter.getFingerprint(molecule1);
				printFP(bs, graphFingerprinter.getSize(),out);

				Fingerprinter fingerprinter = new Fingerprinter(1024,4); //can have size and depth
				bs = fingerprinter.getFingerprint(molecule1);
				printFP(bs, fingerprinter.getSize(),out);



			} catch (InvalidSmilesException e) {
				// TODO Auto-generated catch block
				//e.printStackTrace();
				System.out.println("InvalidSmiles -> Could not generate all FPs for compound with ID: " + ID);
			}
			catch (CDKException e) {
				// TODO Auto-generated catch block
				//e.printStackTrace();
				System.out.println("CDKException -> Could not generate all FPs for compound with ID: " + ID);
			}
			catch (Exception e) {
                                // TODO Auto-generated catch block
                                //e.printStackTrace();
                                System.out.println("Unknown Error -> Could not generate all FPs for compound with ID: " + ID);
            }

			out.close();
			 
		}
		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.println("Time taken: " + elapsedTime +" ms");

	}

	public static void printFP(BitSet bs, int size,PrintWriter out){
		int ssize = size - 1;
		for(int z = 0; z <= ssize; z++){
			if(bs.get(z)==false)
				out.print("0");
			else
				out.print("1");

			if (z != ssize)
				out.print(",");
		}
		out.println();
	}

	/**
	 * read csv file line by line and return as list of list of strings
	 * we can skip 1st line in case it has column names
	 * @param  inFile a file name with full path
	 * @param  skipFstLine true or not
	 * @return  a list of list of strings containing file contents
	 */ 
	public static List<List<String>> readFile(String inFile, boolean skipFstLine){
		List<List<String>> csvList = new ArrayList<List<String>>(); // to store lines after they're split
		try
		{
			BufferedReader infile = new BufferedReader( new FileReader( inFile ) ); // input1.txt

			String line;
			if(skipFstLine)
				infile.readLine(); // skip 1st line - it's the headers - column names!
			//reading in the infile to the ArrayList
			//boolean b = true;			
			while((line = infile.readLine()) != null)
			{
				if (line.isEmpty() || line.trim().equals("") || line.trim().equals("\n")){
					// empty line
				}
				else{
					List<String> csvPieces = splitByCommasNotInQuotes(line);				
					csvList.add(csvPieces);
				}

			}
			infile.close();			
		}
		catch (Exception e)
		{
			StringWriter sw = new StringWriter();
			e.printStackTrace(new PrintWriter(sw));
			System.out.println("EXCEPTION CAUGHT: " + sw.toString() );
			System.exit( 0 );
		}
		return csvList;
	}

	final static private Pattern splitSearchPattern = Pattern.compile("\\s+");

	/**
	 * Splits a csv line into a list of strings, 
	 * avoids splitting if the comma is in double quotes
	 * @param  s a csv text line as a string
	 * @return   A list of strings containing all text between commas in the text line
	 */ 
	private static List<String> splitByCommasNotInQuotes(String s) {
		if (s == null)
			return Collections.emptyList();

		List<String> list = new ArrayList<String>();
		Matcher m = splitSearchPattern.matcher(s);
		int pos = 0;
		boolean quoteMode = false;
		while (m.find())
		{
			String sep = m.group();
			if ("\"".equals(sep))
			{
				quoteMode = !quoteMode;
			}
			else if (!quoteMode && " ".equals(sep))
			{
				int toPos = m.start();
				list.add((s.substring(pos, toPos)).trim());
				pos = m.end();
			}
		}
		if (pos < s.length())
			list.add((s.substring(pos)).trim());
		return list;
	}

}
