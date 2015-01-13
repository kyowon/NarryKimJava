package msgfpipeline.parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import suffixarray.MatchSet;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.Peptide;
import net.sf.samtools.util.BufferedLineReader;

public class MSGFPlusParser {
	
	public static class MSGFPlusPSM implements Comparable<MSGFPlusPSM>{
	//	private SuffixArray sa; 
		static private Enzyme enzyme;
		static private boolean listAllProteins = true;
		

		private String specFile;
		private String specID;
		private int scanNumber;
		private String fragMethod;
		private float precursor;
		private int isotopeErr;
		private float precursorErr;
		private int charge;
		private String peptide;
		//private String[] proteins;
		private int[] startInPro;
		private int[] endInPro;
		private int[] proteinLength;
		
		private String[] proteinID = null;
		private String[] proteinName = null;
		private String[] species = null;
		private String[] geneName = null;
		private float[] proteinMass = null;
		
		private String[] preAAs;
		private String[] postAAs;
		private int deNovoScore;
		private int msgfScore;
		private float specEvalue;
		private float evalue;
		private float qvalue = -1; // optional
		private float pepQvalue = -1; // optional
	//	private String misc = "";
		
		private float recalibratedPrecursor;
		private float recalibratedPrecursorErr;
		
		public static void setListAllProteins(boolean s){
			listAllProteins = s;
		}
		
		public void setQvalue(float q){
			this.qvalue  = q;
		}
		
		public void setPepQvalue(float q){
			this.pepQvalue = q;
		}
		
		public void setRecalibratedPrecursor(float p){
			this.recalibratedPrecursor = p;
			this.recalibratedPrecursorErr = precursorErr - (precursor - p);
		}
		
		public float getRecalibratedPrecursorErr(){
			return recalibratedPrecursorErr;
		}		

		public int[] getStartInPro() {
			return startInPro;
		}

		public int[] getEndInPro() {
			return endInPro;
		}

		public int[] getProteinLength() {
			return proteinLength;
		}

		public String[] getProteinIDs() {
			return proteinID;
		}

		public String[] getProteinNames() {
			return proteinName;
		}

		public String[] getGeneNames() {
			return geneName;
		}

		public String[] getSpecies() {
			return species;
		}
		
		public float[] getProteinMass() {
			return proteinMass;
		}
		
		public String getSpecFile() {
			return specFile;
		}
		
		public String getSpecID() {
			return specID;
		}

		public int getScanNumber() {
			return scanNumber;
		}

		public String getFragMethod() {
			return fragMethod;
		}

		public float getPrecursor() {
			return precursor;
		}

		public int getIsotopeErr() {
			return isotopeErr;
		}

		public float getPrecursorErr() {
			return precursorErr;
		}

		public int getCharge() {
			return charge;
		}

		public String getPeptide() {
			return peptide;
		}

		public int getDeNovoScore() {
			return deNovoScore;
		}

		public int getMsgfScore() {
			return msgfScore;
		}

		public float getSpecEvalue() {
			return specEvalue;
		}

		public float getEvalue() {
			return evalue;
		}

		public float getQvalue() {
			return qvalue;
		}

		public float getPepQvalue() {
			return pepQvalue;
		}
		public String[] getPreAAs(){
			return preAAs;
		}
		public String[] getPostAAs(){
			return postAAs;
		}
		
		
		
		public MSGFPlusPSM(String s, SuffixArray sa){
			//this.sa = sa;
			//this.enzyme = enzyme;
			String[] token = s.split("\t");
			specFile = token[0];
			specID = token[1];
			scanNumber = Integer.parseInt(token[2]);
			fragMethod = token[3];
			precursor = Float.parseFloat(token[4]);
			isotopeErr = Integer.parseInt(token[5]);
			precursorErr = Float.parseFloat(token[6]);
			charge = Integer.parseInt(token[7]);
			peptide = token[8];
			
			deNovoScore = Integer.parseInt(token[10]);
			msgfScore = Integer.parseInt(token[11]);
			specEvalue = Float.parseFloat(token[12]);
			evalue = Float.parseFloat(token[13]);
			if(token.length > 14){
				qvalue = Float.parseFloat(token[14]);
				pepQvalue = Float.parseFloat(token[15]);
			}
			
			SuffixArraySequence sequence = sa.getSequence();
			MatchSet matchSet = sa.findAll(Peptide.getUnmodifiedSequence(peptide).toString());
			int numProtein = matchSet.getSize();
						
			preAAs = new String[numProtein];
			postAAs = new String[numProtein];
			startInPro = new int[numProtein];
			endInPro = new int[numProtein];
			proteinID = new String[numProtein];
			proteinName = new String[numProtein];
			geneName = new String[numProtein];
			species = new String[numProtein];
			proteinMass = new float[numProtein];
			proteinLength = new int[numProtein];
			String[] proteinSeq = new String[numProtein];
			
			for(int i=0;i<numProtein;i++){
				int start = matchSet.getStart(i);
				int end = matchSet.getEnd(i);
				preAAs[i] = sequence.getSubsequence(Math.max(0, start-1), start);
				postAAs[i] = sequence.getSubsequence(end, Math.min(end+1, sequence.getSize())); // bug fixed by kyowon
				int proStart = (int) sequence.getStartPosition(start);
				startInPro[i] = start -proStart;
				endInPro[i] = end - proStart;
				
				String annotation = sequence.getAnnotation(start).split(";")[0];
				int spaceIndex = annotation.indexOf(' ');
				if(spaceIndex < 0) spaceIndex = annotation.length();
				proteinID[i] = annotation.substring(0, spaceIndex);
				
				int pIndex = -1;
				int gnIndex = annotation.indexOf("GN=");
				if(gnIndex > 0){
					geneName[i] = annotation.substring(gnIndex + 3);
					int j = geneName[i].indexOf(' ');
					geneName[i] = j>0? geneName[i].substring(0, j) : geneName[i];
					pIndex = gnIndex;
					//proteinName[i] = annotation.substring(spaceIndex, gnIndex);
				}else{
					//proteinName[i] = annotation.substring(spaceIndex);
					geneName[i] = "N/A";
				}
				
				int osIndex = annotation.indexOf("OS=");
				if(osIndex > 0){
					species[i] = annotation.substring(osIndex + 3);
					int j = species[i].indexOf(' ');
					species[i] = j>0? species[i].substring(0, j) : species[i];
					pIndex = pIndex<0? osIndex : Math.min(osIndex, pIndex);
					//proteinName[i] = annotation.substring(spaceIndex, osIndex);
				}else{
					//proteinName[i] = annotation.substring(spaceIndex);
					species[i] = "N/A";
				}
				
				if(pIndex > 0) proteinName[i] = annotation.substring(spaceIndex, pIndex);
				else proteinName[i] = annotation.substring(spaceIndex);
				
				proteinSeq[i] = sequence.getMatchingEntry(start);
				proteinLength[i] = proteinSeq[i].length();
				proteinMass[i] = new Peptide(proteinSeq[i]).getMass();
			}
		}
		
		
		
		@Override
		public String toString(){
		
			StringBuffer sb = new StringBuffer();
			sb.append(specFile);sb.append('\t');
			sb.append(specID);sb.append('\t');
			sb.append(scanNumber);sb.append('\t');
			sb.append(fragMethod);sb.append('\t');
			sb.append(precursor);sb.append('\t');			
			sb.append(recalibratedPrecursor);sb.append('\t');
			sb.append(new Peptide(peptide).getMass() + Composition.H2O + Composition.PROTON);sb.append('\t');
			sb.append(isotopeErr);sb.append('\t');
			sb.append(precursorErr);sb.append('\t');
			sb.append(recalibratedPrecursorErr);sb.append('\t');
			sb.append(charge);sb.append('\t');
			
			int n = proteinID.length;
			if(listAllProteins) n = Math.min(1, n);
			
			for(int i=0;i<n;i++){
				sb.append(preAAs[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			sb.append(peptide);sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(postAAs[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(startInPro[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(endInPro[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			if(enzyme != null){
				for(int i=0;i<n;i++){
					sb.append(enzyme.getNumCleavedTermini(preAAs[i]+ "." + Peptide.getUnmodifiedSequence(peptide) + "." + postAAs[i], AminoAcidSet.getStandardAminoAcidSet()));
					if(i<n-1) sb.append(';');
				}
				sb.append('\t');
				sb.append(new Peptide(peptide).getNumMissedCleavageSites(enzyme));sb.append('\t');
			}else{
				sb.append("N/A");sb.append('\t');
				sb.append("N/A");sb.append('\t');
			}
			for(int i=0;i<n;i++){
				sb.append(proteinID[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(proteinName[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(geneName[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(species[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			sb.append(proteinID.length);sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(proteinMass[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			for(int i=0;i<n;i++){
				sb.append(proteinLength[i]);
				if(i<n-1) sb.append(';');
			}
			sb.append('\t');
			sb.append(deNovoScore);sb.append('\t');
			sb.append(msgfScore);sb.append('\t');
			sb.append(specEvalue);sb.append('\t');
			sb.append(evalue);
			sb.append('\t');sb.append(qvalue);
			sb.append('\t');sb.append(pepQvalue);
			
					
			return sb.toString();
		}

		public int compareTo(MSGFPlusPSM o) {
			return new Integer(this.getScanNumber()).compareTo(o.getScanNumber());
		}
	}
	
	private ArrayList<MSGFPlusPSM> psms;
	public MSGFPlusParser(String fileName, String fasta, Enzyme enzyme){
		psms = new ArrayList<MSGFPlusPSM>();
		MSGFPlusPSM.enzyme = enzyme;
		SuffixArraySequence sequence = new SuffixArraySequence(fasta);
		SuffixArray sa = new SuffixArray(sequence);
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(fileName));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("#") || s.startsWith("ID")){
					continue;
				}
				psms.add(new MSGFPlusPSM(s, sa));
			}in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}	
	
	public ArrayList<MSGFPlusPSM> getPsms() {
		return psms;
	}
	
	public static String getHeader(){
		String s= "#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tRecalibratedPrecursor\tMH+(theor)\tIsotopeError\t"
				+ "PrecursorError(Da)\tRecalibratedPrecursorError(Da)\tCharge\tPre\tPeptide\tPost\tStartPosition\tEndPosition\tNTT\tNumMC\tProteinID\tProteinName\tGeneName\tSpecies\tMultiProt\t"
				+ "ProtMw\tProtLength\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\tQValue\tPepQValue";
		return s;			
	}
	
	public static int getSpecEvalueCol(){
		String[] token = getHeader().split("\t");
		for(int i=0;i<token.length;i++){
			if(token[i].equals("SpecEValue")) return i;
		}
		return - 1; 
	}
	
	public static int getScanNumberCol(){
		String[] token = getHeader().split("\t");
		for(int i=0;i<token.length;i++){
			if(token[i].equals("ScanNum")) return i;
		}
		return - 1; 
	}
	
	public static int getPeptideCol(){
		String[] token = getHeader().split("\t");
		for(int i=0;i<token.length;i++){
			if(token[i].equals("Peptide")) return i;
		}
		return - 1; 
	}
	
	public static int getProteinIDCol(){
		String[] token = getHeader().split("\t");
		for(int i=0;i<token.length;i++){
			if(token[i].equals("ProteinID")) return i;
		}
		return - 1; 
	}
/*	static public void main(String[] args){
		
		if(args.length < 3){
			//System.out.println(Enzyme.TRYPSIN.getResidues().get(0).getResidue() + " " + Enzyme.TRYPSIN.isCTerm());
			
			System.out.println("Usage : \n"
					+ "java -jar -Xmx1g Reformatter.jar [msgf tsv file name] [reformatted output file name] [fasta file name] [enzyme option]\n"
					+ "\n Enzyme : 0: unspecific cleavage, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Arg-C, 6: Asp-N");
			return;
		}
		
		
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		if(args.length >= 4){
			int ei = Integer.parseInt(args[3]);
			if(ei == 0) enzyme = null;
			else if(ei == 2) enzyme = Enzyme.CHYMOTRYPSIN;
			else if(ei == 3) enzyme = Enzyme.LysC;
			else if(ei == 4) enzyme = Enzyme.LysN;
			else if(ei == 5) enzyme = Enzyme.ArgC;
			else if(ei == 6) enzyme = Enzyme.AspN;
		}
		//[-e EnzymeID] (0: unspecific cleavage, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N,
				//5: : Arg-C, 6: Asp-N,)
		MSGFPlusParser parser = new MSGFPlusParser(args[0]);
		try {
			PrintStream out = new PrintStream(args[1]);
			SuffixArraySequence sequence = new SuffixArraySequence(args[2]);
			SuffixArray sa = new SuffixArray(sequence);
			boolean start = true;
			for(MSGFPlusPSM psm : parser.getPsms()){
				if(start){
					start = false;
					out.println(psm.toFormatChangedHeader());
				}
				out.println(psm.toFormatChangedString(sa, enzyme));
			}
			
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}*/
}
