package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
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
	private String header = "";
	
	public class MSGFPlusPSM{
		private String specFile;
		private String specID;
		private int scanNumber;
		private String fragMethod;
		private float precursor;
		private int isotopeErr;
		private float precursorErr;
		private int charge;
		private String peptide;
		private String[] proteins;
		private char[] preAAs;
		private char[] postAAs;
		private int deNovoScore;
		private int msgfScore;
		private float specEvalue;
		private float evalue;
		private float qvalue = -1; // optional
		private float pepQvalue = -1; // optional
		private String misc = "";
		
		public String getSpecFile() {
			return specFile;
		}
		
		public String getMiscString(){
			return misc;
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

		public String[] getProteins() {
			return proteins;
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
		public char[] getPreAAs(){
			return preAAs;
		}
		public char[] getPostAAs(){
			return postAAs;
		}
		
		@Override
		public String toString(){
			StringBuffer ret = new StringBuffer();
			ret.append(specFile);ret.append('\t');
			ret.append(specID);ret.append('\t');
			ret.append(scanNumber);ret.append('\t');
			ret.append(fragMethod);ret.append('\t');
			ret.append(precursor);ret.append('\t');
			ret.append(isotopeErr);ret.append('\t');
			ret.append(precursorErr);ret.append('\t');
			ret.append(charge);ret.append('\t');
			ret.append(peptide);ret.append('\t');
			for(int i=0; i<proteins.length;i++){
				ret.append(proteins[i]);ret.append(i<proteins.length-1 ? ";" : "\t");
			}
			
			ret.append(deNovoScore);ret.append('\t');
			ret.append(msgfScore);ret.append('\t');
			ret.append(specEvalue);ret.append('\t');
			ret.append(evalue);
			if(qvalue > -1){
				ret.append('\t'); ret.append(qvalue);
				ret.append('\t'); ret.append(pepQvalue);
			}
			ret.append('\t');
			ret.append(misc);
			return ret.toString();
		}	
		
		public MSGFPlusPSM(String s){
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
			proteins = token[9].split(";");
			deNovoScore = Integer.parseInt(token[10]);
			msgfScore = Integer.parseInt(token[11]);
			specEvalue = Float.parseFloat(token[12]);
			evalue = Float.parseFloat(token[13]);
			if(token.length > 14){
				qvalue = Float.parseFloat(token[14]);
				pepQvalue = Float.parseFloat(token[15]);
			}
			
			for(int i=16;i<token.length;i++){
				misc += token[i] + (i == token.length-1 ? "" : "\t");
			}
			
			preAAs = new char[proteins.length];
			postAAs = new char[proteins.length];
			for(int i=0;i<proteins.length;i++){
				preAAs[i] = proteins[i].charAt(proteins[i].lastIndexOf("(pre=")+5);
				postAAs[i] = proteins[i].charAt(proteins[i].lastIndexOf("post=")+5);
			}
		}
		
		public String toFormatChangedHeader(){
			String s= "#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tMH+(theor)\tIsotopeError\t"
					+ "PrecursorError(Da)\tCharge\tPre\tPeptide\tPost\tStartAA\tEndAA\tNTT\t#MC\tProteinID\tProteinName\tMultiProt\t"
					+ "ProtMw\tProtLeng\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue";
			
			if(qvalue >= 0){
				s +="\tQValue\tPepQVale";
			}
			if(!misc.isEmpty()){
				String[] token = misc.split("\t");
				for(int i=0;i<token.length-1;i++){
					s += "misc"+i+"\t";
				}
				s += "misc" + (token.length-1);
			}
			
			return s;			
		}
		
		public String toFormatChangedString(SuffixArray sa, Enzyme enzyme){
			StringBuffer sb = new StringBuffer();
			sb.append(specFile);sb.append('\t');
			sb.append(specID);sb.append('\t');
			sb.append(scanNumber);sb.append('\t');
			sb.append(fragMethod);sb.append('\t');
			sb.append(precursor);sb.append('\t');			
			sb.append(new Peptide(peptide).getMass() + Composition.H2O + Composition.PROTON);sb.append('\t');
			sb.append(isotopeErr);sb.append('\t');
			sb.append(precursorErr);sb.append('\t');
			sb.append(charge);sb.append('\t');
			
			SuffixArraySequence sequence = sa.getSequence();
			MatchSet matchSet = sa.findAll(Peptide.getUnmodifiedSequence(peptide).toString());
			
			int start = matchSet.getStart(0);
			int end = matchSet.getEnd(0);
			String leftStr = sequence.getSubsequence(Math.max(0, start-1), start);
			String rightStr = sequence.getSubsequence(end, Math.min(end+1, sequence.getSize())); // bug fixed by kyowon
			
			int proStart = (int) sequence.getStartPosition(start);
			int startInPro = start -proStart;
			int endInPro = end - proStart;
			
			String annotation = sequence.getAnnotation(start).split(";")[0];
			int spaceIndex = annotation.indexOf(' ');
			if(spaceIndex < 0) spaceIndex = annotation.length();
			String proteinID = annotation.substring(0, spaceIndex);
			String proteinName = annotation.substring(spaceIndex);
			//System.out.println();
			String proteinSeq = sequence.getMatchingEntry(start);
			float proteinMass = new Peptide(proteinSeq).getMass();
				//System.out.println(startInPro + " " + endInPro + " " + proteinID + " " + proteinName + " " + proteinMass + " " + proteinSeq.length() + " " + matchSet.getSize());	
			
			sb.append(leftStr);sb.append('\t');
			sb.append(peptide);sb.append('\t');
			sb.append(rightStr);sb.append('\t');
			sb.append(startInPro);sb.append('\t');
			sb.append(endInPro);sb.append('\t');
		//	System.out.println(leftStr+ "." + peptide + "." + rightStr);
			if(enzyme != null){
				sb.append(enzyme.getNumCleavedTermini(leftStr+ "." + Peptide.getUnmodifiedSequence(peptide) + "." + rightStr, AminoAcidSet.getStandardAminoAcidSet())); sb.append('\t');
				sb.append(new Peptide(peptide).getNumMissedCleavageSites(enzyme));sb.append('\t');
			}else{
				sb.append("N/A");sb.append('\t');
				sb.append("N/A");sb.append('\t');
			}
			
			sb.append(proteinID);sb.append('\t');
			sb.append(proteinName);sb.append('\t');
			sb.append(matchSet.getSize());sb.append('\t');
			sb.append(proteinMass);sb.append('\t');
			sb.append(proteinSeq.length());sb.append('\t');
			sb.append(deNovoScore);sb.append('\t');
			sb.append(msgfScore);sb.append('\t');
			sb.append(specEvalue);sb.append('\t');
			sb.append(evalue);
			
			if(qvalue >= 0){
				sb.append('\t');sb.append(qvalue);
				sb.append('\t');sb.append(pepQvalue);
			}
			if(!misc.isEmpty()){
				sb.append('\t');sb.append(misc);
			}
			
			return sb.toString();
		}
	
		
	
	}
	
	private ArrayList<MSGFPlusPSM> psms;
	public MSGFPlusParser(String fileName){
		psms = new ArrayList<MSGFPlusPSM>();
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(fileName));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")){
					header = s;
					continue;
				}
				psms.add(new MSGFPlusPSM(s));
			}in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}	
	
	public ArrayList<MSGFPlusPSM> getPsms() {
		return psms;
	}
	public String getHeader(){ return header;}
	
	static public void main(String[] args){
		
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
		
		
		
		/*SuffixArraySequence sequence = new SuffixArraySequence("/media/kyowon/Data1/MassSpec/test.fasta");
		SuffixArray sa = new SuffixArray(sequence);
		System.out.println(Peptide.getUnmodifiedSequence("P+20EPTIDE"));
		MatchSet matchSet = sa.findAll("PEPTIDE");
		//sa.getAllMatchingEntries("PEPTIDE");
		//ArrayList<String> matches = new ArrayList<String>();
		
		System.out.println(new Peptide("KPEPKTI+20DEK").getNumMissedCleavageSites(Enzyme.LysN));
		for(int i=0; i<matchSet.getSize(); i++)
		{
			int start = matchSet.getStart(i);
			int end = matchSet.getEnd(i);
			String leftStr = sequence.getSubsequence(Math.max(0, start-1), start);
			String rightStr = sequence.getSubsequence(end, Math.min(end+1, sequence.getSize())); // bug fixed by kyowon
			
			int proStart = (int) sequence.getStartPosition(start);
			int startInPro = start -proStart;
			int endInPro = end - proStart;
			
			String annotation = sequence.getAnnotation(start).split(";")[0];
			int spaceIndex = annotation.indexOf(' ');
			String proteinID = annotation.substring(0, spaceIndex);
			String proteinName = annotation.substring(spaceIndex);
			//System.out.println();
			String proteinSeq = sequence.getMatchingEntry(start);
			float proteinMass = new Peptide(proteinSeq).getMass();
			System.out.println(startInPro + " " + endInPro + " " + proteinID + " " + proteinName + " " + proteinMass + " " + proteinSeq.length() + " " + matchSet.getSize());
			
			
			
			//sequence.getMatchingEntry(matchSet.getStart(i))
			//System.out.println(sequence.getMatchingE);
			//matches.add(leftStr+"."+sequence.getSubsequence(start, end)+"."+rightStr);
			
		}*/
		//System.out.println(matches);
	}
}
