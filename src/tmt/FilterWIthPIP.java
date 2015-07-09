package tmt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Peptide;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class FilterWIthPIP {
	
	public static class TmtPSM implements Comparable<TmtPSM>{
		private String str;
		private Composition composition;
		private double score;
		private String peptide;
		private int scanNumber;
		private int tmtTagCount;
		private int charge;
		private float mz;
		
		public int getCharge() {
			return charge;
		}

		public float getMz(){
			return mz;
		}
		
		public void setCharge(int charge) {
			this.charge = charge;
		}

		static private String header;
		
		public static void setHeader(String header) {
			TmtPSM.header = header;
		}

		static public String getHeader(){
			return header + "\tPIP Score"; 
		}

		public String toString(){
			return str + "\t" + score;
		}
		
		public String getReadStr() {
			return str;
		}

		public void setReadStr(String str) {
			this.str = str;
		}

		public String getPeptide() {
			return peptide;
		}

		public void setPeptide(String peptide) {
			this.peptide = peptide;
		}

		public int getScanNumber() {
			return scanNumber;
		}

		public void setScanNumber(int scanNumber) {
			this.scanNumber = scanNumber;
		}

		public int getTmtTagCount() {
			return tmtTagCount;
		}

		public void setTmtTagCount(int tmtTagCount) {
			this.tmtTagCount = tmtTagCount;
		}
		
		// constructors...
		
		public TmtPSM(String str){
			this.str = str;
			String[] token = str.split("\t");
			this.scanNumber = Integer.parseInt(token[2]);
			this.charge = Integer.parseInt(token[7]);
			Pattern pattern = Pattern.compile("\\+229.163");
			Matcher  matcher = pattern.matcher(token[8]);
			this.tmtTagCount = 0;
			while (matcher.find())
				tmtTagCount++;
			this.peptide = token[8];//.replaceAll("\\+229.163", "");
			this.mz = (float)( new Peptide(this.peptide, AminoAcidSet.getStandardAminoAcidSet()).getParentMass()/this.charge + Composition.PROTON);
			//float mz = (float)(annotation.getPeptide().getParentMass()/c + Composition.PROTON);
			
			//double pmz = Double.parseDouble(token[4]);
			//double perror = Math.abs(Double.parseDouble(token[6]));
			//if(perror < 5 && Math.abs(pmz - mz) >  2)
			//	System.out.println(scanNumber + " " +  token[8] + "  " + token[4] + " "  + mz);
			//System.out.println(new Peptide(this.peptide));			
			this.composition = new Peptide(new Peptide(this.peptide).toString().toUpperCase()).getComposition(); // TODO sucks but I do not care..
		}
	
	
		TmtPSM(int scanNumber){
			this.str = new Integer(scanNumber).toString();
			this.scanNumber = scanNumber;
		}
		
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((str == null) ? 0 : str.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			TmtPSM other = (TmtPSM) obj;
			if (str == null) {
				if (other.str != null)
					return false;
			} else if (!str.equals(other.str))
				return false;
			return true;
		}

		public int compareTo(TmtPSM o) {
			return new Integer(this.scanNumber).compareTo(o.scanNumber);
		}

		public Composition getComposition() {
			return composition;
		}

		public void setComposition(Composition composition) {
			this.composition = composition;
		}

		public double getScore() {
			return score;
		}

		public void setScore(double score) {
			this.score = score;
		}
	}	
	
	private ArrayList<TmtPSM> psms; // sorted
	private String spectrumFileName;
	private Tolerance tolerance;
	private String outFileName;
	
	public FilterWIthPIP(String outFileName, String inFileName, String spectrumFileName, Tolerance tolerance){
		this.outFileName = outFileName;
		this.spectrumFileName = spectrumFileName;
		this.tolerance = tolerance;
		read(inFileName);
	}
	
	
	private void read(String inFileName){
		psms = new ArrayList<FilterWIthPIP.TmtPSM>();
		try {
			BufferedLineReader in = new BufferedLineReader(inFileName);
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")){
					TmtPSM.setHeader(s);
					continue;
				}
				psms.add(new TmtPSM(s));
			}
			in.close();
			Collections.sort(psms);
			//System.out.println(psms.size());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private ArrayList<TmtPSM> recruitTmtPsms(Spectrum ms1Spectrum, Spectrum nextMs1Spectrum){
		int scanNumber = ms1Spectrum.getScanNum();
		int nextScanNumber = nextMs1Spectrum.getScanNum();
		
		ArrayList<TmtPSM> tmtPsms = new ArrayList<TmtPSM>();
		int index = Collections.binarySearch(psms, new TmtPSM(scanNumber));
		index = index < 0 ? -index - 1 : index;
	//	if(index < 0) System.out.println(scanNumber);
		for(int i=index;i<psms.size();i++){
			TmtPSM psm = psms.get(i);
			if(psm.getScanNumber() >= nextScanNumber) break;
			tmtPsms.add(psm);
		}
		return tmtPsms;
	}
	
	public void run(int maxBin){
		Iterator<Spectrum> iterator = new MzXMLSpectraIterator(spectrumFileName, 1, 1);
		
		Spectrum ms1Spectrum = null;
		Spectrum nextMs1Spectrum = null;
		while(iterator.hasNext()){
			ms1Spectrum = nextMs1Spectrum;
			nextMs1Spectrum = iterator.next();
			if(ms1Spectrum == null) continue;
			
			ArrayList<TmtPSM> recruitedPsms = recruitTmtPsms(ms1Spectrum, nextMs1Spectrum);
		//	System.out.println("* " + ms1Spectrum.getScanNum());
			
			for(TmtPSM psm : recruitedPsms){
				//System.out.println(psm.getScanNumber());
				PIPMS1Scorer scorer = new PIPMS1Scorer(ms1Spectrum, psm, tolerance, maxBin);
				psm.setScore(scorer.getScore());
			}
		}
		
		try {
			PrintStream out = new PrintStream(outFileName);
			out.println(TmtPSM.getHeader());
			for(TmtPSM psm : psms){
				out.println(psm);
			}
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
	}
	
	public static void main(String[] args) {
		String specFolder = "/media/kyowon/Data1/TMT/MzXml_phospho/";
		String tsvFolder = "/media/kyowon/Data1/TMT/TSV_phospho/";
		String outFolder = "/media/kyowon/Data1/TMT/TSV_phospho_out/";
		Tolerance tolerance = new Tolerance(10f, true);
		
		if(!new File(outFolder).exists()) new File(outFolder).mkdir();
		
		for(File spec : new File(specFolder).listFiles()){
			String specFilenName = spec.getAbsolutePath();
			
			String tsvFileName = tsvFolder + spec.getName().replace(".mzXML", "Merged.tsv");
			String outFileName = outFolder + spec.getName().replace(".mzXML", "Merged.PIP.tsv");
			FilterWIthPIP test = new FilterWIthPIP(outFileName, tsvFileName, specFilenName, tolerance);
			System.out.println("Processing " + specFilenName);
			test.run(6);
		}
		
		
		

	}
}
