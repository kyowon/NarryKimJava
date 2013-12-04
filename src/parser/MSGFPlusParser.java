package parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import net.sf.samtools.util.BufferedLineReader;

public class MSGFPlusParser {
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
		private String protein;
		private int deNovoScore;
		private int msgfScore;
		private float specEvalue;
		private float evalue;
		private float qvalue = -1; // optional
		private float pepQvalue = -1; // optional
		
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

		public String getProtein() {
			return protein;
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
			protein = token[9];
			deNovoScore = Integer.parseInt(token[10]);
			msgfScore = Integer.parseInt(token[11]);
			specEvalue = Float.parseFloat(token[12]);
			evalue = Float.parseFloat(token[13]);
			if(token.length > 14){
				qvalue = Float.parseFloat(token[14]);
				pepQvalue = Float.parseFloat(token[15]);
			}
		}
		
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
			ret.append(protein);ret.append('\t');
			ret.append(deNovoScore);ret.append('\t');
			ret.append(msgfScore);ret.append('\t');
			ret.append(specEvalue);ret.append('\t');
			ret.append(evalue);
			if(qvalue > -1){
				ret.append('\t'); ret.append(qvalue);
				ret.append('\t'); ret.append(pepQvalue);
			}
			return ret.toString();
		}	
	}
	
	private ArrayList<MSGFPlusPSM> psms;
	public MSGFPlusParser(String fileName){
		psms = new ArrayList<MSGFPlusPSM>();
		try {
			BufferedLineReader in = new BufferedLineReader(new FileInputStream(fileName));
			String s;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				psms.add(new MSGFPlusPSM(s));
			}in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}	
	
	public ArrayList<MSGFPlusPSM> getPsms() {
		return psms;
	}
}
