package msgfpipeline;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;
import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;

public class ProteinListFormatter {
	private static class ProteinSpectraMatch{
		private String proteinID;//
		private String proteinName;//
		private String geneName;//
		private String species;//
		private ArrayList<String> specIDsWithSharedPeptides;
		private ArrayList<String> specIDsWithoutSharedPeptides;
		private ArrayList<int[]> startEndsWithSharedPeptides;
		private ArrayList<int[]> startEndsWithoutSharedPeptides;
		private double coverageWithSharedPeptides;//
		private double coverageWithoutSharedPeptides;//
		private float proteinMass;
		private int proteinLength;
		private int startPositioninSuffixArray;
		
		static public String getHeader(){
			return "ProteinID\tProteinName\tProteinLength\tProteinMass\tGeneName\tSpecies\t#PSMs\t#Unshared PSMs\tCoveredPositions(w/ shared pep)\tCoveredPositions(w/o shared pep)\tCoverage(w/ shared pep)\tCoverage(w/o shared pep)\tMatchingSpecFile&ScanNumber(w/ shared pep)\tMatchingSpecFile&ScanNumber(w/o shared pep)";
		}
		
		public String toFastaString(SuffixArray sa, boolean dispShared){
			StringBuilder sb = new StringBuilder();
			sb.append(">");
			sb.append(sa.getAnnotation(startPositioninSuffixArray));
			sb.append('\n');
			
			SuffixArraySequence sequence = sa.getSequence();
			String fullSeq = sequence.getMatchingEntry(startPositioninSuffixArray);
			char[] fs = fullSeq.toLowerCase().toCharArray();
			ArrayList<int[]> ses = dispShared ? startEndsWithSharedPeptides : startEndsWithoutSharedPeptides;
			if(ses != null){
				for(int[] se : ses){
					for(int l = se[0]; l<se[1];l++){
						fs[l-1] = Character.toUpperCase(fs[l-1]);
					}
				}
			}
			for(char c : fs){
				sb.append(c);
			}
			return sb.toString();
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(proteinID);sb.append('\t');
			sb.append(proteinName);sb.append('\t');
			sb.append(proteinLength);sb.append('\t');
			sb.append(proteinMass);sb.append('\t');
			sb.append(geneName);sb.append('\t');
			sb.append(species);sb.append('\t');
			sb.append(specIDsWithSharedPeptides.size());sb.append('\t');
			sb.append(specIDsWithoutSharedPeptides == null ? 0 : specIDsWithoutSharedPeptides.size());sb.append('\t');
			for(int[] interval : startEndsWithSharedPeptides){
				sb.append(interval[0]);sb.append('-');
				sb.append(interval[1]-1);sb.append(';');
			}
			sb.append('\t');
			if(startEndsWithoutSharedPeptides != null){
				for(int[] interval : startEndsWithoutSharedPeptides){
					sb.append(interval[0]);sb.append('-');
					sb.append(interval[1]-1);sb.append(';');
				}
			}
			sb.append('\t');
			sb.append(coverageWithSharedPeptides);sb.append('\t');
			sb.append(coverageWithoutSharedPeptides);sb.append('\t');
			for(String specID : specIDsWithSharedPeptides){
				sb.append(specID);sb.append(';');
			}
			sb.append('\t');
			if(specIDsWithoutSharedPeptides != null){
				for(String specID : specIDsWithoutSharedPeptides){
					sb.append(specID);sb.append(';');
				}
			}
			return sb.toString();
		}
		
		@Override
		public int hashCode(){
			return proteinID.hashCode();	
		}
		
		@Override
		public boolean equals(Object o){
			return proteinID.equals(o);
		}
		
		private ProteinSpectraMatch(String proteinID, ArrayList<MSGFPlusPSM> psms, ArrayList<MSGFPlusPSM> unsharedPsms){ // unsharedPsms may be null
			this.proteinID = proteinID;
			MSGFPlusPSM psm = psms.get(0);
			int i = psm.getProteinIndex(proteinID);
			this.startPositioninSuffixArray = psm.getStartInSuffixArray()[i];
			this.proteinName = psm.getProteinNames()[i];
			this.proteinLength = psm.getProteinLength()[i];
			this.proteinMass = psm.getProteinMass()[i];
			this.geneName = psm.getGeneNames()[i];
			this.species = psm.getSpecies()[i];
			this.specIDsWithSharedPeptides = new ArrayList<String>();
			for(MSGFPlusPSM p : psms){
				this.specIDsWithSharedPeptides.add(p.getSpecFile()+"-"+p.getScanNumber());
			}
			if(unsharedPsms != null){
				this.specIDsWithoutSharedPeptides = new ArrayList<String>();
				for(MSGFPlusPSM p : unsharedPsms){
					this.specIDsWithoutSharedPeptides.add(p.getSpecFile()+"-"+p.getScanNumber());
				}
			}
			this.startEndsWithSharedPeptides = getPSMStartEnds(proteinID, proteinLength, psms);
			this.startEndsWithoutSharedPeptides = getPSMStartEnds(proteinID, proteinLength, unsharedPsms);
			this.coverageWithSharedPeptides = getCoverage(proteinLength, this.startEndsWithSharedPeptides);
			this.coverageWithoutSharedPeptides = getCoverage(proteinLength, this.startEndsWithoutSharedPeptides);
		}
		
		
		private ArrayList<int[]> getPSMStartEnds(String proteinID, int length, ArrayList<MSGFPlusPSM> psms){
			if(psms == null) return null;
			ArrayList<int[]> ses = new ArrayList<int[]>();
			BitSet bs = new BitSet(length);
			for(MSGFPlusPSM psm : psms){
				int i = psm.getProteinIndex(proteinID);
				int s = psm.getStartInPro()[i];
				int e = psm.getEndInPro()[i];
				for(int j=s;j<e;j++) bs.set(j);
			}
			int j=0;
			for (int i = bs.nextSetBit(0); i >= 0 && j >=0; i = bs.nextSetBit(j)) {
				j = bs.nextClearBit(i);
				int[] se = new int[2];
				se[0] = i; se[1] = j;
				ses.add(se);
			 }
			return ses;
			
		}
		
		private double getCoverage(int length, ArrayList<int[]> startEnds){//startEnds should be exclusive each other TODO check if start ends are inclusive or not
			if(startEnds == null) return 0;
			int cov = 0;
			for(int[] interval : startEnds){
				cov += interval[1] - interval[0];
			}
			return (100.0*cov)/length;
		}
		
	}
	
	// protein seqs..
	public ProteinListFormatter(ArrayList<MSGFPlusPSM> psms, String originalFasta, String outFileName, String outFastaNameShared, String outFastaNameUnshared){
		PrintStream out = null;
		try {
			out = new PrintStream(outFileName);
			out.println(ProteinSpectraMatch.getHeader());
			HashMap<String, ArrayList<MSGFPlusPSM>> proteinMap = new HashMap<String, ArrayList<MSGFPlusPSM>>();
			HashMap<String, ArrayList<MSGFPlusPSM>> proteinUnsharedMap = new HashMap<String, ArrayList<MSGFPlusPSM>>();
			
			for(int i=0;i<psms.size();i++){
				MSGFPlusPSM psm = psms.get(i);
				for(String proteinID : psm.getProteinIDs()){
					if(!proteinMap.containsKey(proteinID)) proteinMap.put(proteinID, new ArrayList<MSGFPlusPSM>());
					proteinMap.get(proteinID).add(psm);
				}
				if(psm.getProteinIDs().length == 1){ // unshared peptide
					for(String proteinID : psm.getProteinIDs()){
						if(!proteinUnsharedMap.containsKey(proteinID)) proteinUnsharedMap.put(proteinID, new ArrayList<MSGFPlusPSM>());
						proteinUnsharedMap.get(proteinID).add(psm);
					}
				}
			}
			
			SuffixArraySequence sequence = new SuffixArraySequence(originalFasta);
			SuffixArray sa = new SuffixArray(sequence);
			
			PrintStream outfastashared = new PrintStream(outFastaNameShared);
			PrintStream outfastaunshared  = new PrintStream(outFastaNameUnshared);
			
			for(String proteinID : proteinMap.keySet()){
				ProteinSpectraMatch prosm = new ProteinSpectraMatch(proteinID, proteinMap.get(proteinID), proteinUnsharedMap.get(proteinID));
				out.println(prosm);
				outfastashared.println(prosm.toFastaString(sa, true));
				outfastaunshared.println(prosm.toFastaString(sa, false));
			}
			outfastashared.close();
			outfastaunshared.close();
			out.close();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
}
