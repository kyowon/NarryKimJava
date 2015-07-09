package msgfpipeline;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;

public class ProteinListFormatter {
	private static class ProteinSpectraMatch{
		private String proteinID;//
		private String proteinName;//
		private String geneName;//
		private String species;//
		private int length;//
		private int numPSMs;//
		private int numUnsharedPSMs;//
		private double coverageWithSharedPeptides;//
		private double coverageWithoutSharedPeptides;//
		private ArrayList<Integer> rowIndices;//
		
		static public String getHeader(){
			return "ProteinID\tProteinName\tGeneName\tSpecies\t#PSMs\t#Unshared PSMs\tCoverage(w/o shared pep)\tCoverage(w/ shared pep)\tMatchingSpecRowNumbers";
		}
		
		@Override
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(proteinID);sb.append('\t');
			sb.append(proteinName);sb.append('\t');
			sb.append(geneName);sb.append('\t');
			sb.append(species);sb.append('\t');
			sb.append(numPSMs);sb.append('\t');
			sb.append(numUnsharedPSMs);sb.append('\t');
			sb.append(coverageWithoutSharedPeptides);sb.append('\t');
			sb.append(coverageWithSharedPeptides);sb.append('\t');
			for(int i=0; i<rowIndices.size(); i++){
				sb.append(rowIndices.get(i));
				if(i<rowIndices.size()-1) sb.append(';');
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
		
		private ProteinSpectraMatch(String proteinInfo, ArrayList<Integer> rowIndices, ArrayList<int[]> startEnds){
			String[] token = proteinInfo.split("\t");
			this.proteinID = token[0];
			this.proteinName = token[1];
			this.geneName = token[2];
			this.species = token[3];
			this.length = Integer.parseInt(token[4]);
			
			ArrayList<int[]> unsharedStartEnds = getUnsharedPSMStartEnds(startEnds);
			this.coverageWithoutSharedPeptides = getCoverage(length, unsharedStartEnds);
			this.coverageWithSharedPeptides = getCoverage(length, startEnds);
			
			this.numPSMs = startEnds.size();
			this.numUnsharedPSMs = unsharedStartEnds.size();
			this.rowIndices = rowIndices;
		}
		
		
		private ArrayList<int[]> getUnsharedPSMStartEnds(ArrayList<int[]> startEnds){
			return null;
		}
		
		private double getCoverage(int length, ArrayList<int[]> startEnds){
			
			return 0;
		}
		
	}
	
	// protein seqs..
	public ProteinListFormatter(ArrayList<MSGFPlusPSM> psms, String outFileName){
		PrintStream out = null;
		try {
			out = new PrintStream(outFileName);
			out.println(ProteinSpectraMatch.getHeader());
			//HashMap<String, ArrayList<MSGFPlusPSM>> proteinMap = new HashMap<String, ArrayList<MSGFPlusPSM>>();
			HashMap<String, ArrayList<Integer>> rowIndexMap = new HashMap<String, ArrayList<Integer>>();
			HashMap<String, ArrayList<int[]>> startEndMap = new HashMap<String, ArrayList<int[]>>();
			
			for(int i=0;i<psms.size();i++){
				MSGFPlusPSM psm = psms.get(i);
				String[] proteinIDs = psm.getProteinIDs();
				String[] proteinNames = psm.getProteinNames();
				String[] geneNames = psm.getGeneNames();
				String[] species = psm.getSpecies();
				int[] lengths = psm.getProteinLength();
				
				int[] starts = psm.getStartInPro();
				int[] ends = psm.getEndInPro();
				
				for(int j=0;j<proteinIDs.length;j++){
					String proteinInfo = proteinIDs[j] + "\t" + proteinNames[j] + "\t" + geneNames[j] + "\t" + species[j] + "\t" + lengths[j];
					
					if(!rowIndexMap.containsKey(proteinInfo)){
					//	proteinMap.put(proteinInfo, new ArrayList<MSGFPlusPSM>());
						rowIndexMap.put(proteinInfo, new ArrayList<Integer>());
						startEndMap.put(proteinInfo, new ArrayList<int[]>());
					}
				//	proteinMap.get(proteinInfo).add(psm);
					rowIndexMap.get(proteinInfo).add(i+2);
					int[] se = new int[2];
					se[0] = starts[j];
					se[1] = ends[j];
					startEndMap.get(proteinInfo).add(se);
				}
			}
			
			
			for(String proteinInfo : rowIndexMap.keySet()){
				out.println(new ProteinSpectraMatch(proteinInfo, rowIndexMap.get(proteinInfo), startEndMap.get(proteinInfo)));
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
