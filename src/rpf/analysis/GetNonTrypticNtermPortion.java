package rpf.analysis;

import java.io.File;
import java.util.HashMap;

import parser.MSGFPlusParser;
import parser.MSGFPlusParser.MSGFPlusPSM;

public class GetNonTrypticNtermPortion {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		HashMap<String, Character> peptideNPreAAMap = new HashMap<String, Character>();
		HashMap<Character, int[]> numMap = new HashMap<Character, int[]>();
		File dir = new File(args[0]);
		
		for(File file : dir.listFiles()){
			if(!file.getName().endsWith(".tsv")) continue;
			MSGFPlusParser parser = new MSGFPlusParser(file.getAbsolutePath());
			for(MSGFPlusPSM psm : parser.getPsms()){
				if(psm.getPepQvalue() > 0.01) continue;
				if(psm.getQvalue() > 0.05) continue;
				peptideNPreAAMap.put(psm.getPeptide(), psm.getPreAA());
				//System.out.println(psm.getPreAA());
			}				
		}
		
		for(String peptide : peptideNPreAAMap.keySet()){
			char nAA = peptide.charAt(0);
			char preAA = peptideNPreAAMap.get(peptide);
			if(!numMap.containsKey(nAA)) numMap.put(nAA, new int[2]);
			numMap.get(nAA)[0]++;
			if(preAA != 'K' && preAA != 'R' && preAA != '-') numMap.get(nAA)[1]++;
		}
		
		for(char aa : numMap.keySet()){
			int[] nn = numMap.get(aa);
			System.out.println(aa + " " + nn[0] + " " + nn[1] + " " + ((double)nn[1]/nn[0]*100));
		}
		
	}

}
