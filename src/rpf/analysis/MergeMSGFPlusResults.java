package rpf.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import parser.MSGFPlusParser;
import parser.MSGFPlusParser.MSGFPlusPSM;

public class MergeMSGFPlusResults {
	public static void main(String[] args){
		String key = "uORF";
		File dir = new File("/media/kyowon/Data1/RPF_Project/data/MSMS/HeLa_RP");
		ArrayList<MSGFPlusPSM> allPsms = new ArrayList<MSGFPlusPSM>();
		for(File file : dir.listFiles()){
			if(!file.getAbsolutePath().endsWith(".tsv")) continue;
			MSGFPlusParser parser = new MSGFPlusParser(file.getAbsolutePath());
			ArrayList<MSGFPlusPSM> psms = parser.getPsms();
			for(MSGFPlusPSM psm : psms){
				if(psm.getPepQvalue() > 0.01 || psm.getSpecEvalue() > 1e-11) continue;
				if(psm.getProtein().contains(key))
					allPsms.add(psm);// || NPLPSKETIEQEK
			}
		}
		
		HashSet<String> peptides = new HashSet<String>();
		for(MSGFPlusPSM psm : allPsms){
			System.out.println(psm);
			peptides.add(psm.getPeptide());
		}
		
		int i=0;
		for(String peptide : peptides){
			System.out.println(">" + i++ + "_" + key + "\n" + peptide.replace("C+57.021", "C"));
		}	
	}
}
