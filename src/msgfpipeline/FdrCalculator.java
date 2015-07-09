package msgfpipeline;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import fdr.TargetDecoyPSMSet;
import msgfpipeline.parser.MSGFPlusParser;
import msgfpipeline.parser.MSGFPlusParser.MSGFPlusPSM;

public class FdrCalculator {
	private ArrayList<MSGFPlusPSM> allPSMs;
	private String tmpFilename;
	
	public FdrCalculator(ArrayList<MSGFPlusPSM> allPSMs, String tmpFilename){
		this.allPSMs = allPSMs;
		this.tmpFilename = tmpFilename;
	}
	
	public void calculate(){
		PrintStream out;
		try {
			out = new PrintStream(tmpFilename);
			for(MSGFPlusPSM psm : allPSMs){
				out.println(psm);
			}
			out.close();
			TargetDecoyPSMSet psmSet = new TargetDecoyPSMSet(new File(tmpFilename), 
					"\t", false,
					MSGFPlusParser.getSpecEvalueCol(), false, MSGFPlusParser.getScanNumberCol(), 
					MSGFPlusParser.getPeptideCol(), null, MSGFPlusParser.getProteinIDCol(), "XXX");
			for(MSGFPlusPSM psm : allPSMs){
				psm.setPepQvalue(psmSet.getPepFDR(psm.getSpecEvalue()));
				psm.setQvalue(psmSet.getPSMFDR(psm.getSpecEvalue()));
			//	System.out.println(psm.getSpecEvalue());
			}
		//	System.out.println(MSGFPlusParser.getSpecEvalueCol() + " " +  MSGFPlusParser.getScanNumberCol() + " " +MSGFPlusParser.getPeptideCol() +" " + MSGFPlusParser.getProteinIDCol() );
			new File(tmpFilename).delete();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
}
