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
	private String[] excludeList;
			
	public FdrCalculator(ArrayList<MSGFPlusPSM> allPSMs, String tmpFilename){
		this.allPSMs = allPSMs;
		this.tmpFilename = tmpFilename;
	}
	
	public FdrCalculator(ArrayList<MSGFPlusPSM> allPSMs, String tmpFilename, String[] excludeList){
		this.allPSMs = allPSMs;
		this.tmpFilename = tmpFilename;
		this.excludeList = excludeList;		
	}
	
	
	public void calculate(){
		PrintStream out;
		try {
			out = new PrintStream(tmpFilename);
			ArrayList<MSGFPlusPSM> tmpPSMs = new ArrayList<MSGFPlusPSM>(allPSMs);
			for(MSGFPlusPSM psm : tmpPSMs){
				if(excludeList != null){
					ArrayList<String> newProteinIDs = new ArrayList<String>();
					for(String pid : psm.getProteinIDs()){
						boolean exclude = false;
						for(String ep : excludeList){
							if(pid.startsWith(ep) || pid.startsWith("XXX_"+ep)){
								//System.out.println(pid);
								exclude = true;
								break;
							}
						}
						if(!exclude) newProteinIDs.add(pid);
					}
					if(newProteinIDs.isEmpty()) continue;
					else{
						String[] newPIDs = new String[newProteinIDs.size()];
						for(int i=0;i<newPIDs.length;i++) newPIDs[i] = newProteinIDs.get(i);
						psm.setProteinIDs(newPIDs);
					}
				}
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
			}
				new File(tmpFilename).deleteOnExit();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
}
