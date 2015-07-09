package fCLIP.analysis;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class FilterResults {

	public static void main(String[] args) {
		String cis = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/cis/merged.complete.csv";
		String trans = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/trans/merged.trans.complete.M.csv";
		
		String cisOutFile = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/cis/merged.complete.filtered.csv";
		String transOutFile = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/trans/merged.trans.complete.M.filtered.csv";
		
		boolean overhangFilter = true;
		boolean classificationFilter = true;
		boolean filterIntergenic = true;
		double transScoreThreshold = 0;
		
		ScoredPositionOutputParser cisParser = new ScoredPositionOutputParser(cis);
		try {
			PrintStream cisOut = new PrintStream(cisOutFile);
			cisOut.println(cisParser.getHeader());
			for(ScoredPosition sp : cisParser.getPositions()){
				if(overhangFilter && sp.getOverHang() !=2) continue;
				if(classificationFilter && !sp.getClassification().equals("M")) continue;
				if(filterIntergenic && (sp.getContainingGeneAccessions() == null || sp.getContainingGeneAccessions().isEmpty())) continue;
				cisOut.println(sp);
			}			
			cisOut.close();
			
			ScoredPairOutputParser transParser = new ScoredPairOutputParser(trans);
			PrintStream transOut = new PrintStream(transOutFile);
			transOut.println(transParser.getHeader());
			for(ScoredPair pair : transParser.getPairs()){
				if(overhangFilter && pair.getOverHang() !=2) continue;
				if(filterIntergenic){
					if(pair.getPairedScoredPositions()[0].getContainingGeneAccessions() == null || 
							pair.getPairedScoredPositions()[0].getContainingGeneAccessions().isEmpty()) continue;
					if(pair.getPairedScoredPositions()[1].getContainingGeneAccessions() == null || 
							pair.getPairedScoredPositions()[1].getContainingGeneAccessions().isEmpty()) continue;
				}
				if(pair.getFivePScore() < transScoreThreshold || pair.getThreePScore() < transScoreThreshold) continue;
				transOut.println(pair);
			}
			
			transOut.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 

	}

}
