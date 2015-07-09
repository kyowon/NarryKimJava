package fCLIP.analysis;

import java.util.ArrayList;
import java.util.HashSet;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class CheckList {

	public static void main(String[] args) {
		String cis = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/cis/merged.beforeBlatCsv";
		String trans = "/media/kyowon/Data1/fCLIP/samples/sample4/results/lists/trans/merged.trans.complete.M.csv";
		int numMiRNAs = 0; //
		int novelCis = 0; //
		int novelTrans = 0; //
		int cisOverhang2 = 0; //
		int numAlu = 0;
		double transThreshold = 0.0;
		double predictionThreshold = .7;
		
		int transOverhang2 = 0; //
		ArrayList<ScoredPosition> snoRNAs = new ArrayList<ScoredPosition>();//
		ArrayList<ScoredPosition> rRNAs = new ArrayList<ScoredPosition>();//
		String[] cisKnownGenes = {"AURKB", "STAT6", "MIR17HG", "SH2B1", "SPHK", "ANKS6", "DROSHA", "CBX6"};//
		String[] transKnownGenes = {"CERS5", "SCAF4", "CBX6", "ACSL3", "PCBP2", "SPECC1L", "MAPK8IP3", "KNTC1", "TRAPPC3", "IRF3", "DROSHA"};//		
		
		HashSet<String> cisExpGenes = new HashSet<String>();
		HashSet<String> transExpGenes = new HashSet<String>();
		
		
		ScoredPositionOutputParser cisParser = new ScoredPositionOutputParser(cis);
		for(ScoredPosition sp : cisParser.getPositions()){
			if(!sp.getClassification().equals("M")) continue;
			if(sp.getPredictionScore() < predictionThreshold) continue;
			if(sp.hasMatchingMiRNA()){
				numMiRNAs++;
				//continue;
			}else{
				novelCis++;
				if(sp.getOverHang() == 2) cisOverhang2 ++;
				else continue;
				
			}
			ArrayList<String> genes = sp.getContainingGeneNames();
			if(genes == null || genes.isEmpty()) continue;
			boolean isSnoRNA = false;
			boolean isRRNA = false;
			boolean known = false;
			for(String gene : genes){
				if(gene.toLowerCase().startsWith("sno")) isSnoRNA = true;
				if(gene.toLowerCase().startsWith("rrna")) isRRNA = true;
				for(String kgene : cisKnownGenes){
					if(gene.toUpperCase().startsWith(kgene)){
						known = true;
						cisExpGenes.add(gene);
						//break;
					}
				}
			}
			
			if(isSnoRNA) snoRNAs.add(sp);
			if(isRRNA) rRNAs.add(sp);
			//if(known) System.out.println(sp);
			
		}
		System.out.println("# miRNAs: " + numMiRNAs);
		System.out.println("# novel cis duplexes: " + novelCis);
		System.out.println("# overhang == 2 cis duplexes: " + (cisOverhang2));
		System.out.println("After overhang == 2 filtration...");
		//System.out.println("# novel Cis: " + novelCis);
		System.out.println("# snoRNAs: " + snoRNAs.size());
		System.out.println("# rRNAs: " + rRNAs.size());
		System.out.println("Detected Genes: " + cisExpGenes);
	//	int top = 0;
	//	int alu = 0;
	//	int topalu = 0;
	//	int topoverhang2 = 0;
		//System.exit(1);
		ScoredPairOutputParser transParser = new ScoredPairOutputParser(trans);
		for(ScoredPair pair : transParser.getPairs()){
			if(pair.getFivePScore() < transThreshold || pair.getThreePScore() < transThreshold) continue;
			if(pair.getPredictionScore() < predictionThreshold) continue;
			novelTrans++;
			if(pair.getOverHang() == 2) transOverhang2++;
			else continue;
			
			if(pair.getMisc().contains("Alu")) numAlu++;
			ScoredPosition[] pairedPositions = pair.getPairedScoredPositions();
			boolean known = false;
			for(ScoredPosition pp : pairedPositions){
				if(known) break;
				ArrayList<String> genes = pp.getContainingGeneNames();
				if(genes == null || genes.isEmpty()) continue;
				for(String gene : genes){
					if(known) break;
					for(String kgene : transKnownGenes){
						if(gene.toUpperCase().startsWith(kgene)){
							known = true;
							//break;
							transExpGenes.add(gene);
						}
					}
				}
			}
			//if(known) System.out.println(pair);			
		}
		
		System.out.println();
		System.out.println("Trans duplex score threshold: " + transThreshold);
		System.out.println("# novel trans duplexes: " + novelTrans);		
		System.out.println("# overhang == 2 trans duplexes: " + (transOverhang2));
		System.out.println("After overhang == 2 filtration...");
		
		System.out.println("Detected Genes: " + transExpGenes);
		System.out.println("# Alu: " + numAlu);
	//	System.out.println(top + " " + alu + " " + topalu + " " + topoverhang2);
		
	}

}
