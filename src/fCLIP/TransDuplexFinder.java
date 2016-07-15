package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import launcher.RNAcofoldLauncher;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;

public class TransDuplexFinder {
	
	public static void main(String[] args) { // give options for num 5p 3p paired... then analyze ... 
		if(args.length == 10){
			String outFileNameM = args[0];
			//String outFileNameU = args[1];
			String cisOutFileName = args[1];
			String arffTrainFileName = args[2];
			int blatHitThreshold = Integer.parseInt(args[3]);
			
			int flankingLength = Integer.parseInt(args[4]);
			int minpre = Integer.parseInt(args[5]);
			int maxpre = Integer.parseInt(args[6]);
			
			FCLIP_Scorer.setFlankingNTNumber(flankingLength);
			FCLIP_Scorer.setMaxReadDiff(maxpre);
			FCLIP_Scorer.setMinReadDiff(minpre);
			RNAcofoldLauncher.setSeqLength(maxpre);			
			
			int seqLength = (80 + flankingLength * 2); 
			
			boolean sameDirection = Boolean.parseBoolean(args[7]);
			int numPositions = Integer.parseInt(args[8]);
			int numThreads = Integer.parseInt(args[9]);
			
			ScoredPositionOutputParser sparser = new ScoredPositionOutputParser(cisOutFileName);
			System.out.println("Scoring " + outFileNameM + " pairs");	
			Classifier classifier = new Classifier(arffTrainFileName);
			String tmpM = outFileNameM + ".tmp.csv";
			//String tmpU = outFileNameU + ".tmp.csv";
			if(!new File(tmpM).exists())
				FCLIP_ScorePairs.getScoredPairs(sparser, classifier, blatHitThreshold, seqLength, sameDirection, tmpM, numPositions, numThreads);
			FCLIP_ScorePairs.setMatchedNumPositions(tmpM, outFileNameM);
			//FCLIP_ScorePairs.setMatchedNumPositions(tmpU, outFileNameU);
			
			new File(tmpM).delete();
			//new File(tmpU).delete();
		}else if(args.length == 4){
			String outCsv = args[0];
			String pairCsv = args[1];
			int num3pPaired = Integer.parseInt(args[2]);
			int num5pPaired = Integer.parseInt(args[3]);
			
			ScoredPairOutputParser parser = new ScoredPairOutputParser(pairCsv);
			PrintStream out;
			int sum = 0;
			int o2 = 0;
			try {
				out = new PrintStream(outCsv);
				out.println(ScoredPair.getHeader());
				for(ScoredPair p : parser.getPairs()){
					if(p.getMatched3pNum() > num3pPaired && p.getMatched5pNum() > num5pPaired){
						out.println(p);
						sum ++;
						if(p.getOverHang() == 2) o2 ++;
					}
				}
				System.out.println("Total pairs : " + sum + " overhang portion: " + (double)o2/sum * 100);
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}	
		}
	}
}
