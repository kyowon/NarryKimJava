package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;

public class TransDuplexFinder {
	
	public static void main(String[] args) { // give options for num 5p 3p paired... then analyze ... 
		if(args.length == 9){
			String outFileNameM = args[0];
			String outFileNameU = args[1];
			String cisOutFileName = args[2];
			String arffTrainFileName = args[3];
			int blatHitThreshold = Integer.parseInt(args[4]);
			int seqLength = Integer.parseInt(args[5]);
			boolean sameDirection = Boolean.parseBoolean(args[6]);
			int numPositions = Integer.parseInt(args[7]);
			int numThreads = Integer.parseInt(args[8]);
			
			ScoredPositionOutputParser sparser = new ScoredPositionOutputParser(cisOutFileName);
			System.out.println("Scoring " + outFileNameM + " pairs");	
			Classifier classifier = new Classifier(arffTrainFileName);
			String tmpM = outFileNameM + ".tmp.csv";
			String tmpU = outFileNameU + ".tmp.csv";
			//if(!new File(tmpM).exists() || !new File(tmpU).exists() )
			FCLIP_ScorePairs.getScoredPairs(sparser, classifier, blatHitThreshold, seqLength, sameDirection, tmpM, tmpU, numPositions, numThreads);
			FCLIP_ScorePairs.setMatchedNumPositions(tmpM, outFileNameM);
			FCLIP_ScorePairs.setMatchedNumPositions(tmpU, outFileNameU);
			
			new File(tmpM).delete();
			new File(tmpU).delete();
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
