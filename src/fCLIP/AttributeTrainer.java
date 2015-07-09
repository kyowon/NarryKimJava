package fCLIP;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import launcher.RNAcofoldLauncher;
import parser.AnnotationFileParser;
import parser.Bed12Parser;
import parser.MirGff3FileParser;
import parser.ZeroBasedFastaParser;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class AttributeTrainer {
	static public void main(String[] args){
		String outFileName = args[0];
		String arffTrainOutFileName = args[1];
		String bedFileName = args[2];
		ZeroBasedFastaParser fastaParser = new ZeroBasedFastaParser(args[3]);
		AnnotationFileParser annotationParser = new AnnotationFileParser(args[4]);
		MirGff3FileParser mirParser = new MirGff3FileParser(args[5]);
		String parameterFileName = args[6];
		double unpairedScoreThreshold = Double.parseDouble(args[7]);
		double pairedScoreThreshold = Double.parseDouble(args[8]);
		FCLIP_Scorer.setFlankingNTNumber(Integer.parseInt(args[9]));
		int minpre =  Integer.parseInt(args[10]);
		int maxpre =  Integer.parseInt(args[11]);
		RNAcofoldLauncher.setSeqLength(minpre);
		
		FCLIP_Scorer.setMaxReadDiff(maxpre);
		FCLIP_Scorer.setMinReadDiff(minpre);
		
		System.out.println("Training Classifier for " + bedFileName);
		HashSet<ScoredPosition> trainingPositions = new HashSet<ScoredPosition>();
		int annotated = 0;
		for(String contig : fastaParser.getContigs()){
			if(contig.length() > 5 || contig.equals("chrM")) continue;
			Bed12Parser bedParser = new Bed12Parser(bedFileName, annotationParser, contig, true);
			FCLIP_Scorer scorer = new FCLIP_Scorer(bedParser, fastaParser, mirParser, parameterFileName);
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, pairedScoreThreshold, true)){
				trainingPositions.add(sp);
			}
			for(ScoredPosition sp : scorer.getScoredPositions(annotationParser, null, unpairedScoreThreshold, pairedScoreThreshold, false)){
				trainingPositions.add(sp);
			}		
		}	
		
		try {
			PrintStream outTrainArff = new PrintStream(arffTrainOutFileName);
			PrintStream outTrain = new PrintStream(outFileName);
			outTrainArff.println(Classifier.getArffHeader());
			outTrain.println(ScoredPositionOutputParser.ScoredPosition.getHeader());
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() != null) 	annotated++;
			}
			
			System.out.println("Training size : " + trainingPositions.size() + " Annotated : " + annotated);
			int numUnannotated = 300;
			for(ScoredPosition sp : trainingPositions){
				if(sp.getMiRNAs() == null) numUnannotated--;
				if(numUnannotated < 0 && sp.getMiRNAs() == null) continue;
			
				outTrainArff.println(Classifier.toTrainingArffString(sp));
				outTrain.println(sp);
			}
			outTrainArff.close(); // training done		
			outTrain.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
}