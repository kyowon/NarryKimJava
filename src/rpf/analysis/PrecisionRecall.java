package rpf.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import parser.AnnotationFileParser;
import parser.AnnotationFileParser.AnnotatedGene;
import parser.BedCovFileParser;
import parser.ScoringOutputParser;
import parser.ScoringOutputParser.ScoredPosition;
import rpf.Scorer;

public class PrecisionRecall {
	private static int getTrueSetNumber(AnnotationFileParser annotationParser, BedCovFileParser covPlusParser, BedCovFileParser covMinusParser, int leftWindowSize, int rightWindowSize, int numberOfNonZeroElements){
		HashMap<String, HashMap<Boolean, HashSet<Integer>>> map = new HashMap<String, HashMap<Boolean, HashSet<Integer>>>();
		Iterator<AnnotatedGene> iterator = annotationParser.getAnnotatedGeneIterator();
		while(iterator.hasNext()){
			AnnotatedGene gene = iterator.next();
			int position = gene.isPlusStrand()? gene.getCdsStart() : gene.getCdsEnd();
			BedCovFileParser covParser = gene.isPlusStrand()? covPlusParser : covMinusParser;
			double[] cov = covParser.getCoverages(gene.getContig(), position, leftWindowSize, rightWindowSize, gene.isPlusStrand());
			if(Scorer.numberOfNonZeroElements(cov) < numberOfNonZeroElements) continue;
			
			if(!map.containsKey(gene.getContig())) map.put(gene.getContig(), new HashMap<Boolean, HashSet<Integer>>());
			HashMap<Boolean, HashSet<Integer>> sub = map.get(gene.getContig());
			if(!sub.containsKey(gene.isPlusStrand())) sub.put(gene.isPlusStrand(), new HashSet<Integer>());
			sub.get(gene.isPlusStrand()).add(position);
		}
		
		int sum = 0;
		for(HashMap<Boolean, HashSet<Integer>> sub : map.values()){
			for(HashSet<Integer> ssub : sub.values()){
				sum += ssub.size();
			}
		}		
		return sum;
	}
	
	private static ArrayList<Double> getPrecisionNRecall(ScoringOutputParser parser, int trueSetNumber, double scoreThreshold){
		int totalPositive = 0;
		int truePositive = 0;
		for(ScoredPosition position : parser.getPositions()){
			if(position.getScore() < scoreThreshold) continue;
			totalPositive++;
			if(position.isAnnotated()) truePositive++;
		}
		ArrayList<Double> ret = new ArrayList<Double>();
		ret.add((double)truePositive/totalPositive);
		ret.add((double)truePositive/trueSetNumber);		
		return ret;
	}

	/*
	public static void main(String[] args){
		AnnotationFileParser annotationParser = new AnnotationFileParser("/media/kyowon/Data1/RPF_Project/data/refFlatHuman.txt");
		BedCovFileParser covPlusParser = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10mNew.sorted.plus.cov"); 
		BedCovFileParser covMinusParser = new BedCovFileParser("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10mNew.sorted.minus.cov");
		ScoringOutputParser scoringOutputParser = 
				new ScoringOutputParser("/media/kyowon/Data1/RPF_Project/data/Samfiles/Uncollapsed/Noco_Harr10mNew.sorted.plus.cov.score.tsv");
		int leftWindowSize = 0;
		int rightWindowSize = 50;
		int numberOfNonZeroElements = 5;
		
		int trueSetNumber = PrecisionRecall.getTrueSetNumber(annotationParser, covPlusParser, covMinusParser, leftWindowSize, rightWindowSize, numberOfNonZeroElements);
		System.out.println(trueSetNumber);
		for(double scoreThreshold = 1.3; scoreThreshold <2.5;scoreThreshold+=.05){
			//System.out.println(scoreThreshold + "\t" + PrecisionRecall.getPrecisionNRecall(scoringOutputParser, trueSetNumber, scoreThreshold));
			System.out.println(scoreThreshold + "\t" + PrecisionRecall.getPrecisionNRecall(scoringOutputParser, trueSetNumber, scoreThreshold).get(1));
		}
		
	}*/
	
}
