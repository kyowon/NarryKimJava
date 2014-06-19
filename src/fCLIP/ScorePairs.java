package fCLIP;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;

import parser.MafParser;
import fCLIP.parser.PairScoringOutputParser.ScoredPair;
import fCLIP.parser.ScoringOutputParser;
import fCLIP.parser.ScoringOutputParser.ScoredPosition;

public class ScorePairs {

	static double energyThreshold = -100;
	static int depthThreshold = 20;
	static int hairpinThreshold = 100;
	static double scoreThreshold = .0;
	static int seqLength = 50 + Scorer.flankingNTNumber*2;
	static Classifier classifier;
	
	static String getRandomSequence(int length){
		char[] nu = {'A', 'C', 'G', 'T'};
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<length;i++){
			sb.append(nu[new Random().nextInt(4)]);
		}
		return sb.toString();
	}
	
	static String randomSubstitute(String a, int numSubstitutes){
		String[] nu = {"A", "C", "G", "T"};
		StringBuilder sb = new StringBuilder(a);
		for(int i=0;i<numSubstitutes;i++){
			int index = new Random().nextInt(a.length());
			sb.replace(index, index+1, nu[new Random().nextInt(4)]);
		}
		return sb.toString();
	}
	
	static String permutate(String a){
		ArrayList<Integer> index = new ArrayList<Integer>();
		for(int i=0;i<a.length();i++) index.add(i);
		Collections.shuffle(index);
		StringBuilder sb = new StringBuilder();
		for(int i: index) sb.append(a.charAt(i));
		return sb.toString();
	}
	
	public static void run(String sample, ScoringOutputParser sparser, MafParser mafParser, String arffTrainFileName, boolean randomize) throws IOException {
		//ScoringOutputParser sparser = new ScoringOutputParser( "/media/kyowon/Data1/fCLIP/samples/sample1/bed/" + sample + ".sorted.out.csv");
		PrintStream out = new PrintStream("/media/kyowon/Data1/Dropbox/"+ sample + ".sorted.out.pair"+ (randomize? ".random" : "") + ".csv");
		String arfTestFileName = "/media/kyowon/Data1/Dropbox/"+ sample + ".sorted.out.pair"+ (randomize? ".random" : "") + ".arff";
		PrintStream arffOut = new PrintStream(arfTestFileName);
		
		System.out.println(sample + " pairs");	
		
		classifier = new Classifier(arffTrainFileName);
		ArrayList<ScoredPair> pairs = getScoredPairs(sparser, mafParser, hairpinThreshold, energyThreshold, depthThreshold, randomize);
		arffOut.println(ScoredPosition.getArffHeader());
		out.println(ScoredPair.getHeader());
		
		for(ScoredPair pair : pairs){
			classifier.setClassification(pair);
			out.println(pair);
			arffOut.println(pair.toArffString());
		}
				
		arffOut.close();
		out.close();
	}
	
	private static ArrayList<ScoredPair> getScoredPairs(ScoringOutputParser sparser, MafParser mafParser, int hairpinThreshold, double energyThreshold, int depthThreshold, boolean randomize){
		ArrayList<ScoredPair> allPairs = new ArrayList<ScoredPair>();
		ArrayList<ScoredPosition> allPositions = new ArrayList<ScoredPosition>();
		HashSet<String> seqs = new HashSet<String>();
		
		int three = 0, five =0;
		for(ScoredPosition sp : sparser.getPositions()){
			if(sp.isPaired()) continue;
			if(sp.getThreePScore() < scoreThreshold && sp.getFivePScore() < scoreThreshold) continue;
			if(sp.getSeq() == null || sp.getSeq().length() < Scorer.flankingNTNumber) continue;
			if(seqs.contains(sp.getSeq())) continue;
			if(sp.getHairpinNumber() > hairpinThreshold) continue;
			if(sp.getEnergy() < energyThreshold) continue;
			if(sp.getDepth() > depthThreshold) continue;
			if(sp.getMiRNAs() != null && !sp.getMiRNAs().isEmpty()) continue;
			if(sp.getClassification().equals("M")) continue;
			
			allPositions.add(sp);
			seqs.add(sp.getSeq());
			
			if(sp.is3pScored()) three++;
			else five++;
		}		
		
		if(randomize) RandomizeSequences(allPositions, energyThreshold, depthThreshold, hairpinThreshold);
		System.out.println(sparser.getPositions().size() + " " + allPositions.size() + " " + five + " " + three);
		
		HashSet<String> mergedSeqs = new HashSet<String>();
		for(int i=0;i<allPositions.size();i++){
			System.out.println(i + " " + allPositions.size());
			ScoredPosition p1 = allPositions.get(i);
			boolean is3pScored = p1.is3pScored();
			double minEnergy = 10000;
			ScoredPair maxScoredPair = null;
			
			for(int j=0;j<allPositions.size();j++){
				if(i == j) continue;
				ScoredPosition p2 = allPositions.get(j);
				if(is3pScored == p2.is3pScored()) continue;
				
				if(is3pScored){
					if(p1.getContig().equals(p2.getContig()) && Math.abs(p1.getThreePposition() - p2.getFivePposition())<1000) continue;
				}else{
					if(p1.getContig().equals(p2.getContig()) && Math.abs(p2.getThreePposition() - p1.getFivePposition())<1000) continue;
				}
				
				ScoredPair pair = new ScoredPair(p1, p2, seqLength);
				
				if(pair.getEnergy() < minEnergy){ //TDOD
					minEnergy = pair.getEnergy();
					maxScoredPair = pair;
				}
			}			
			if(maxScoredPair == null) continue;
			if(mergedSeqs.contains(maxScoredPair.getSeq())) continue;		
			
			if(!randomize){
				maxScoredPair.setDnDs(mafParser);
				maxScoredPair.setRNAzScores(mafParser);
			}
			allPairs.add(maxScoredPair);		
			mergedSeqs.add(maxScoredPair.getSeq());
			
		}		
		return allPairs;
		
	}
	
	static private void RandomizeSequences(ArrayList<ScoredPosition> allPositions, double energyThreshold, int depthThreshold, int hairpinThreshold){
		int m = 0;
		for(ScoredPosition sp : allPositions){
			boolean isSet = false;
			int cntr = 0;
			while(!isSet){
				if(cntr++ > 5000) break;
				String rs = permutate(sp.getSeq());
				rs = randomSubstitute(rs, 10);
				RNAfoldLauncher rnafold = new RNAfoldLauncher(rs, Scorer.flankingNTNumber);
				if(rnafold.getNumberOfHairPins() > hairpinThreshold) continue;
				if(rnafold.getEnergy() < energyThreshold) continue;
				if(rnafold.getDepth() > depthThreshold) continue;
				
			//	double prediction = classifier.classify(rnafold.toInstance(classifier.getDataset()));
			//	if(prediction > -0.99) continue;
				
				sp.setSeq(rs);
				isSet = true;
			}
			System.out.println("Randomizing " + m++ + " out of " + allPositions.size());
		}
	}
	
	

}
