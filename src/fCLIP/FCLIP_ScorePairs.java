package fCLIP;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.HashSet;
import java.util.Random;
import java.util.Vector;

import launcher.RNAfoldLauncher;
import fCLIP.parser.ScoredPairOutputParser;
import fCLIP.parser.ScoredPairOutputParser.ScoredPair;
import fCLIP.parser.ScoredPositionOutputParser;
import fCLIP.parser.ScoredPositionOutputParser.ScoredPosition;

public class FCLIP_ScorePairs {
		
	public static class ScorePairRunner implements Runnable{
		private Vector<ScoredPosition> allPositions;
		private Vector<ScoredPosition> compairngPositions;
		private boolean sameDirection;
		private int seqLength;
		private String mf;
		private Classifier classifier;
		
		ScorePairRunner(String mf, ArrayList<ScoredPosition> compairngPositions, ArrayList<ScoredPosition> allPositions, 
				int seqLength, Classifier classifier, boolean sameDirection){
	//		this.uf = uf;
			this.mf = mf;
			this.compairngPositions = new Vector<ScoredPosition>(compairngPositions);
			this.allPositions = new Vector<ScoredPosition>(allPositions);
			this.seqLength = seqLength;
			this.classifier = classifier;
			this.sameDirection = sameDirection;
		}
		
		public void run() {
			try {
				//PrintStream outU = new PrintStream(uf);
				PrintStream outM = new PrintStream(mf);
				//HashSet<String> mergedSeqs = new HashSet<String>();
				int n = 0;
						
				for(int i=0;i<compairngPositions.size();i++){
					ScoredPosition p1 = compairngPositions.get(i);
					boolean is3pScored = p1.is3pScored();
					if(!is3pScored) continue;  // p1 is 3pscored
					//System.out.println(n++ + " " + compairngPositions.size() + " 3p");
					
					for(int j=0;j<allPositions.size();j++){
						if(i == j) continue;
						ScoredPosition p2 = allPositions.get(j);
						if(!sameDirection){
							if(!p2.is5pScored()) continue;
						}else{
							if(!p2.is3pScored()) continue;
						}
						
						if(!sameDirection){
							if(p1.getContig().equals(p2.getContig()) && Math.abs(p1.getThreePPosition() - p2.getFivePPosition())<=150) continue;
						}else{
							if(p1.getContig().equals(p2.getContig()) && Math.abs(p1.getThreePPosition() - p2.getThreePPosition())<=150) continue;
						}
						
						ScoredPair pair = new ScoredPair(p1, sameDirection? ScoredPosition.get3p5pSwappedPosition(p2, FCLIP_Scorer.getFlankingNTNumber()) : p2, seqLength);
						//if(!pair.isFeasibleFold()) continue;
						//if(pair.getOverHang() < 0 || pair.getOverHang() > 5) continue;
						//if(mergedSeqs.contains(pair.getSeq())) continue;	
						
						classifier.setClassification(pair);
						if(pair.getClassification().equals("M")) outM.println(pair);
						//else outU.println(pair);
						//mergedSeqs.add(pair.getSeq());
					}	
				}				
			    n = 0;
				for(int i=0;i<compairngPositions.size();i++){
					ScoredPosition p1 = compairngPositions.get(i);
					boolean is5pScored = p1.is5pScored();
					if(!is5pScored) continue;  // p1 is 5pscored
					//System.out.println(n++ + " " + compairngPositions.size() + " 5p");
					
					for(int j=0;j<allPositions.size();j++){
						if(i == j) continue;
						ScoredPosition p2 = allPositions.get(j);
						if(!sameDirection){
							if(!p2.is3pScored()) continue;
						}else{
							if(!p2.is5pScored()) continue;
						}
						
						if(!sameDirection){
							if(p1.getContig().equals(p2.getContig()) && Math.abs(p2.getThreePPosition() - p1.getFivePPosition())<1000) continue;
						}else{
							if(p1.getContig().equals(p2.getContig()) && Math.abs(p1.getFivePPosition() - p2.getFivePPosition())<1000) continue;
						}
						
						ScoredPair pair = new ScoredPair(p2, sameDirection? ScoredPosition.get3p5pSwappedPosition(p1, FCLIP_Scorer.getFlankingNTNumber()) : p1, seqLength);
						if(pair.getOverHang() < RNAfoldLauncher.overhanglimit[0] || pair.getOverHang() > RNAfoldLauncher.overhanglimit[1]) continue;
						//if(mergedSeqs.contains(pair.getSeq())) continue;	
						
						classifier.setClassification(pair);
						if(pair.getClassification().equals("M")) outM.println(pair);
						//else outU.println(pair);
						//mergedSeqs.add(pair.getSeq());
					}	
				}
			//	outU.close();
				outM.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
	}
	
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
	

	static void setMatchedNumPositions(String inFile, String outFile){
		ScoredPairOutputParser parser = new ScoredPairOutputParser(inFile);
		try{
			PrintStream out = new PrintStream(outFile);
			out.println(ScoredPair.getHeader());
			Hashtable<String, ArrayList<ScoredPair>> threePMap = new Hashtable<String, ArrayList<ScoredPair>>();
			for(ScoredPair pair : parser.getPairs()){
				ScoredPosition tp = pair.getPairedScoredPositions()[0];
			//	ScoredPosition fp = pair.getPairedScoredPositions()[1];
				String key = tp.getContig()+tp.getThreePPosition()+tp.isPlusStrand();
				if(!threePMap.containsKey(key)) threePMap.put(key, new ArrayList<ScoredPair>());
				threePMap.get(key).add(pair);
				//if(!fivePMap.containsKey(fp)) fivePMap.put(fp, new ArrayList<ScoredPair>());
				//f//ivePMap.get(fp).add(pair);
			}
			
			for(ArrayList<ScoredPair> pairs : threePMap.values()){
				int n = pairs.size();
				for(ScoredPair pair : pairs){
					pair.setMatched3pNum(n);
				}
			}
			
			threePMap = null;
			
			Hashtable<String, ArrayList<ScoredPair>> fivePMap = new Hashtable<String, ArrayList<ScoredPair>>();
			for(ScoredPair pair : parser.getPairs()){
				//ScoredPosition tp = pair.getPairedScoredPositions()[0];
				ScoredPosition fp = pair.getPairedScoredPositions()[1];
				//if(!threePMap.containsKey(tp)) threePMap.put(tp, new ArrayList<ScoredPair>());
				//threePMap.get(tp).add(pair);
				String key = fp.getContig()+fp.getThreePPosition()+fp.isPlusStrand();
				if(!fivePMap.containsKey(key)) fivePMap.put(key, new ArrayList<ScoredPair>());
				fivePMap.get(key).add(pair);
			}
			
			for(ArrayList<ScoredPair> pairs : fivePMap.values()){
				int n = pairs.size();
				for(ScoredPair pair : pairs){
					pair.setMatched5pNum(n);
				}
			}
			
			for(ArrayList<ScoredPair> pairs : fivePMap.values()){
				for(ScoredPair pair : pairs){
					out.println(pair);
				}
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	static void getScoredPairs(ScoredPositionOutputParser sparser, Classifier classifier, int blatHitThreshold, int seqLength, boolean sameDirection, 
			String outMfileName, int readThreshold, int numberThreads){
		try{
			ArrayList<ScoredPosition> allPositions = new ArrayList<ScoredPosition>();
			
			HashSet<String> seqs = new HashSet<String>();
			
			//int three = 0, five =0;
			ArrayList<Double> tscores = new ArrayList<Double>();
			ArrayList<Double> fscores = new ArrayList<Double>();
			//
			for(ScoredPosition sp : sparser.getPositions()){
			//	if(allPositions.size() > 5) break; //TODO
				if(sp.getClassification().toUpperCase().equals("M")) continue;
				if(sp.isPaired()) continue;
				if(sp.hasMatchingMiRNA()) continue;
				if(sp.getBlatHits() > blatHitThreshold) continue;
				if(sp.getSeq() == null || sp.getSeq().length() < FCLIP_Scorer.getFlankingNTNumber()) continue;
				if(seqs.contains(sp.getSeq())) continue;
			
				if(sp.is3pScored()){
					if(sp.getThreePReads()[0] < readThreshold || sp.getThreePReads()[1] < readThreshold) continue;// TODO
				}
				if(sp.is5pScored()){
					if(sp.getFivePReads()[0] < readThreshold || sp.getFivePReads()[1] < readThreshold) continue;// TODO
				}
				
				allPositions.add(sp);
				seqs.add(sp.getSeq());
				
				if(sp.is3pScored()){
			//		three++;
					tscores.add(sp.getThreePScore());
				}
				if(sp.is5pScored()){
			//		five++;
					fscores.add(sp.getFivePScore());
				}
			}		
			
			Collections.sort(tscores);
			Collections.sort(fscores);
			
		/*	if(positionNumber > 0){
				ArrayList<ScoredPosition> filteredPositions = new ArrayList<ScoredPosition>();
				
				double threshold3p = tscores.size() > positionNumber? tscores.get(tscores.size() - positionNumber - 1) : 0;
				double threshold5p = fscores.size() > positionNumber? fscores.get(fscores.size() - positionNumber - 1) : 0;
				
				for(ScoredPosition sp : allPositions){
					if(sp.is3pScored()){
						if(sp.getThreePScore() >= threshold3p){
							filteredPositions.add(sp);
							continue;
						}
					}
					if(sp.is5pScored()){
						if(sp.getFivePScore() >= threshold5p){
							filteredPositions.add(sp);
							continue;
						}
					}
				}
			
				allPositions = filteredPositions;
			}*/
			
			//
			/*if(tscores.size() > fscores.size()){
				double threshold = tscores.get(tscores.size() - fscores.size() - 1);
				for(ScoredPosition sp : allPositions){
					if(sp.is3pScored() && sp.getThreePScore() <= threshold) continue;
					filteredPositions.add(sp);
				}
			}else{
				double threshold = fscores.get(fscores.size() - tscores.size() - 1);
				for(ScoredPosition sp : allPositions){
					if(sp.is5pScored() && sp.getFivePScore() <= threshold) continue;
					filteredPositions.add(sp);
				}
			}*/
			
			System.out.println(sparser.getPositions().size() + " " + allPositions.size());
			
			ArrayList<Thread> threads = new ArrayList<Thread>();
			//ArrayList<String> tmpU = new ArrayList<String>();
			ArrayList<String> tmpM = new ArrayList<String>();
			
			for(int nt = 0; nt<numberThreads; nt++){
				ArrayList<ScoredPosition> compairngPositions = new ArrayList<ScoredPositionOutputParser.ScoredPosition>();
				// form compared / numberThreads /
				for(int j=0;j<allPositions.size();j++){
					if(j%numberThreads == nt) compairngPositions.add(allPositions.get(j));
				}
				ScorePairRunner runner = new ScorePairRunner(outMfileName + nt, 
						compairngPositions, allPositions, seqLength, classifier, sameDirection);
				//tmpU.add(outUfileName + nt);
				tmpM.add(outMfileName + nt);
				
				Thread thread = new Thread(runner, "Batch " + nt); //Thread created       
				thread.start();
				threads.add(thread);
			}			
			
			for(Thread rt : threads){
				try {
					rt.join();
				} catch (InterruptedException e) {						
					e.printStackTrace();
				}
			}
			
			PrintStream outM = new PrintStream(outMfileName);
		//	PrintStream outU = new PrintStream(outUfileName);
			
			outM.println(ScoredPair.getHeader());
			//outU.println(ScoredPair.getHeader());
			
			HashSet<String> mergedSeqs = new HashSet<String>();
			for(String tm : tmpM){
				for(ScoredPair pair : new ScoredPairOutputParser(tm).getPairs()){
					if(mergedSeqs.contains(pair.getSeq())) continue;
					outM.println(pair);
					mergedSeqs.add(pair.getSeq());
				}
				new File(tm).delete();
			}
			outM.close();
			
			
		/*	 mergedSeqs = new HashSet<String>();
			for(String tu : tmpU){
				for(ScoredPair pair : new ScoredPairOutputParser(tu).getPairs()){
					if(mergedSeqs.contains(pair.getSeq())) continue;
					outU.println(pair);
					mergedSeqs.add(pair.getSeq());
				}
				new File(tu).delete();				
			}
			outU.close();*/
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	static private void RandomizeSequences(ArrayList<ScoredPosition> allPositions, double energyThreshold, int depthThreshold, int hairpinThreshold){
		int m = 0;
		for(ScoredPosition sp : allPositions){
			boolean isSet = false;
			double originalEnergy = sp.getEnergy();
			int cntr = 0;
			while(!isSet){
				if(cntr++ > 50000) break;
				String rs = permutate(sp.getSeq());
				rs = randomSubstitute(rs, 8);
				RNAfoldLauncher rnafold = new RNAfoldLauncher(rs, FCLIP_Scorer.getFlankingNTNumber());
				if(rnafold.getNumberOfHairPins() > hairpinThreshold) continue;
				if(Math.abs(rnafold.getEnergyPerNT() - originalEnergy) > 0.1) continue;
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
	

